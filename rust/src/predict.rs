//! Pocket prediction from scored SAS points (PocketPredictor, SingleLinkageClustering,
//! score transformers) and residue score labeling (ModelBasedResidueLabeler).

use crate::atoms::{PointSet, P3};
use crate::features::pow2;
use crate::json::Json;
use crate::pdb::ResNum;
use crate::protein::Protein;
use std::collections::HashMap;

pub const PRED_POINT_THRESHOLD: f64 = 0.35;
pub const PRED_CLUSTERING_DIST: f64 = 3.0;
pub const PRED_MIN_CLUSTER_SIZE: usize = 3;
pub const EXTENDED_POCKET_CUTOFF: f64 = 3.5;
pub const PRED_PROTEIN_SURFACE_CUTOFF: f64 = 3.5;

pub struct LabeledPoint {
    pub pos: P3,
    pub score: f64,
    pub transformed: f64,
    pub predicted: bool,
    pub pocket: i32,
}

/// PointScoreCalculator.transformScore with point_score_pow = 2
pub fn transform_point_score(score: f64) -> f64 {
    let s = if score < 0.0 { 0.0 } else { score };
    pow2(s)
}

#[allow(dead_code)] // mirrors the Java data model
pub struct Pocket {
    pub name: String,
    pub rank: i32,
    pub new_rank: i32,
    pub score: f64,
    pub centroid: P3,
    /// extended pocket points (indices into labeled points)
    pub sas_points: Vec<usize>,
    pub cluster_count: usize,
    /// pocket surface atoms (atom indices, in P2Rank's order)
    pub surface_atoms: Vec<usize>,
    pub zscore: f64,
    pub proba: f64,
}

// ------------------------------------------------------------------ score transformers

pub enum ScoreTransformer {
    Zscore { mean: f64, stdev: f64 },
    Probability { min: f64, max: f64, nbins: i32, tp: Vec<i32>, fp: Vec<i32> },
}

impl ScoreTransformer {
    pub fn load(path: &std::path::Path) -> Result<ScoreTransformer, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {}: {}", path.display(), e))?;
        let j = crate::json::parse(&text)?;
        let name = j.get("name").and_then(|n| n.as_str()).ok_or("missing transformer name")?;
        let p = j.get("params").ok_or("missing transformer params")?;
        let num = |k: &str| -> f64 { p.get(k).and_then(Json::as_f64).unwrap_or(0.0) };
        let ints = |k: &str| -> Vec<i32> {
            p.get(k)
                .and_then(Json::as_array)
                .map(|a| a.iter().map(|v| v.as_f64().unwrap_or(0.0) as i32).collect())
                .unwrap_or_default()
        };
        match name {
            "ZscoreTpTransformer" => Ok(ScoreTransformer::Zscore { mean: num("mean"), stdev: num("stdev") }),
            "ProbabilityScoreTransformer" => Ok(ScoreTransformer::Probability {
                min: num("min"),
                max: num("max"),
                nbins: num("nbins") as i32,
                tp: ints("tp_cumul_hist"),
                fp: ints("fp_cumul_hist"),
            }),
            _ => Err(format!("Invalid score transformer name: {}", name)),
        }
    }

    pub fn transform(&self, raw: f64) -> f64 {
        match self {
            ScoreTransformer::Zscore { mean, stdev } => (raw - mean) / stdev,
            ScoreTransformer::Probability { min, max, nbins, tp, fp } => {
                let interp = |hist: &Vec<i32>, x: f64| -> f64 {
                    let step = (max - min) / (*nbins as f64);
                    let mut idx = ((x - min) / step) as i32;
                    if idx < 0 {
                        idx = 0;
                    }
                    if idx >= *nbins {
                        idx = nbins - 1;
                    }
                    if idx == nbins - 1 {
                        return hist[idx as usize] as f64;
                    }
                    let m = x - min - step * (idx as f64);
                    let add = (hist[idx as usize + 1] - hist[idx as usize]) as f64 * (m / step);
                    hist[idx as usize] as f64 + add
                };
                let tpx = interp(tp, raw);
                let fpx = interp(fp, raw);
                let mut res = tpx / (tpx + fpx);
                if res <= 0.0 {
                    res = 0.0;
                }
                res
            }
        }
    }
}

pub struct Transformers {
    pub zscore: Option<ScoreTransformer>,
    pub proba: Option<ScoreTransformer>,
    pub res_zscore: Option<ScoreTransformer>,
    pub res_proba: Option<ScoreTransformer>,
}

// ------------------------------------------------------------------ clustering

/// SingleLinkageClustering: connected components of ligandable points (sqr dist <= 9),
/// clusters ordered by first member, members in point order; then filtered by size.
fn cluster_points(points: &[LabeledPoint]) -> Vec<Vec<usize>> {
    let ligandable: Vec<usize> = (0..points.len()).filter(|&i| points[i].predicted).collect();
    let n = ligandable.len();
    if n == 0 {
        return Vec::new();
    }
    let set = PointSet::new(ligandable.iter().map(|&i| points[i].pos).collect(), (0..n as u32).collect());
    let min_dist = PRED_CLUSTERING_DIST * PRED_CLUSTERING_DIST;
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let kd = set.kd();
    let mut found = Vec::new();
    for i in 0..n {
        found.clear();
        let p = set.pts[i];
        kd.for_each_within(p[0], p[1], p[2], min_dist, |j| found.push(j as usize));
        for &j in &found {
            if j <= i {
                continue;
            }
            // exact predicate of SLinkClustererV2: Struct.sqrDist(a, b) <= minDist
            let q = set.pts[j];
            let dx = p[0] - q[0];
            let dy = p[1] - q[1];
            let dz = p[2] - q[2];
            if dx * dx + dy * dy + dz * dz <= min_dist {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }
    let mut order: Vec<usize> = Vec::new();
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        members.entry(r).or_insert_with(|| {
            order.push(r);
            Vec::new()
        });
        members.get_mut(&r).unwrap().push(ligandable[i]);
    }
    order.into_iter().map(|r| members.remove(&r).unwrap()).filter(|c| c.len() >= PRED_MIN_CLUSTER_SIZE).collect()
}

/// PocketPredictor.predictPockets
pub fn predict_pockets(points: &mut [LabeledPoint], all_points: &PointSet, exposed: &PointSet, tr: &Transformers) -> Vec<Pocket> {
    let clusters = cluster_points(points);
    let mut pockets: Vec<Pocket> = Vec::new();
    for cluster in clusters {
        let cluster_set = PointSet::new(cluster.iter().map(|&i| points[i].pos).collect(), cluster.iter().map(|&i| i as u32).collect());
        let pocket_points = all_points.cutout_shell(&cluster_set, EXTENDED_POCKET_CUTOFF);
        let surface = exposed.cutout_shell(&pocket_points, PRED_PROTEIN_SURFACE_CUTOFF);

        let mut point_scores: Vec<f64> = pocket_points.ids.iter().map(|&i| points[i as usize].transformed).collect();
        // Arrays.sort(order, Double.compare(b, a)): descending, stable
        point_scores.sort_by(|a, b| b.total_cmp(a));
        let mut score = 0.0f64;
        for s in &point_scores {
            score += *s;
        }

        let centroid = cluster_set.centroid();
        let zscore = tr.zscore.as_ref().map(|t| t.transform(score)).unwrap_or(0.0);
        let proba = tr.proba.as_ref().map(|t| t.transform(score)).unwrap_or(0.0);
        pockets.push(Pocket {
            name: "pocket".to_string(),
            rank: 0,
            new_rank: 0,
            score,
            centroid,
            sas_points: pocket_points.ids.iter().map(|&i| i as usize).collect(),
            cluster_count: cluster.len(),
            surface_atoms: surface.ids.iter().map(|&a| a as usize).collect(),
            zscore,
            proba,
        });
    }
    // descending by score, stable (Double.compare semantics)
    pockets.sort_by(|a, b| java_double_compare(b.score, a.score));

    for (i, p) in pockets.iter_mut().enumerate() {
        let r = (i + 1) as i32;
        p.rank = r;
        p.name = format!("pocket{}", r);
        p.new_rank = r;
    }
    // finalize: assign lp.pocket, best pocket wins for shared points
    for p in &pockets {
        for &i in &p.sas_points {
            points[i].pocket = 0;
        }
    }
    for p in &pockets {
        for &i in &p.sas_points {
            if points[i].pocket == 0 {
                points[i].pocket = p.new_rank;
            }
        }
    }
    pockets
}

/// java.lang.Double.compare
pub fn java_double_compare(a: f64, b: f64) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    if a < b {
        return Less;
    }
    if a > b {
        return Greater;
    }
    let ab = if a.is_nan() { 0x7ff8000000000000u64 as i64 } else { a.to_bits() as i64 };
    let bb = if b.is_nan() { 0x7ff8000000000000u64 as i64 } else { b.to_bits() as i64 };
    ab.cmp(&bb)
}

// ------------------------------------------------------------------ residues

pub struct ResidueRow {
    pub score: f64,
    pub zscore: f64,
    pub proba: f64,
    pub pocket: i32,
}

/// ResidueLabelings.calculate: per residue score, zscore, probability and pocket reference
pub fn label_residues(p: &Protein, points: &[LabeledPoint], all_points: &PointSet, pockets: &[Pocket], tr: &Transformers) -> Vec<ResidueRow> {
    let radius = (1.6f64 + 2.5f64) + 2.0f64;
    let mut scores: Vec<f64> = Vec::with_capacity(p.residues.len());
    for r in &p.residues {
        let atoms = p.residue_heavy_atoms(r);
        let res_set = p.point_set(&atoms);
        let shell = all_points.cutout_shell(&res_set, radius);
        let mut sum = 0.0f64;
        for &i in &shell.ids {
            sum += transform_point_score(points[i as usize].score);
        }
        // base = pow(n, residue_score_sum_to_avg = 0) = 1
        let score = if shell.ids.is_empty() { 0.0 } else { sum / 1.0 };
        scores.push(score);
    }

    // pocket reference labeling, keyed by residue key (last pocket in reversed order wins)
    let mut pocket_labels: HashMap<&ResNum, i32> = HashMap::new();
    for r in &p.residues {
        pocket_labels.insert(p.residue_key(r), 0);
    }
    for pk in pockets.iter().rev() {
        for &a in &pk.surface_atoms {
            let g = p.s.atoms[a].group;
            pocket_labels.insert(&p.s.groups[g].resnum, pk.rank);
        }
    }

    // labelings are maps keyed by residue key: duplicate keys resolve to the last residue
    let mut score_by_key: HashMap<&ResNum, f64> = HashMap::new();
    for (i, r) in p.residues.iter().enumerate() {
        score_by_key.insert(p.residue_key(r), scores[i]);
    }

    p.residues
        .iter()
        .map(|r| {
            let key = p.residue_key(r);
            let score = score_by_key[key];
            let zscore = tr.res_zscore.as_ref().map(|t| t.transform(score)).unwrap_or(f64::NAN);
            let proba = tr.res_proba.as_ref().map(|t| t.transform(score)).unwrap_or(f64::NAN);
            ResidueRow { score, zscore, proba, pocket: pocket_labels[key] }
        })
        .collect()
}
