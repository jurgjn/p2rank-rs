//! The `predict` pipeline for a single structure.

use crate::atoms::PointSet;
use crate::elements::surface_vdw_radius;
use crate::features::{FeatureContext, FVec};
use crate::forest::Forest;
use crate::output;
use crate::predict::{self, LabeledPoint, Transformers};
use crate::protein::Protein;
use std::path::{Path, PathBuf};

pub struct Config {
    pub install_dir: PathBuf,
    pub threads: usize,
    pub visualizations: bool,
    /// use the full CCD (BioJava online behaviour) instead of BioJava's reduced fallback set
    pub full_ccd: bool,
    /// random forest early-exit stages (empty = exact)
    pub rf_stages: Vec<(usize, f64)>,
    /// fast mode (approximations allowed)
    pub fast: bool,
    /// faster mode (implies fast): score a sparse subset of SAS points, fill in near pockets
    pub faster: Option<FasterParams>,
    /// also write per-pocket descriptors (`<label>_pockets_descriptors.csv`)
    pub descriptors: bool,
}

/// `-mode faster` parameters (see `faster_scores`)
#[derive(Clone, Debug)]
pub struct FasterParams {
    /// seed only points on the 12 level-1 directions (of the 42 level-2 ones)
    pub seed_coarse: bool,
    /// seed only points with at least this protrusion (deep-layer atoms within 10 A); 0 = off
    pub prot_gate: f64,
    /// points scoring at least this trigger refinement of their neighborhood
    pub refine_score: f64,
    /// seed points scoring at least this also trigger refinement (<= refine_score)
    pub seed_trigger: f64,
    /// refinement neighborhood radius
    pub refine_radius: f64,
    /// early-exit stages used when refining (empty = full forest)
    pub refine_stages: Vec<(usize, f64)>,
    /// use only the first `trees` trees of the forest
    pub trees: usize,
}

pub struct Log {
    lines: std::sync::Mutex<Vec<String>>,
    pub quiet: bool,
}

impl Log {
    pub fn new(quiet: bool) -> Log {
        Log { lines: std::sync::Mutex::new(Vec::new()), quiet }
    }
    pub fn info(&self, src: &str, msg: &str) {
        let line = format!("[INFO] {} - {}", src, msg);
        if !self.quiet {
            println!("{}", line);
        }
        self.lines.lock().unwrap().push(line);
    }
    pub fn error(&self, src: &str, msg: &str) {
        let line = format!("[ERROR] {} - {}", src, msg);
        eprintln!("{}", line);
        self.lines.lock().unwrap().push(line);
    }
    pub fn take(&self) -> Vec<String> {
        std::mem::take(&mut self.lines.lock().unwrap())
    }
}

/// Futils.inputStream: transparently decompress .gz files
pub fn read_structure_bytes(path: &Path) -> Result<Vec<u8>, ItemError> {
    let raw = std::fs::read(path).map_err(|e| file_open_error(path, &e))?;
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if name.ends_with(".gz") {
        use std::io::Read;
        let mut d = flate2::read::MultiGzDecoder::new(&raw[..]);
        let mut out = Vec::new();
        d.read_to_end(&mut out).map_err(|e| ItemError::new("ZipException", e.to_string()))?;
        Ok(out)
    } else if name.ends_with(".zst") || name.ends_with(".zstd") {
        use std::io::Read;
        let mut d = ruzstd::decoding::StreamingDecoder::new(&raw[..]).map_err(|e| ItemError::new("IOException", format!("{:?}", e)))?;
        let mut out = Vec::new();
        d.read_to_end(&mut out).map_err(|e| ItemError::new("IOException", e.to_string()))?;
        Ok(out)
    } else if name.ends_with(".lzma") {
        Err(ItemError::new("UnsupportedOperationException", format!("lzma input is not supported by prank-rs: {}", name)))
    } else {
        Ok(raw)
    }
}

/// Futils.realExtension
fn real_extension(name: &str) -> String {
    let mut n = name.to_string();
    for ext in ["gz", "zst", "zip", "bz2"] {
        if let Some(stripped) = n.strip_suffix(&format!(".{}", ext)) {
            n = stripped.to_string();
            break;
        }
    }
    n.rsplit_once('.').map(|(_, e)| e.to_string()).unwrap_or_default()
}

fn parallel_map<T: Send, R: Send + Default + Clone>(items: &[T], threads: usize, f: impl Fn(&T, &mut Vec<u32>) -> R + Sync) -> Vec<R>
where
    T: Sync,
{
    let n = items.len();
    let mut out: Vec<R> = vec![R::default(); n];
    // at least ~512 items per thread: thread start-up is not free
    let threads = threads.max(1).min(n.div_ceil(512).max(1));
    if threads <= 1 {
        let mut scratch = Vec::new();
        for (i, it) in items.iter().enumerate() {
            out[i] = f(it, &mut scratch);
        }
        return out;
    }
    let chunk = n.div_ceil(threads);
    std::thread::scope(|s| {
        for (items_c, out_c) in items.chunks(chunk).zip(out.chunks_mut(chunk)) {
            let f = &f;
            s.spawn(move || {
                let mut scratch = Vec::new();
                for (it, o) in items_c.iter().zip(out_c.iter_mut()) {
                    *o = f(it, &mut scratch);
                }
            });
        }
    });
    out
}

fn parallel_map_blocks<T: Sync, R: Send>(items: &[T], threads: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let n = items.len();
    let threads = threads.min(n);
    if threads <= 1 || n <= 1 {
        return items.iter().map(&f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<R>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.min(n) {
            s.spawn(|| loop {
                let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if k >= n {
                    break;
                }
                let r = f(&items[k]);
                *slots[k].lock().unwrap() = Some(r);
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
}

/// Fast mode, after the early-exit cascade: early-stopped points near (likely) ligandable points
/// get the full forest score, repeated until stable, so pocket regions are scored exactly. Only
/// background points far from any pocket keep their partial (early-exit) score.
fn refine_fast_scores(sas: &[[f64; 3]], vectors: &[V], sc: &mut [f64], used: &mut [u16], forest: &Forest, cfg: &Config) {
    let full = forest.num_trees as u16;
    loop {
        let hi: Vec<[f64; 3]> = (0..sas.len()).filter(|&i| sc[i] >= FAST_REFINE_SCORE).map(|i| sas[i]).collect();
        let grid = match crate::atoms::CountGrid::new(&hi, FAST_REFINE_RADIUS) {
            Some(g) => g,
            None => break,
        };
        let todo: Vec<usize> = (0..sas.len()).filter(|&i| used[i] < full && grid.any_within(&sas[i], FAST_REFINE_RADIUS)).collect();
        if todo.is_empty() {
            break;
        }
        let blocks: Vec<&[usize]> = todo.chunks(256).collect();
        let res: Vec<Vec<f64>> = parallel_map_blocks(&blocks, cfg.threads, |b| {
            let xs: Vec<&[f64]> = b.iter().map(|&i| &vectors[i].0[..]).collect();
            let mut out = vec![0.0f64; b.len()];
            forest.predict_block_cascade(&xs, &[], &mut out);
            out
        });
        for (b, r) in blocks.iter().zip(res) {
            for (&i, v) in b.iter().zip(r) {
                sc[i] = v;
                used[i] = full;
            }
        }
    }
}

/// Faster mode. Feature vectors and random forest scores are computed lazily, for three sets of
/// points, keeping per-point class sums so that more trees can be added later:
/// 1. seeds: points on the coarse (level-1) tessellation directions (optionally only those with
///    enough protrusion), scored with the early-exit cascade;
/// 2. expansion: every point within `refine_radius` of a trigger point (score >= `refine_score`,
///    or a seed with score >= `seed_trigger`) is scored with the `refine_stages` cascade;
/// 3. completion: every point within the extended pocket cutoff (3.5 A) of a point scoring
///    >= 0.35 gets all `trees` trees, so pockets (clustering and pocket scores) only see full scores.
/// Steps 2 and 3 repeat until stable. Early-exited points with a (near) zero partial score get
/// the typical full score of such points, and points never scored get a distance-weighted mean of
/// nearby scored points (capped below the refine trigger), or that same typical score.
fn faster_scores(sas: &[[f64; 3]], tags: &[u32], ctx: &FeatureContext, forest: &Forest, cfg: &Config, fp: &FasterParams) -> Vec<f64> {
    let n = sas.len();
    let full = fp.trees.clamp(1, forest.num_trees) as u16;
    let coarse = crate::surface::nested_direction_mask(1, 2).expect("tessellation tables");
    let mut is_seed: Vec<bool> = (0..n).map(|i| !fp.seed_coarse || coarse[crate::surface::tag_dir(tags[i])]).collect();
    if fp.prot_gate > 0.0 {
        let idx: Vec<u32> = (0..n as u32).filter(|&i| is_seed[i as usize]).collect();
        let keep: Vec<bool> = parallel_map(&idx, cfg.threads, |&i, _| ctx.protrusion(&sas[i as usize]) >= fp.prot_gate);
        for (&i, k) in idx.iter().zip(keep) {
            is_seed[i as usize] = k;
        }
    }

    let mut vectors: Vec<V> = vec![V::default(); n];
    let mut has_vec = vec![false; n];
    let mut s0 = vec![0.0f64; n];
    let mut s1 = vec![0.0f64; n];
    // trees in the class sums (0 = not scored)
    let mut used = vec![0u16; n];
    let mut sc = vec![0.0f64; n];

    // scores each job's points up to `full` trees (resuming from the current sums) with the job's
    // early-exit stages, computing missing feature vectors on the way (one parallel pass)
    let run = |jobs: &mut [(&mut Vec<u32>, &[(usize, f64)])], vectors: &mut Vec<V>, has_vec: &mut Vec<bool>, s0: &mut Vec<f64>, s1: &mut Vec<f64>, used: &mut Vec<u16>, sc: &mut Vec<f64>| {
        let total: usize = jobs.iter().map(|j| j.0.len()).sum();
        // blocks of points with the same number of trees done; small blocks if there are few
        // points, so that all threads get work (scores do not depend on the blocking)
        let block = (total / (2 * cfg.threads)).clamp(32, 256);
        let mut blocks: Vec<(&[u32], &[(usize, f64)])> = Vec::new();
        for (todo, stages) in jobs.iter_mut() {
            todo.sort_by_key(|&i| (used[i as usize], i));
            let todo: &[u32] = todo;
            let mut k = 0;
            while k < todo.len() {
                let u = used[todo[k] as usize];
                let mut e = k;
                while e < todo.len() && e - k < block && used[todo[e] as usize] == u {
                    e += 1;
                }
                blocks.push((&todo[k..e], *stages));
                k = e;
            }
        }
        let (vectors_r, has_vec_r, s0r, s1r, usedr) = (&*vectors, &*has_vec, &*s0, &*s1, &*used);
        let res: Vec<(Vec<(u32, V)>, Vec<f64>, Vec<f64>, Vec<u16>)> = parallel_map_blocks(&blocks, cfg.threads, |(b, stages)| {
            let mut scratch = Vec::new();
            let new: Vec<(u32, V)> = b.iter().filter(|&&i| !has_vec_r[i as usize]).map(|&i| (i, V(ctx.sas_vector(&sas[i as usize], &mut scratch)))).collect();
            let mut k = 0;
            let xs: Vec<&[f64]> = b
                .iter()
                .map(|&i| {
                    if has_vec_r[i as usize] {
                        &vectors_r[i as usize].0[..]
                    } else {
                        k += 1;
                        &new[k - 1].1 .0[..]
                    }
                })
                .collect();
            let mut a: Vec<f64> = b.iter().map(|&i| s0r[i as usize]).collect();
            let mut c: Vec<f64> = b.iter().map(|&i| s1r[i as usize]).collect();
            let u = forest.advance(&xs, &mut a, &mut c, usedr[b[0] as usize] as usize, stages, full as usize);
            (new, a, c, u)
        });
        for ((b, _), (new, a, c, u)) in blocks.iter().zip(res) {
            for (i, v) in new {
                vectors[i as usize] = v;
                has_vec[i as usize] = true;
            }
            for (j, &i) in b.iter().enumerate() {
                let i = i as usize;
                s0[i] = a[j];
                s1[i] = c[j];
                used[i] = u[j];
                sc[i] = Forest::normalized(a[j], c[j]);
            }
        }
    };

    let mut seeds: Vec<u32> = (0..n as u32).filter(|&i| is_seed[i as usize]).collect();
    let n_seeds = seeds.len();
    run(&mut [(&mut seeds, &cfg.rf_stages)], &mut vectors, &mut has_vec, &mut s0, &mut s1, &mut used, &mut sc);

    // all points, for neighborhood queries around trigger points
    let all = crate::atoms::CountGrid::new_indexed(sas, fp.refine_radius.max(predict::EXTENDED_POCKET_CUTOFF));
    let mut expanded = vec![false; n]; // neighborhood expanded (as a trigger)
    let mut completed = vec![false; n]; // neighborhood completed (as a pocket point)
    let mut in_expansion = vec![false; n];
    loop {
        let mut grow: Vec<u32> = Vec::new();
        let mut complete: Vec<u32> = Vec::new();
        let mut mark = vec![false; n];
        if let Some(all) = &all {
            for i in 0..n {
                if used[i] == 0 {
                    continue;
                }
                if !completed[i] && sc[i] >= predict::PRED_POINT_THRESHOLD {
                    completed[i] = true;
                    all.for_each_within(&sas[i], predict::EXTENDED_POCKET_CUTOFF, |j| {
                        let j = j as usize;
                        if used[j] < full && !mark[j] {
                            mark[j] = true;
                            complete.push(j as u32);
                        }
                    });
                }
                if !expanded[i] && (sc[i] >= fp.refine_score || (is_seed[i] && sc[i] >= fp.seed_trigger)) {
                    expanded[i] = true;
                    all.for_each_within(&sas[i], fp.refine_radius, |j| {
                        let j = j as usize;
                        if !in_expansion[j] && used[j] < full {
                            in_expansion[j] = true;
                            grow.push(j as u32);
                        }
                    });
                }
            }
        }
        grow.retain(|&j| !mark[j as usize]);
        if grow.is_empty() && complete.is_empty() {
            break;
        }
        run(&mut [(&mut complete, &[]), (&mut grow, &fp.refine_stages)], &mut vectors, &mut has_vec, &mut s0, &mut s1, &mut used, &mut sc);
    }

    // early-exited points mostly have a partial score of 0 after 10 trees: use their typical
    // full-forest score instead, which keeps residue scores (sums of squared scores) unbiased
    for i in 0..n {
        if used[i] > 0 && used[i] < full && sc[i] < FASTER_EARLY_EXIT_SCORE {
            sc[i] = FASTER_EARLY_EXIT_SCORE;
        }
    }
    // impute never-scored points: inverse-distance weighted mean of scored points within
    // FASTER_IMPUTE_RADIUS (capped below the refine trigger), or the early-exit score
    let scored: Vec<u32> = (0..n as u32).filter(|&i| used[i as usize] > 0).collect();
    let scored_pos: Vec<[f64; 3]> = scored.iter().map(|&i| sas[i as usize]).collect();
    let grid = crate::atoms::CountGrid::new_indexed(&scored_pos, FASTER_IMPUTE_RADIUS);
    let unscored: Vec<u32> = (0..n as u32).filter(|&i| used[i as usize] == 0).collect();
    let cap = fp.refine_score - 0.01;
    let fill: Vec<f64> = parallel_map(&unscored, cfg.threads, |&i, _| {
        let q = &sas[i as usize];
        let (mut ws, mut wv) = (0.0, 0.0);
        if let Some(grid) = &grid {
            grid.for_each_within(q, FASTER_IMPUTE_RADIUS, |k| {
                let j = scored[k as usize] as usize;
                let w = 1.0 / (crate::atoms::dist(q, &sas[j]) + 0.5);
                ws += w;
                wv += w * sc[j];
            });
        }
        if ws > 0.0 { (wv / ws).min(cap) } else { FASTER_EARLY_EXIT_SCORE }
    });
    for (&i, v) in unscored.iter().zip(fill) {
        sc[i as usize] = v;
    }
    if std::env::var("P2RS_DEV_STATS").is_ok() {
        let nf = has_vec.iter().filter(|&&b| b).count();
        let nfull = used.iter().filter(|&&u| u == full).count();
        eprintln!("faster: points {} seeds {} vectors {} full-forest {}", n, n_seeds, nf, nfull);
    }
    sc
}

/// RMS full-forest score of points stopped after 10 trees (0.0297 on ah4h, which is neither
/// training nor evaluation data; on the training set chen11 it is lower, 0.0179)
const FASTER_EARLY_EXIT_SCORE: f64 = 0.03;

/// neighborhood radius for imputing scores of never-scored points (coarse points are ~3.5 A apart)
const FASTER_IMPUTE_RADIUS: f64 = 3.0;

/// points with at least this (partial or full) score trigger exact rescoring of their neighbors
const FAST_REFINE_SCORE: f64 = 0.2;
/// neighborhood radius for exact rescoring (covers pocket extension and clustering distances)
const FAST_REFINE_RADIUS: f64 = 4.0;

#[derive(Clone, Copy)]
struct V(FVec);
impl Default for V {
    fn default() -> Self {
        V([0.0; crate::features::NUM_FEATURES])
    }
}

pub struct Item {
    pub path: PathBuf,
    /// reduce the structure to these chains (author ids), None = all
    pub chains: Option<Vec<String>>,
}

#[allow(dead_code)] // mirrors the Java data model
pub struct ItemResult {
    pub label: String,
}

/// An item error, modeled after the Java exception cause chain
/// (list of (exception class simple name, message)) used in P2Rank's error reports.
#[derive(Debug)]
pub struct ItemError(pub Vec<(String, String)>);

impl ItemError {
    pub fn new(class: &str, msg: impl Into<String>) -> ItemError {
        ItemError(vec![(class.to_string(), msg.into())])
    }
    pub fn wrap(self, class: &str, msg: impl Into<String>) -> ItemError {
        let mut v = vec![(class.to_string(), msg.into())];
        v.extend(self.0);
        ItemError(v)
    }
    /// Dataset.formatErrorMessage
    pub fn format(&self) -> String {
        let s: Vec<String> = self.0.iter().map(|(c, m)| format!("{}: {}", c, m)).collect();
        s.join(" | ").replace('\n', " ").replace(',', " ")
    }
}

impl From<String> for ItemError {
    fn from(s: String) -> Self {
        ItemError::new("PrankException", s)
    }
}

fn file_open_error(path: &Path, e: &std::io::Error) -> ItemError {
    let reason = if path.is_dir() {
        "Is a directory".to_string()
    } else {
        match e.kind() {
            std::io::ErrorKind::NotFound => "No such file or directory".to_string(),
            std::io::ErrorKind::PermissionDenied => "Permission denied".to_string(),
            _ => e.to_string(),
        }
    };
    ItemError::new("FileNotFoundException", format!("{} ({})", path.display(), reason))
}

/// Runs prediction for one structure file and writes all output files into `outdir`.
pub fn predict_file(item: &Item, outdir: &Path, cfg: &Config, forest: &Forest, tr: &Transformers, log: &Log) -> Result<ItemResult, ItemError> {
    let timing = std::env::var("P2RS_TIMING").is_ok();
    let t0 = std::time::Instant::now();
    let tick = |what: &str| {
        if timing {
            eprintln!("timing: {:>8.1} ms  {}", t0.elapsed().as_secs_f64() * 1000.0, what);
        }
    };
    let file: &Path = &item.path;
    let label = file.file_name().and_then(|n| n.to_str()).ok_or_else(|| ItemError::new("IllegalArgumentException", "invalid file name"))?.to_string();
    let abs = std::fs::canonicalize(file).unwrap_or_else(|_| file.to_path_buf());
    log.info("Protein", &format!("loading protein [{}]", abs.display()));
    log.info("PdbUtils", &format!("loading file [{}]", file.display()));

    let ext = real_extension(&label);
    let mut bytes = read_structure_bytes(file)?;
    let dict = if cfg.full_ccd { crate::chemcomp::ChemCompDictionary::full() } else { crate::chemcomp::ChemCompDictionary::reduced() };
    let load_err = |class: String, msg: String| {
        ItemError::new(&class, msg).wrap("PrankException", format!("Failed to load structure from '{}'", file.display()))
    };
    let split_err = |e: &str| -> (String, String) {
        match e.split_once(": ") {
            Some((c, m)) if !c.contains(' ') => (c.to_string(), m.to_string()),
            _ => ("RuntimeException".to_string(), e.to_string()),
        }
    };
    let s = if ext == "cif" || ext == "bcif" {
        // CifIO guesses the format from the first byte: gzip, BinaryCIF (msgpack map) or text
        if bytes.first() == Some(&0x1f) {
            use std::io::Read;
            let mut out = Vec::new();
            flate2::read::MultiGzDecoder::new(&bytes[..])
                .read_to_end(&mut out)
                .map_err(|e| load_err("ZipException".into(), e.to_string()))?;
            bytes = out;
        }
        crate::cif::parse_cif(&bytes, dict).map_err(|e| {
            let (c, m) = split_err(&e.0);
            load_err(c, m)
        })?
    } else {
        crate::pdb::parse_pdb(&bytes, dict).map_err(|e| {
            let (c, m) = split_err(&e.0);
            load_err(c, m)
        })?
    };
    tick("parsed");
    let n_struct_atoms: usize = {
        let mut n = 0;
        if let Some(m) = s.models.first() {
            for &c in m {
                for &g in &s.chains[c].groups {
                    n += s.groups[g].atoms.len();
                }
            }
        }
        n
    };
    let p = Protein::from_structure(label.clone(), s, item.chains.as_deref())?;
    tick("protein");
    log.info("Protein", &format!("structure atoms: {}", n_struct_atoms));
    log.info("Protein", &format!("protein   atoms: {}", p.protein_atoms.len()));
    if p.protein_atoms.is_empty() {
        log.error("Protein", &format!("Structure with no protein chain atoms! [{}]", label));
    }
    log.info("Protein", "ignoring ligands");

    // surface
    let centers: Vec<[f64; 3]> = p.protein_atoms.iter().map(|&a| p.atom_pos(a)).collect();
    let mut vdw = Vec::with_capacity(centers.len());
    for &a in &p.protein_atoms {
        vdw.push(surface_vdw_radius(p.s.atoms[a].element, &p.s.atoms[a].name).map_err(|e| ItemError::new("IllegalArgumentException", e))?);
    }
    let (sas, sas_tags) = crate::surface::compute_surface_points_tagged(&centers, &vdw, 1.6, 2, cfg.threads).map_err(|e| ItemError::new("IllegalArgumentException", e))?;
    tick("surface");
    log.info("Protein", &format!("SAS points (solventRadius=1.6, tessellation=2): {}", sas.len()));

    let protein_set = p.point_set(&p.protein_atoms);
    let sas_set = PointSet::new(sas.clone(), (0..sas.len() as u32).collect());
    let exposed = crate::atoms::cutout_shell_par(&protein_set, &sas_set, 1.6 + 2.5, cfg.threads);
    log.info("Protein", &format!("exposed protein atoms: {} of {}", exposed.len(), p.protein_atoms.len()));
    let deep = crate::atoms::cutout_shell_par(&protein_set, &exposed, 10.0, cfg.threads);
    if deep.len() >= crate::atoms::USE_KDTREE_CUTOUT_SPHERE_THRESHOLD {
        deep.kd();
    }
    tick("exposed+deep");
    let ctx = FeatureContext::new(&p, exposed, deep, cfg.fast);

    if ctx.exposed.len() >= crate::atoms::USE_KDTREE_CUTOUT_SPHERE_THRESHOLD {
        ctx.exposed.kd();
    }

    // features + classification
    let scores: Vec<f64> = if let Some(fp) = &cfg.faster {
        faster_scores(&sas, &sas_tags, &ctx, forest, cfg, fp)
    } else {
        let vectors: Vec<V> = parallel_map(&sas, cfg.threads, |pt, scratch| V(ctx.sas_vector(pt, scratch)));
        tick("features");
        // classify in blocks (tree-major traversal keeps each tree cache resident)
        const BLOCK: usize = 256;
        let blocks: Vec<&[V]> = vectors.chunks(BLOCK).collect();
        let block_scores: Vec<(Vec<f64>, Vec<u16>)> = parallel_map_blocks(&blocks, cfg.threads, |b| {
            let xs: Vec<&[f64]> = b.iter().map(|v| &v.0[..]).collect();
            let mut out = vec![0.0f64; b.len()];
            let used = forest.predict_block_cascade(&xs, &cfg.rf_stages, &mut out);
            (out, used)
        });
        let mut sc: Vec<f64> = Vec::with_capacity(sas.len());
        let mut used: Vec<u16> = Vec::with_capacity(sas.len());
        for (o, u) in block_scores {
            sc.extend(o);
            used.extend(u);
        }
        if !cfg.rf_stages.is_empty() {
            refine_fast_scores(&sas, &vectors, &mut sc, &mut used, forest, cfg);
        }
        sc
    };
    tick("features+rf");
    let mut points: Vec<LabeledPoint> = sas
        .iter()
        .zip(scores.iter())
        .map(|(pt, &sc)| LabeledPoint {
            pos: *pt,
            score: sc,
            transformed: predict::transform_point_score(sc),
            predicted: sc >= predict::PRED_POINT_THRESHOLD,
            pocket: 0,
        })
        .collect();

    log.info("PocketPredictor", "PREDICTING POCKETS.... ====================================");
    log.info("PocketPredictor", &format!("SAS POINTS: {}", points.len()));
    let all_points = PointSet::new(sas.clone(), (0..sas.len() as u32).collect());
    all_points.kd();
    let pockets = predict::predict_pockets(&mut points, &all_points, &ctx.exposed, tr);
    for (i, pk) in pockets.iter().enumerate() {
        log.info(
            "PocketPredictor",
            &format!(
                "pocket{} -  surf_atoms: {}   points: {}   score: {}",
                crate::jfmt::format_wd((i + 1) as i64, 2),
                crate::jfmt::format_wd(pk.surface_atoms.len() as i64, 3),
                crate::jfmt::format_wd(pk.cluster_count as i64, 3),
                crate::jfmt::format_wf(pk.score, 6, 1)
            ),
        );
    }

    tick("pockets");
    // SAS points file (formatting + gzip) is written concurrently with residue labeling
    let vis = outdir.join("visualizations");
    let points_rel = format!("data/{}_points.pdb.gz", label);
    let points_job = |points: &[LabeledPoint]| -> Result<(), ItemError> {
        if !cfg.visualizations {
            return Ok(());
        }
        let data = vis.join("data");
        std::fs::create_dir_all(&data).map_err(|e| ItemError::new("IOException", format!("cannot create {}: {}", data.display(), e)))?;
        let pts_text = output::points_pdb(points, cfg.threads);
        output::write_gzip_java(&vis.join(&points_rel), pts_text.as_bytes(), if cfg.fast { 1 } else { 6 })
            .map_err(|e| ItemError::new("IOException", format!("cannot write points: {}", e)))
    };
    let (points_res, residue_rows) = std::thread::scope(|sc| {
        let h = sc.spawn(|| points_job(&points));
        let rows = predict::label_residues(&p, &points, &all_points, &pockets, tr, cfg.fast);
        (h.join().unwrap(), rows)
    });
    points_res?;
    tick("residues+points");

    // output files
    let w = |name: String, content: String| -> Result<(), ItemError> {
        let path = outdir.join(&name);
        std::fs::write(&path, content).map_err(|e| ItemError::new("IOException", format!("cannot write {}: {}", path.display(), e)))
    };
    w(format!("{}_predictions.csv", label), output::predictions_csv(&p, &pockets))?;
    w(format!("{}_residues.csv", label), output::residues_csv(&p, &residue_rows))?;
    if cfg.descriptors {
        w(format!("{}_pockets_descriptors.csv", label), crate::descriptors::csv(&p, &pockets, &points, &ctx))?;
    }

    if cfg.visualizations {
        let data = vis.join("data");
        let copy_to = data.join(&label);
        log.info("PredictionVisualizer", &format!("copying [{}] to [{}]", abs.display(), copy_to.display()));
        std::fs::copy(file, &copy_to).map_err(|e| ItemError::new("IOException", format!("cannot copy protein file: {}", e)))?;
        let protein_rel = format!("data/{}", label);
        std::fs::write(vis.join(format!("{}_pymol.pml", label)), output::pymol_script(&p, &pockets, &protein_rel, &points_rel))
            .map_err(|e| ItemError::new("IOException", format!("cannot write pymol script: {}", e)))?;
        std::fs::write(vis.join(format!("{}_chimerax.cxc", label)), output::chimerax_script(&p, &pockets, &protein_rel, &points_rel))
            .map_err(|e| ItemError::new("IOException", format!("cannot write chimerax script: {}", e)))?;
    }
    tick("outputs");
    Ok(ItemResult { label })
}
