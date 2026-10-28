//! Per-pocket descriptors (`-descriptors 1`): cheap global pocket properties in the spirit of
//! fpocket / DoGSiteScorer descriptors and FuzCav / VolSite pharmacophore counts, computed from
//! what P2Rank already has for each pocket (SAS points, surface atoms, residues, point scores).
//! Written to `<label>_pockets_descriptors.csv`; does not change any other output.

use crate::atoms::P3;
use crate::elements::{ELEMENT_N, ELEMENT_O};
use crate::features::FeatureContext;
use crate::predict::{LabeledPoint, Pocket};
use crate::protein::Protein;

pub struct Descriptors {
    /// radius of gyration of the pocket SAS points (A)
    pub rgyr: f64,
    /// standard deviations along the principal axes of the pocket SAS points, largest first (A)
    pub axes: [f64; 3],
    pub shape: &'static str,
    /// mean protrusion (deep-layer protein atoms within 10 A) of the pocket SAS points
    pub buriedness: f64,
    /// fractions of the surface atoms with each VolSite pharmacophore property
    /// [aromatic, cation, anion, hydrophobic, acceptor, donor]
    pub pharm: [f64; 6],
    /// fraction of N and O surface atoms
    pub polar: f64,
    /// mean Kyte-Doolittle hydropathy of the pocket residues
    pub hydropathy: f64,
    /// ARG + LYS - ASP - GLU over the pocket residues
    pub net_charge: i32,
    pub residues: usize,
    /// mean B-factor of the surface atoms
    pub bfactor: f64,
    /// max and mean ligandability score of the pocket's ligandable (>= 0.35) points
    pub score_max: f64,
    pub score_mean: f64,
    pub flavour: &'static str,
}

const PHARM_NAMES: [&str; 6] = ["aromatic", "cation", "anion", "hydrophobic", "acceptor", "donor"];

fn kyte_doolittle(code: &str) -> Option<f64> {
    Some(match code {
        "ALA" => 1.8,
        "ARG" => -4.5,
        "ASN" => -3.5,
        "ASP" => -3.5,
        "CYS" => 2.5,
        "GLN" => -3.5,
        "GLU" => -3.5,
        "GLY" => -0.4,
        "HIS" => -3.2,
        "ILE" => 4.5,
        "LEU" => 3.8,
        "LYS" => -3.9,
        "MET" => 1.9,
        "PHE" => 2.8,
        "PRO" => -1.6,
        "SER" => -0.8,
        "THR" => -0.7,
        "TRP" => -0.9,
        "TYR" => -1.3,
        "VAL" => 4.2,
        _ => return None,
    })
}

/// Eigenvalues of a symmetric 3x3 matrix, descending (closed form, Smith 1961).
fn sym3_eigenvalues(m: [[f64; 3]; 3]) -> [f64; 3] {
    let p1 = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
    let mut e = if p1 == 0.0 {
        [m[0][0], m[1][1], m[2][2]]
    } else {
        let q = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
        let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
        let p = (p2 / 6.0).sqrt();
        let mut b = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                b[i][j] = (m[i][j] - if i == j { q } else { 0.0 }) / p;
            }
        }
        let det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
            + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
        let r = (det / 2.0).clamp(-1.0, 1.0);
        let phi = r.acos() / 3.0;
        let e1 = q + 2.0 * p * phi.cos();
        let e3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
        [e1, 3.0 * q - e1 - e3, e3]
    };
    e.sort_by(|a, b| b.total_cmp(a));
    e
}

fn shape_stats(pts: &[P3]) -> (f64, [f64; 3]) {
    let n = pts.len() as f64;
    if pts.is_empty() {
        return (0.0, [0.0; 3]);
    }
    let mut c = [0.0; 3];
    for p in pts {
        for k in 0..3 {
            c[k] += p[k] / n;
        }
    }
    let mut cov = [[0.0; 3]; 3];
    for p in pts {
        let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j] / n;
            }
        }
    }
    let ev = sym3_eigenvalues(cov);
    let rgyr = (cov[0][0] + cov[1][1] + cov[2][2]).sqrt();
    (rgyr, [ev[0].max(0.0).sqrt(), ev[1].max(0.0).sqrt(), ev[2].max(0.0).sqrt()])
}

/// Shape label from the principal axes of the pocket SAS points (thresholds: see README).
fn shape_label(axes: [f64; 3]) -> &'static str {
    let elong = axes[0] / axes[1].max(1e-9);
    let flat = axes[1] / axes[2].max(1e-9);
    if elong >= SHAPE_ELONGATED {
        "tunnel"
    } else if flat >= SHAPE_FLAT {
        "saucer"
    } else {
        "bowl"
    }
}

/// Pocket "flavour": the pharmacophore property most enriched over the average P2Rank pocket,
/// as a z-score against the mean and SD of pockets predicted on ah4h (see README).
fn flavour_label(pharm: [f64; 6]) -> &'static str {
    const LABELS: [&str; 6] = ["aromatic", "basic", "acidic", "greasy", "acceptor-rich", "donor-rich"];
    let mut best = (FLAVOUR_MIN_Z, "balanced");
    for k in 0..6 {
        let z = (pharm[k] - PHARM_MEAN[k]) / PHARM_SD[k];
        if z > best.0 {
            best = (z, LABELS[k]);
        }
    }
    best.1
}

/// principal-axis ratios from which a pocket counts as elongated (largest / middle) or flat
/// (middle / smallest): about the upper quartile of each ratio over the 3575 pockets predicted
/// on ah4h (quartiles 1.29 / 1.50 / 1.83 for both ratios)
const SHAPE_ELONGATED: f64 = 1.8;
const SHAPE_FLAT: f64 = 1.8;
/// mean and SD of the pharmacophore fractions over the 3575 pockets predicted on ah4h
/// (neither training nor evaluation data)
const PHARM_MEAN: [f64; 6] = [0.1482, 0.0292, 0.0691, 0.4570, 0.2024, 0.1445];
const PHARM_SD: [f64; 6] = [0.1386, 0.0400, 0.0735, 0.1317, 0.0887, 0.0840];
/// minimum z-score for a flavour (otherwise "balanced")
const FLAVOUR_MIN_Z: f64 = 1.0;

pub fn compute(p: &Protein, pk: &Pocket, points: &[LabeledPoint], ctx: &FeatureContext) -> Descriptors {
    let pts: Vec<P3> = pk.sas_points.iter().map(|&i| points[i].pos).collect();
    let (rgyr, axes) = shape_stats(&pts);
    let buriedness = if pts.is_empty() { 0.0 } else { pts.iter().map(|q| ctx.protrusion(q)).sum::<f64>() / pts.len() as f64 };

    let na = pk.surface_atoms.len().max(1) as f64;
    let mut pharm = [0.0; 6];
    let (mut polar, mut bfactor) = (0.0, 0.0);
    for &a in &pk.surface_atoms {
        let v = crate::features::atom_pharmacophore(p, a);
        for k in 0..6 {
            pharm[k] += v[k] / na;
        }
        let e = p.s.atoms[a].element;
        if e == ELEMENT_N || e == ELEMENT_O {
            polar += 1.0 / na;
        }
        bfactor += p.s.atoms[a].temp_factor as f64 / na;
    }

    let mut groups: Vec<usize> = pk.surface_atoms.iter().map(|&a| p.s.atoms[a].group).collect();
    groups.sort_unstable();
    groups.dedup();
    let (mut kd_sum, mut kd_n, mut net_charge) = (0.0, 0usize, 0i32);
    for &g in &groups {
        let code = crate::features::correct_residue_code(&p.s.groups[g].pdb_name);
        if let Some(h) = kyte_doolittle(&code) {
            kd_sum += h;
            kd_n += 1;
        }
        net_charge += match code.as_str() {
            "ARG" | "LYS" => 1,
            "ASP" | "GLU" => -1,
            _ => 0,
        };
    }

    // ligandable points of the pocket (its cluster, plus any extended point above the threshold)
    let cluster: Vec<f64> = pk.sas_points.iter().map(|&i| points[i].score).filter(|&s| s >= crate::predict::PRED_POINT_THRESHOLD).collect();
    let score_max = cluster.iter().cloned().fold(0.0, f64::max);
    let score_mean = if cluster.is_empty() { 0.0 } else { cluster.iter().sum::<f64>() / cluster.len() as f64 };

    Descriptors {
        rgyr,
        axes,
        shape: shape_label(axes),
        buriedness,
        pharm,
        polar,
        hydropathy: if kd_n > 0 { kd_sum / kd_n as f64 } else { 0.0 },
        net_charge,
        residues: groups.len(),
        bfactor,
        score_max,
        score_mean,
        flavour: flavour_label(pharm),
    }
}

pub fn csv(p: &Protein, pockets: &[Pocket], points: &[LabeledPoint], ctx: &FeatureContext) -> String {
    let mut s = String::from("name, rank, residues, rgyr, axis1, axis2, axis3, shape, buriedness");
    for n in PHARM_NAMES {
        s.push_str(", ");
        s.push_str(n);
    }
    s.push_str(", polar_atoms, hydropathy, net_charge, bfactor, score_max, score_mean, flavour\n");
    for pk in pockets {
        let d = compute(p, pk, points, ctx);
        s.push_str(&format!(
            "{}, {}, {}, {:.2}, {:.2}, {:.2}, {:.2}, {}, {:.1}",
            pk.name, pk.new_rank, d.residues, d.rgyr, d.axes[0], d.axes[1], d.axes[2], d.shape, d.buriedness
        ));
        for v in d.pharm {
            s.push_str(&format!(", {:.3}", v));
        }
        s.push_str(&format!(
            ", {:.3}, {:.2}, {}, {:.1}, {:.3}, {:.3}, {}\n",
            d.polar, d.hydropathy, d.net_charge, d.bfactor, d.score_max, d.score_mean, d.flavour
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigenvalues_of_known_matrices() {
        let e = sym3_eigenvalues([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(e, [3.0, 2.0, 1.0]);
        // [[2,1,0],[1,2,0],[0,0,5]] has eigenvalues 5, 3, 1
        let e = sym3_eigenvalues([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 5.0]]);
        for (a, b) in e.iter().zip([5.0, 3.0, 1.0]) {
            assert!((a - b).abs() < 1e-9, "{:?}", e);
        }
    }

    #[test]
    fn line_of_points_is_a_tunnel() {
        let pts: Vec<P3> = (0..20).map(|i| [i as f64, (i % 2) as f64 * 0.5, 0.0]).collect();
        let (_, axes) = shape_stats(&pts);
        assert_eq!(shape_label(axes), "tunnel");
    }
}
