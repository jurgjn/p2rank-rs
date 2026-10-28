//! End-to-end regression tests: for every directory in tests/cases/<label>/ prank-rs is run on the
//! structure file <label> and its outputs must be byte-identical to the outputs of the Java P2Rank
//! (stored next to it as expected_*).

use std::path::{Path, PathBuf};

fn run_case(case: &Path) {
    let label = case.file_name().unwrap().to_str().unwrap().to_string();
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let distro = root.join("..").join("distro");
    let out = std::env::temp_dir().join(format!("prank-rs-test-{}-{}", std::process::id(), label));
    let _ = std::fs::remove_dir_all(&out);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_prank-rs"))
        .args(["predict", "-f"])
        .arg(case.join(&label))
        .arg("-o")
        .arg(&out)
        .env("P2RANK_HOME", &distro)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run prank-rs");
    assert!(status.success(), "prank-rs failed on {}", label);
    for (suffix, sub) in [
        ("predictions.csv", ""),
        ("residues.csv", ""),
        ("pymol.pml", "visualizations/"),
        ("chimerax.cxc", "visualizations/"),
        ("points.pdb.gz", "visualizations/data/"),
    ] {
        let name = format!("{}_{}", label, suffix);
        let expected = std::fs::read(case.join(format!("expected_{}", name))).unwrap();
        let got = std::fs::read(out.join(sub).join(&name)).unwrap_or_default();
        assert!(expected == got, "{}{} differs from the Java output", sub, name);
    }
    let _ = std::fs::remove_dir_all(&out);
}

#[test]
fn outputs_identical_to_java() {
    let cases = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("cases");
    let mut n = 0;
    for e in std::fs::read_dir(&cases).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            run_case(&p);
            n += 1;
        }
    }
    assert!(n > 0);
}

/// Runs `prank-rs predict -mode <mode>` on 1fbl; returns (expected, got) rows of the
/// predictions and residues CSV files.
type Rows = Vec<Vec<String>>;
fn run_1fbl(mode: &str) -> ((Rows, Rows), (Rows, Rows)) {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let case = root.join("tests").join("cases").join("1fbl.pdb");
    let out = std::env::temp_dir().join(format!("prank-rs-test-{}-{}", mode, std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_prank-rs"))
        .args(["predict", "-mode", mode, "-f"])
        .arg(case.join("1fbl.pdb"))
        .arg("-o")
        .arg(&out)
        .env("P2RANK_HOME", root.join("..").join("distro"))
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run prank-rs");
    assert!(status.success());
    let rows = |text: String| -> Rows { text.lines().skip(1).map(|l| l.split(',').map(|c| c.trim().to_string()).collect()).collect() };
    let read = |dir: &std::path::Path, name: &str| rows(std::fs::read_to_string(dir.join(name)).unwrap());
    let res = (
        (read(&case, "expected_1fbl.pdb_predictions.csv"), read(&out, "1fbl.pdb_predictions.csv")),
        (read(&case, "expected_1fbl.pdb_residues.csv"), read(&out, "1fbl.pdb_residues.csv")),
    );
    let _ = std::fs::remove_dir_all(&out);
    res
}

fn center_dist(e: &[String], g: &[String]) -> f64 {
    (6..9).map(|k| (e[k].parse::<f64>().unwrap() - g[k].parse::<f64>().unwrap()).powi(2)).sum::<f64>().sqrt()
}

fn residue_diffs(exp: &Rows, got: &Rows) -> Vec<f64> {
    assert_eq!(exp.len(), got.len());
    exp.iter().zip(got).map(|(e, g)| (e[3].parse::<f64>().unwrap() - g[3].parse::<f64>().unwrap()).abs()).collect()
}

/// `-mode fast` is approximate: check that it agrees closely with the Java output on 1fbl
/// (same pockets at the same places, residue scores close).
#[test]
fn fast_mode_close_to_java() {
    let ((exp, got), (exp_r, got_r)) = run_1fbl("fast");
    assert_eq!(exp.len(), got.len(), "pocket count");
    for (e, g) in exp.iter().zip(&got) {
        let d = center_dist(e, g);
        assert!(d < 1.0, "pocket {} center moved by {} A", e[0], d);
    }
    for (e, d) in exp_r.iter().zip(residue_diffs(&exp_r, &got_r)) {
        assert!(d < 0.5, "residue {} score differs by {}", e[1], d);
    }
}

/// `-mode faster` trades accuracy: check that the top pockets are at the same places, the pocket
/// count is about the same and residue scores are close on average.
#[test]
fn faster_mode_close_to_java() {
    let ((exp, got), (exp_r, got_r)) = run_1fbl("faster");
    assert!((exp.len() as i64 - got.len() as i64).abs() <= 1, "pocket count {} vs {}", got.len(), exp.len());
    for (e, g) in exp.iter().zip(&got).take(3) {
        let d = center_dist(e, g);
        assert!(d < 2.0, "pocket {} center moved by {} A", e[0], d);
    }
    let diffs = residue_diffs(&exp_r, &got_r);
    let mean = diffs.iter().sum::<f64>() / diffs.len() as f64;
    assert!(mean < 0.1, "mean residue score difference {}", mean);
}

/// `-descriptors 1` adds `<label>_pockets_descriptors.csv` (one row per pocket) and leaves the
/// other output files identical to the Java output.
#[test]
fn descriptors_add_a_file_only() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let case = root.join("tests").join("cases").join("1fbl.pdb");
    let out = std::env::temp_dir().join(format!("prank-rs-test-descriptors-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_prank-rs"))
        .args(["predict", "-descriptors", "1", "-f"])
        .arg(case.join("1fbl.pdb"))
        .arg("-o")
        .arg(&out)
        .env("P2RANK_HOME", root.join("..").join("distro"))
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run prank-rs");
    assert!(status.success());
    for name in ["1fbl.pdb_predictions.csv", "1fbl.pdb_residues.csv"] {
        let exp = std::fs::read(case.join(format!("expected_{}", name))).unwrap();
        assert!(exp == std::fs::read(out.join(name)).unwrap(), "{} differs", name);
    }
    let pockets = std::fs::read_to_string(out.join("1fbl.pdb_predictions.csv")).unwrap().lines().count() - 1;
    let desc = std::fs::read_to_string(out.join("1fbl.pdb_pockets_descriptors.csv")).unwrap();
    assert_eq!(desc.lines().count() - 1, pockets);
    assert!(desc.lines().nth(1).unwrap().starts_with("pocket1, 1, "));
    let _ = std::fs::remove_dir_all(&out);
}
