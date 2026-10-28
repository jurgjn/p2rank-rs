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

/// `-mode fast` is approximate: check that it agrees closely with the Java output on 1fbl
/// (same pockets at the same places, residue scores close).
#[test]
fn fast_mode_close_to_java() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let case = root.join("tests").join("cases").join("1fbl.pdb");
    let out = std::env::temp_dir().join(format!("prank-rs-test-fast-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_prank-rs"))
        .args(["predict", "-mode", "fast", "-f"])
        .arg(case.join("1fbl.pdb"))
        .arg("-o")
        .arg(&out)
        .env("P2RANK_HOME", root.join("..").join("distro"))
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run prank-rs");
    assert!(status.success());
    let rows = |text: String| -> Vec<Vec<String>> {
        text.lines().skip(1).map(|l| l.split(',').map(|c| c.trim().to_string()).collect()).collect()
    };
    let exp = rows(std::fs::read_to_string(case.join("expected_1fbl.pdb_predictions.csv")).unwrap());
    let got = rows(std::fs::read_to_string(out.join("1fbl.pdb_predictions.csv")).unwrap());
    assert_eq!(exp.len(), got.len(), "pocket count");
    for (e, g) in exp.iter().zip(&got) {
        let d: f64 = (6..9).map(|k| (e[k].parse::<f64>().unwrap() - g[k].parse::<f64>().unwrap()).powi(2)).sum::<f64>().sqrt();
        assert!(d < 1.0, "pocket {} center moved by {} A", e[0], d);
    }
    let exp = rows(std::fs::read_to_string(case.join("expected_1fbl.pdb_residues.csv")).unwrap());
    let got = rows(std::fs::read_to_string(out.join("1fbl.pdb_residues.csv")).unwrap());
    assert_eq!(exp.len(), got.len());
    for (e, g) in exp.iter().zip(&got) {
        let d = (e[3].parse::<f64>().unwrap() - g[3].parse::<f64>().unwrap()).abs();
        assert!(d < 0.5, "residue {} score differs by {}", e[1], d);
    }
    let _ = std::fs::remove_dir_all(&out);
}
