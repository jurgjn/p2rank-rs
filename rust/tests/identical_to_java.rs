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
