//! prank-rs: a Rust port of the P2Rank `predict` command (default config and model)
//! that produces output identical to the Java implementation.

mod atoms;
mod bcif;
mod chemcomp;
mod descriptors;
mod cif;
mod elements;
mod elements_data;
mod features;
mod forest;
mod javaser;
mod jdtoa;
mod jfmt;
mod json;
mod jutil;
mod kdtree;
mod output;
mod pdb;
mod pipeline;
mod predict;
mod protein;
mod surface;

use pipeline::{Config, Log};
use std::path::{Path, PathBuf};
use std::time::Instant;

const VERSION: &str = "2.6-alpha.8";

/// `-mode fast`: random forest early exit. SAS points whose partial ligandability score is below
/// 0.05 after 10 trees, or below 0.1 after 20 trees (far below the 0.35 pocket threshold), keep the
/// partial score; all other points get the full 100-tree score.
const FAST_RF_STAGES: &[(usize, f64)] = &[(10, 0.05), (20, 0.1)];

/// `-mode faster` (see `pipeline::faster_scores`), development overrides in parentheses:
/// seed scoring only on the 12 coarse (level-1) of the 42 tessellation directions (P2RS_DEV_SEED=all|coarse)
const FASTER_SEED_COARSE: bool = true;
/// seed only points with protrusion >= this, 0 = off (P2RS_DEV_PROT_GATE)
const FASTER_PROT_GATE: f64 = 0.0;
/// points scoring >= this trigger scoring of their neighborhood (P2RS_DEV_REFINE_SCORE)
const FASTER_REFINE_SCORE: f64 = 0.2;
/// ... and so do seed points scoring >= this (P2RS_DEV_SEED_TRIGGER)
const FASTER_SEED_TRIGGER: f64 = 0.1;
/// neighborhood radius (P2RS_DEV_REFINE_RADIUS)
const FASTER_REFINE_RADIUS: f64 = 4.0;
/// early-exit stages for neighborhood points (P2RS_DEV_REFINE_STAGES="k:cutoff,...")
const FASTER_REFINE_STAGES: &[(usize, f64)] = &[(10, 0.05), (20, 0.1)];
/// number of trees, i.e. the first n of 100 (P2RS_DEV_TREES)
const FASTER_TREES: usize = 100;

/// "k:cutoff,k:cutoff" early-exit stages
fn parse_stages(s: &str) -> Vec<(usize, f64)> {
    s.split(',')
        .filter(|x| !x.is_empty())
        .map(|kv| {
            let (k, c) = kv.split_once(':').expect("stage k:cutoff");
            (k.parse().expect("stage k"), c.parse().expect("stage cutoff"))
        })
        .collect()
}

fn usage() -> ! {
    eprintln!(
        "prank-rs {VERSION} - Rust port of P2Rank `predict` (default model)

usage:
  prank-rs predict -f <structure> [-o <outdir>] [-mode exact|fast|faster] [-threads <n>] [-visualizations 0|1] [-descriptors 0|1] [-ccd reduced|full]
  prank-rs predict <dataset.ds>  [-o <outdir>] [-mode exact|fast|faster] [-threads <n>] [-visualizations 0|1] [-descriptors 0|1] [-ccd reduced|full]

  -mode exact|fast|faster
                      exact (default): output byte-identical to Java P2Rank;
                      fast: faster, small documented differences (see README);
                      faster: faster still, scores only part of the surface (see README)
  -descriptors 0|1    also write per-pocket descriptors (shape, buriedness, pharmacophore
                      composition, hydropathy, charge, ...) to <label>_pockets_descriptors.csv
  -ccd reduced|full   chemical component dictionary used for residue classification:
                      reduced (default) = BioJava's bundled fallback set (P2Rank without network),
                      full = complete wwPDB CCD (P2Rank with access to RCSB)
  prank-rs convert-model [<model_dir>]    # model.zst -> model.flb (fast loading cache)

The P2Rank installation (distro) directory is located via the P2RANK_HOME environment variable,
or the directory containing the executable (or its parents)."
    );
    std::process::exit(1)
}

fn find_install_dir() -> Result<PathBuf, String> {
    let ok = |p: &Path| p.join("models").join("default").is_dir() && p.join("models").join("_score_transform").is_dir();
    if let Ok(h) = std::env::var("P2RANK_HOME") {
        let p = PathBuf::from(h);
        if ok(&p) {
            return Ok(std::fs::canonicalize(&p).unwrap_or(p));
        }
        return Err(format!("P2RANK_HOME={} does not look like a P2Rank installation", p.display()));
    }
    if let Ok(exe) = std::env::current_exe() {
        let mut dir = exe.parent().map(|p| p.to_path_buf());
        for _ in 0..5 {
            if let Some(d) = dir {
                if ok(&d) {
                    return Ok(d);
                }
                if ok(&d.join("distro")) {
                    return Ok(d.join("distro"));
                }
                dir = d.parent().map(|p| p.to_path_buf());
            }
        }
    }
    Err("cannot locate P2Rank installation directory (set P2RANK_HOME)".to_string())
}

fn load_forest(install_dir: &Path, log: &Log) -> Result<forest::Forest, String> {
    let dir = install_dir.join("models").join("default");
    log.info("Model", &format!("Loading model from directory (v3 format): {}", dir.display()));
    let flb = dir.join("model.flb");
    let zst = dir.join("model.zst");
    if flb.is_file() {
        // use the cache only if it is newer than the original model
        let newer = match (std::fs::metadata(&flb).and_then(|m| m.modified()), std::fs::metadata(&zst).and_then(|m| m.modified())) {
            (Ok(a), Ok(b)) => a >= b,
            _ => true,
        };
        if newer {
            match forest::Forest::load_flb(&flb) {
                Ok(f) => return Ok(f),
                Err(e) => eprintln!("[WARN] ignoring model cache: {}", e),
            }
        }
    }
    let f = javaser::load_legacy_flat_forest_zst(&zst)?;
    // best effort: write the fast-loading cache next to the model (ignored if not writable)
    if f.write_flb(&flb).is_ok() {
        log.info("Model", &format!("model cache written to {}", flb.display()));
    }
    Ok(f)
}

fn check_model_features(install_dir: &Path) -> Result<(), String> {
    let path = install_dir.join("models").join("default").join("features.txt");
    if let Ok(text) = std::fs::read_to_string(&path) {
        let header: Vec<&str> = text.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
        if header != features::FEATURE_HEADER.to_vec() {
            return Err(format!("model feature header {} does not match the features computed by prank-rs", path.display()));
        }
    }
    Ok(())
}

fn load_transformers(install_dir: &Path) -> predict::Transformers {
    let st = install_dir.join("models").join("_score_transform");
    let load = |p: PathBuf| match predict::ScoreTransformer::load(&p) {
        Ok(t) => Some(t),
        Err(e) => {
            eprintln!("[ERROR] Failed to load score transformer '{}': {}", p.display(), e);
            None
        }
    };
    predict::Transformers {
        zscore: load(st.join("default_ZscoreTpTransformer.json")),
        proba: load(st.join("default_ProbabilityScoreTransformer.json")),
        res_zscore: load(st.join("residue").join("default_ZscoreTpTransformer.json")),
        res_proba: load(st.join("residue").join("default_ProbabilityScoreTransformer.json")),
    }
}

/// Dataset file (.ds), parsed like Dataset.loadFromFile: `#` starts a comment, `PARAM.` lines are
/// dataset attributes, `HEADER:` names the columns (default: `protein`). Only the `protein` and
/// `chains` columns influence `predict`.
fn load_dataset(path: &Path) -> Result<Vec<pipeline::Item>, String> {
    let text = std::fs::read_to_string(path).map_err(|_| format!("cannot find dataset file [{}]", path.display()))?;
    let dir = path.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let mut header: Vec<String> = vec!["protein".to_string()];
    let mut items = Vec::new();
    for line in text.lines() {
        let l = match line.find('#') {
            Some(i) => &line[..i],
            None => line,
        };
        let l = jutil::jtrim(l);
        if l.is_empty() || l.starts_with("PARAM.") {
            continue;
        }
        if l.starts_with("HEADER:") {
            header = l.split_whitespace().skip(1).map(|s| s.to_string()).collect();
            continue;
        }
        let cols: Vec<&str> = l.split_whitespace().collect();
        let get = |name: &str| header.iter().position(|c| c == name).and_then(|i| cols.get(i).copied());
        let f = get("protein").ok_or_else(|| format!("invalid dataset line: {}", l))?;
        let p = PathBuf::from(f);
        let path = if p.is_absolute() { p } else { dir.join(p) };
        let chains = match get("chains") {
            None => None,
            Some(c) if c.trim() == "*" => None,
            Some(c) => Some(c.split(',').map(|s| s.to_string()).collect()),
        };
        items.push(pipeline::Item { path, chains });
    }
    if items.is_empty() {
        return Err(format!("Empty dataset [{}]", path.file_name().and_then(|n| n.to_str()).unwrap_or("")));
    }
    if items.iter().any(|it| !it.path.exists()) {
        for it in items.iter().filter(|it| !it.path.exists()) {
            eprintln!("[ERROR] Dataset - file [{}] does not exist", it.path.display());
        }
        return Err("dataset contains invalid files".to_string());
    }
    Ok(items)
}

fn git_head() -> Option<String> {
    match std::process::Command::new("git").args(["rev-parse", "--short", "HEAD"]).stderr(std::process::Stdio::null()).output() {
        Ok(o) => Some(String::from_utf8_lossy(&o.stdout).to_string()),
        Err(_) => None,
    }
}

fn params_txt(install_dir: &Path, threads: usize, visualizations: bool, git: Option<String>) -> String {
    let mut v = format!("version: {}\n", VERSION);
    if let Some(g) = git {
        v.push_str(&format!("git head: {}\n", g));
    }
    let inst = install_dir.display().to_string();
    let mut body = include_str!("data/params_template.txt")
        .replace("{DATASET_BASE_DIR}", &format!("{}/test_data", inst))
        .replace("{OUTPUT_BASE_DIR}", &format!("{}/test_output", inst))
        .replace("{INSTALL_DIR}", &inst)
        .replace("{THREADS}", &threads.to_string());
    if !visualizations {
        body = body.replace("  visualizations = true\n", "  visualizations = false\n");
    }
    v + &body
}

fn format_time(ms: u128) -> String {
    format!("{} hours {} minutes {}.{:03} seconds", ms / 3_600_000, (ms / 60_000) % 60, (ms / 1000) % 60, ms % 1000)
}

fn run_predict(args: &[String]) -> Result<(), String> {
    let start = Instant::now();
    // `git rev-parse` (for params.txt, like P2Rank) can be slow: run it concurrently
    let git_thread = std::thread::spawn(git_head);
    let mut file: Option<String> = None;
    let mut dataset: Option<String> = None;
    let mut outdir: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut visualizations = true;
    let mut full_ccd = false;
    let mut descriptors = false;
    let mut mode = "exact".to_string();
    let mut i = 0;
    while i < args.len() {
        let a = args[i].as_str();
        let val = || -> Result<String, String> { args.get(i + 1).cloned().ok_or_else(|| format!("missing value for {}", a)) };
        match a {
            "-f" => file = Some(val()?),
            "-o" => outdir = Some(val()?),
            "-threads" => threads = Some(val()?.parse().map_err(|_| "invalid -threads value")?),
            "-visualizations" => {
                let v = val()?;
                visualizations = !(v == "0" || v == "false");
            }
            "-mode" => {
                mode = val()?;
                if mode != "exact" && mode != "fast" && mode != "faster" {
                    return Err(format!("invalid -mode '{}' (exact|fast|faster)", mode));
                }
            }
            "-descriptors" => {
                let v = val()?;
                descriptors = v == "1" || v == "true";
            }
            "-ccd" => {
                full_ccd = match val()?.as_str() {
                    "full" => true,
                    "reduced" => false,
                    v => return Err(format!("invalid -ccd value '{}' (full|reduced)", v)),
                };
            }
            "-c" => {
                let v = val()?;
                if v != "default" && v != "default.groovy" && !v.ends_with("config/default.groovy") {
                    return Err(format!("only the default config is supported by prank-rs (got -c {})", v));
                }
            }
            "-m" => {
                let v = val()?;
                if v != "default" {
                    return Err(format!("only the default model is supported by prank-rs (got -m {})", v));
                }
            }
            _ if a.starts_with('-') => return Err(format!("unsupported option: {}", a)),
            _ => {
                dataset = Some(a.to_string());
                i += 1;
                continue;
            }
        }
        i += 2;
    }
    let fast = mode == "fast" || mode == "faster";
    let faster = if mode == "faster" {
        // development overrides, see the FASTER_* constants
        let env = |k: &str| std::env::var(k).ok();
        Some(pipeline::FasterParams {
            seed_coarse: env("P2RS_DEV_SEED").map(|v| v != "all").unwrap_or(FASTER_SEED_COARSE),
            prot_gate: env("P2RS_DEV_PROT_GATE").map(|v| v.parse().expect("P2RS_DEV_PROT_GATE")).unwrap_or(FASTER_PROT_GATE),
            refine_score: env("P2RS_DEV_REFINE_SCORE").map(|v| v.parse().expect("P2RS_DEV_REFINE_SCORE")).unwrap_or(FASTER_REFINE_SCORE),
            seed_trigger: env("P2RS_DEV_SEED_TRIGGER").map(|v| v.parse().expect("P2RS_DEV_SEED_TRIGGER")).unwrap_or(FASTER_SEED_TRIGGER),
            refine_radius: env("P2RS_DEV_REFINE_RADIUS").map(|v| v.parse().expect("P2RS_DEV_REFINE_RADIUS")).unwrap_or(FASTER_REFINE_RADIUS),
            trees: env("P2RS_DEV_TREES").map(|v| v.parse().expect("P2RS_DEV_TREES")).unwrap_or(FASTER_TREES),
            refine_stages: env("P2RS_DEV_REFINE_STAGES").map(|v| parse_stages(&v)).unwrap_or_else(|| FASTER_REFINE_STAGES.to_vec()),
        })
    } else {
        None
    };
    let mut rf_stages: Vec<(usize, f64)> = if fast { FAST_RF_STAGES.to_vec() } else { Vec::new() };
    // development override, e.g. P2RS_DEV_STAGES="10:0.05,20:0.1"
    if let Ok(s) = std::env::var("P2RS_DEV_STAGES") {
        rf_stages = parse_stages(&s);
    }
    let install_dir = find_install_dir()?;
    let threads = threads.unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) + 1);

    let (items, ds_name): (Vec<pipeline::Item>, String) = match (&file, &dataset) {
        (Some(f), None) => {
            let p = PathBuf::from(f);
            let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
            (vec![pipeline::Item { path: p, chains: None }], name)
        }
        (None, Some(d)) => {
            let p = PathBuf::from(d);
            let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
            (load_dataset(&p)?, name)
        }
        _ => return Err("specify exactly one of -f <file> or <dataset.ds>".to_string()),
    };
    let ds_label = match ds_name.rsplit_once('.') {
        Some((a, _)) => a.to_string(),
        None => ds_name.clone(),
    };
    let outdir = PathBuf::from(outdir.unwrap_or_else(|| format!("{}/test_output/predict_{}", install_dir.display(), ds_label)));

    let log = Log::new(false);
    log.info("Writable", &format!("predicting pockets for proteins from dataset [{}]", ds_name));
    std::fs::create_dir_all(&outdir).map_err(|e| format!("cannot create {}: {}", outdir.display(), e))?;
    log.info("PredictPocketsRoutine", &format!("outdir: {}", outdir.display()));
    if visualizations {
        std::fs::create_dir_all(outdir.join("visualizations")).map_err(|e| e.to_string())?;
    }

    check_model_features(&install_dir)?;
    let t0 = Instant::now();
    let forest = load_forest(&install_dir, &log)?;
    if std::env::var("P2RS_TIMING").is_ok() {
        eprintln!("timing: model load {:?}", t0.elapsed());
    }
    log.info("FeatureSetup", "effectively enabled features: [chem, volsite, protrusion, bfactor, atom_table]");
    let tr = load_transformers(&install_dir);
    let cfg = Config { install_dir: install_dir.clone(), threads, visualizations, full_ccd, rf_stages: rf_stages.clone(), fast, faster: faster.clone(), descriptors };

    log.info("Dataset", &format!("processing dataset [{}] using {} threads", ds_name, threads));
    let n = items.len();
    // like P2Rank, process dataset items in parallel; a single item uses all threads internally
    let item_threads = if n > 1 { 1 } else { threads };
    let workers = if n > 1 { threads.min(n) } else { 1 };
    let cfg_item = Config { install_dir: cfg.install_dir.clone(), threads: item_threads, visualizations, full_ccd, rf_stages: rf_stages.clone(), fast, faster: faster.clone(), descriptors };
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results: std::sync::Mutex<Vec<(usize, String, pipeline::ItemError)>> = std::sync::Mutex::new(Vec::new());
    std::thread::scope(|sc| {
        for _ in 0..workers {
            sc.spawn(|| loop {
                let k = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if k >= n {
                    break;
                }
                let item = &items[k];
                let lbl = item.path.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
                log.info("Writable", &format!("processing [{}] ({}/{})", lbl, k + 1, n));
                if let Err(e) = pipeline::predict_file(item, &outdir, &cfg_item, &forest, &tr, &log) {
                    log.error("Dataset", &format!("error processing dataset item [{}]: {}", lbl, e.format()));
                    results.lock().unwrap().push((k, lbl, e));
                }
            });
        }
    });
    let mut errs = results.into_inner().unwrap();
    errs.sort_by_key(|e| e.0);
    let item_errors: Vec<(String, pipeline::ItemError)> = errs.into_iter().map(|(_, l, e)| (l, e)).collect();
    let errors = item_errors.len();
    if errors > 0 {
        write_error_reports(&outdir, &item_errors);
    }
    if std::env::var("P2RS_TIMING").is_ok() {
        eprintln!("timing: items done {:?}", start.elapsed());
    }
    let git = git_thread.join().unwrap_or(None);
    std::fs::write(outdir.join("params.txt"), params_txt(&install_dir, threads, visualizations, git))
        .map_err(|e| format!("cannot write params.txt: {}", e))?;
    log.info("Writable", &format!("predicting pockets finished in {}", format_time(start.elapsed().as_millis())));
    let abs_out = std::fs::canonicalize(&outdir).unwrap_or(outdir.clone());
    log.info("Writable", &format!("results saved to directory [{}]", abs_out.display()));
    let mut lines = log.take();
    if errors > 0 {
        lines.push(format!("[ERROR] Main - {} of {} items failed", errors, n));
    } else {
        let t = format_time(start.elapsed().as_millis());
        lines.push("[INFO] Console - ".to_string());
        lines.push(format!("[INFO] Console - Finished successfully in {}.", t));
        println!("\nFinished successfully in {}.", t);
    }
    let _ = std::fs::write(outdir.join("run.log"), lines.join("\n") + "\n");
    if errors > 0 {
        return Err(format!("{} of {} items failed", errors, n));
    }
    Ok(())
}

/// Dataset.Result error reports (errors.csv, errors_aggregated.csv, errors_full.txt.gz)
fn write_error_reports(outdir: &Path, errors: &[(String, pipeline::ItemError)]) {
    let mut sb = String::from("file (dataset.item.label), error (exception.message)\n");
    for (lbl, e) in errors {
        sb.push_str(&format!("{}, {}\n", lbl, e.format()));
    }
    let _ = std::fs::write(outdir.join("errors.csv"), sb);
    let mut counts: Vec<(String, usize)> = Vec::new();
    for (_, e) in errors {
        let m = e.format();
        match counts.iter_mut().find(|(k, _)| *k == m) {
            Some(c) => c.1 += 1,
            None => counts.push((m, 1)),
        }
    }
    counts.sort_by(|a, b| b.1.cmp(&a.1));
    let mut sb = String::from("count, error\n");
    for (m, c) in &counts {
        sb.push_str(&format!("{}, {}\n", c, m));
    }
    let _ = std::fs::write(outdir.join("errors_aggregated.csv"), sb);
    let mut full = String::new();
    for (lbl, e) in errors {
        full.push_str(&format!("ITEM: {}\n", lbl));
        for (i, (c, m)) in e.0.iter().enumerate() {
            if i > 0 {
                full.push_str("Caused by: ");
            }
            full.push_str(&format!("{}: {}\n", c, m));
        }
        full.push_str("\n\n");
    }
    let _ = output::write_gzip_java(&outdir.join("errors_full.txt.gz"), full.as_bytes(), 6);
}

fn run_convert_model(args: &[String]) -> Result<(), String> {
    let dir = match args.first() {
        Some(d) => PathBuf::from(d),
        None => find_install_dir()?.join("models").join("default"),
    };
    let f = javaser::load_legacy_flat_forest_zst(&dir.join("model.zst"))?;
    let out = dir.join("model.flb");
    f.write_flb(&out).map_err(|e| format!("cannot write {}: {}", out.display(), e))?;
    println!("written {}", out.display());
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }
    let res = match args[1].as_str() {
        "predict" => run_predict(&args[2..]),
        "convert-model" => run_convert_model(&args[2..]),
        "-v" | "--version" | "version" => {
            println!("prank-rs {}", VERSION);
            Ok(())
        }
        _ => usage(),
    };
    if let Err(e) = res {
        eprintln!("ERROR: {}", e);
        std::process::exit(1);
    }
}
