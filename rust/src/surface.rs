//! Solvent accessible surface: port of faster-molecular-surface `DistinctPackedNumericalSurfaceV4`
//! (icosahedral tessellation, one point per surviving distinct direction, no sparsification).

use crate::atoms::P3;

pub struct Directions {
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
}

fn parse_tess(text: &str) -> Vec<P3> {
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let v: Vec<f64> = l.split(' ').map(|h| f64::from_bits(u64::from_str_radix(h, 16).unwrap())).collect();
            [v[0], v[1], v[2]]
        })
        .collect()
}

/// Distinct directions of CDK's icosahedral `Tessellate("ico", level)` in first-occurrence order.
pub fn directions(level: i32) -> Result<Directions, String> {
    let text = match level {
        1 => include_str!("data/tess1.txt"),
        2 => include_str!("data/tess2.txt"),
        3 => include_str!("data/tess3.txt"),
        4 => include_str!("data/tess4.txt"),
        _ => return Err(format!("unsupported tessellation level {}", level)),
    };
    let pts = parse_tess(text);
    let mut seen: Vec<[u64; 3]> = Vec::new();
    let mut d = Directions { dx: Vec::new(), dy: Vec::new(), dz: Vec::new() };
    for p in pts {
        let key = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
        if !seen.contains(&key) {
            seen.push(key);
            d.dx.push(p[0]);
            d.dy.push(p[1]);
            d.dz.push(p[2]);
        }
    }
    Ok(d)
}

/// Computes SAS points for atoms with given centers and van der Waals radii.
/// Returns points in atom order, then direction order.
///
/// Each atom's surviving directions depend only on its own neighbor set (all atoms j != i with
/// d2 < (ri + rj)^2), so atoms are processed in parallel chunks and concatenated in order.
pub fn compute_surface_points(centers: &[P3], vdw: &[f64], solvent_radius: f64, tess_level: i32, threads: usize) -> Result<Vec<P3>, String> {
    let n = centers.len();
    let dirs = directions(tess_level)?;

    let mut radius = vec![0.0f64; n];
    let mut radius2 = vec![0.0f64; n];
    let mut max_radius = 0.0f64;
    for i in 0..n {
        let r = vdw[i] + solvent_radius;
        radius[i] = r;
        radius2[i] = r * r;
        if r > max_radius {
            max_radius = r;
        }
    }

    // uniform cell grid with edge >= any (ri + rj), so all neighbors are in adjacent cells
    let cell = 2.0 * (max_radius + solvent_radius);
    let key = |p: &P3| -> (i64, i64, i64) { ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64, (p[2] / cell).floor() as i64) };
    let mut cells: std::collections::HashMap<(i64, i64, i64), Vec<u32>> = std::collections::HashMap::new();
    for (i, p) in centers.iter().enumerate() {
        cells.entry(key(p)).or_default().push(i as u32);
    }

    let atom_points = |i: usize, out: &mut Vec<P3>, nb: &mut Vec<(f64, f64, f64, f64)>| {
        let total_radius = radius[i];
        let total_radius2 = radius2[i];
        let twice = 2.0 * total_radius;
        let (ax, ay, az) = (centers[i][0], centers[i][1], centers[i][2]);
        nb.clear();
        let (kx, ky, kz) = key(&centers[i]);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = cells.get(&(kx + dx, ky + dy, kz + dz)) {
                        for &j in list {
                            let j = j as usize;
                            if j == i {
                                continue;
                            }
                            let x12 = centers[j][0] - ax;
                            let y12 = centers[j][1] - ay;
                            let z12 = centers[j][2] - az;
                            let d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            let sum_r = radius[i] + radius[j];
                            if d2 < sum_r * sum_r {
                                nb.push((x12, y12, z12, (d2 + total_radius2 - radius2[j]) / twice));
                            }
                        }
                    }
                }
            }
        }
        let mut last: Option<usize> = None;
        for d in 0..dirs.dx.len() {
            let px = dirs.dx[d];
            let py = dirs.dy[d];
            let pz = dirs.dz[d];
            if let Some(l) = last {
                let (x, y, z, t) = nb[l];
                if x * px + y * py + z * pz > t {
                    continue;
                }
            }
            let mut buried = false;
            for (k, &(x, y, z, t)) in nb.iter().enumerate() {
                if x * px + y * py + z * pz > t {
                    last = Some(k);
                    buried = true;
                    break;
                }
            }
            if !buried {
                out.push([total_radius * px + ax, total_radius * py + ay, total_radius * pz + az]);
            }
        }
    };

    let threads = threads.max(1).min(n.div_ceil(1024).max(1));
    if threads <= 1 {
        let mut out = Vec::with_capacity(n * 8);
        let mut nb = Vec::new();
        for i in 0..n {
            atom_points(i, &mut out, &mut nb);
        }
        return Ok(out);
    }
    let chunk = n.div_ceil(threads);
    let parts: Vec<Vec<P3>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let atom_points = &atom_points;
                sc.spawn(move || {
                    let mut out = Vec::new();
                    let mut nb = Vec::new();
                    for i in t * chunk..((t + 1) * chunk).min(n) {
                        atom_points(i, &mut out, &mut nb);
                    }
                    out
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    Ok(parts.concat())
}
