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

/// Point tag: atom index << 10 | direction index (level 4 has 642 directions).
pub fn tag_dir(tag: u32) -> usize {
    (tag & 0x3ff) as usize
}

/// Mask of the `coarse` level directions within the `level` direction list (the icosahedral
/// tessellations are nested: every level-1 direction is bitwise a level-2 direction).
pub fn nested_direction_mask(coarse: i32, level: i32) -> Result<Vec<bool>, String> {
    let c = directions(coarse)?;
    let d = directions(level)?;
    Ok((0..d.dx.len())
        .map(|i| (0..c.dx.len()).any(|j| c.dx[j].to_bits() == d.dx[i].to_bits() && c.dy[j].to_bits() == d.dy[i].to_bits() && c.dz[j].to_bits() == d.dz[i].to_bits()))
        .collect())
}

/// Computes SAS points for atoms with given centers and van der Waals radii.
/// Returns points in atom order, then direction order.
///
/// Each atom's surviving directions depend only on its own neighbor set (all atoms j != i with
/// d2 < (ri + rj)^2), so atoms are processed in parallel chunks and concatenated in order.
/// Also returns a tag (atom index, direction index) per point.
pub fn compute_surface_points_tagged(centers: &[P3], vdw: &[f64], solvent_radius: f64, tess_level: i32, threads: usize) -> Result<(Vec<P3>, Vec<u32>), String> {
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

    // dense uniform cell grid (CSR layout); edge = max(ri + rj) so all neighbors are in adjacent cells.
    // The neighbor *set* is what matters (the buried test is an "any" predicate), so the grid
    // geometry does not influence results.
    let cell = (2.0 * max_radius).max(1e-6);
    let mut lo = [f64::INFINITY; 3];
    for p in centers {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
        }
    }
    let cidx = |p: &P3, k: usize| -> i64 { ((p[k] - lo[k]) / cell).floor() as i64 };
    let mut dims = [1i64; 3];
    for p in centers {
        for k in 0..3 {
            dims[k] = dims[k].max(cidx(p, k) + 1);
        }
    }
    let use_dense = n > 0 && centers.iter().all(|p| p.iter().all(|v| v.is_finite())) && dims[0] * dims[1] * dims[2] <= (64 * n as i64).max(1 << 20);
    let ncell = if use_dense { (dims[0] * dims[1] * dims[2]) as usize } else { 0 };
    let cell_of = |p: &P3| -> usize { ((cidx(p, 2) * dims[1] + cidx(p, 1)) * dims[0] + cidx(p, 0)) as usize };
    let mut cell_start = vec![0u32; ncell + 1];
    let mut cell_atoms = vec![0u32; if use_dense { n } else { 0 }];
    let mut hashed: std::collections::HashMap<(i64, i64, i64), Vec<u32>> = std::collections::HashMap::new();
    if use_dense {
        for p in centers {
            cell_start[cell_of(p) + 1] += 1;
        }
        for c in 0..ncell {
            cell_start[c + 1] += cell_start[c];
        }
        let mut cursor = cell_start.clone();
        for (i, p) in centers.iter().enumerate() {
            let c = cell_of(p);
            cell_atoms[cursor[c] as usize] = i as u32;
            cursor[c] += 1;
        }
    } else {
        // fallback for sparse / non-finite coordinates
        let key = |p: &P3| -> (i64, i64, i64) { ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64, (p[2] / cell).floor() as i64) };
        for (i, p) in centers.iter().enumerate() {
            hashed.entry(key(p)).or_default().push(i as u32);
        }
    }
    let for_candidates = |i: usize, f: &mut dyn FnMut(usize)| {
        let p = &centers[i];
        if use_dense {
            let (cx, cy, cz) = (cidx(p, 0), cidx(p, 1), cidx(p, 2));
            for z in (cz - 1).max(0)..=(cz + 1).min(dims[2] - 1) {
                for y in (cy - 1).max(0)..=(cy + 1).min(dims[1] - 1) {
                    let row = (z * dims[1] + y) * dims[0];
                    let c0 = (row + (cx - 1).max(0)) as usize;
                    let c1 = (row + (cx + 1).min(dims[0] - 1)) as usize;
                    for &j in &cell_atoms[cell_start[c0] as usize..cell_start[c1 + 1] as usize] {
                        f(j as usize);
                    }
                }
            }
        } else {
            let k = ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64, (p[2] / cell).floor() as i64);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = hashed.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                            for &j in list {
                                f(j as usize);
                            }
                        }
                    }
                }
            }
        }
    };

    let atom_points = |i: usize, out: &mut Vec<P3>, tags: &mut Vec<u32>, nb: &mut Vec<(f64, f64, f64, f64)>| {
        let total_radius = radius[i];
        let total_radius2 = radius2[i];
        let twice = 2.0 * total_radius;
        let (ax, ay, az) = (centers[i][0], centers[i][1], centers[i][2]);
        nb.clear();
        for_candidates(i, &mut |j| {
            if j == i {
                return;
            }
            let x12 = centers[j][0] - ax;
            let y12 = centers[j][1] - ay;
            let z12 = centers[j][2] - az;
            let d2 = x12 * x12 + y12 * y12 + z12 * z12;
            let sum_r = radius[i] + radius[j];
            if d2 < sum_r * sum_r {
                nb.push((x12, y12, z12, (d2 + total_radius2 - radius2[j]) / twice));
            }
        });
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
                tags.push((i as u32) << 10 | d as u32);
            }
        }
    };

    let threads = threads.max(1).min(n.div_ceil(1024).max(1));
    if threads <= 1 {
        let mut out = Vec::with_capacity(n * 8);
        let mut tags = Vec::with_capacity(n * 8);
        let mut nb = Vec::new();
        for i in 0..n {
            atom_points(i, &mut out, &mut tags, &mut nb);
        }
        return Ok((out, tags));
    }
    let chunk = n.div_ceil(threads);
    let parts: Vec<(Vec<P3>, Vec<u32>)> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let atom_points = &atom_points;
                sc.spawn(move || {
                    let mut out = Vec::new();
                    let mut tags = Vec::new();
                    let mut nb = Vec::new();
                    for i in t * chunk..((t + 1) * chunk).min(n) {
                        atom_points(i, &mut out, &mut tags, &mut nb);
                    }
                    (out, tags)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let (pts, tags): (Vec<Vec<P3>>, Vec<Vec<u32>>) = parts.into_iter().unzip();
    Ok((pts.concat(), tags.concat()))
}
