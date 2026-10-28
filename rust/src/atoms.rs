//! Point/atom set operations mirroring P2Rank's `Atoms` class (cutoutSphere, cutoutShell, ...),
//! including its kd-tree vs. serial-scan decisions, which determine result order.

use crate::kdtree::{jmax, jmin, KdTree3D};
use std::sync::OnceLock;

/// Params.use_kdtree_cutout_sphere_thrashold
pub const USE_KDTREE_CUTOUT_SPHERE_THRESHOLD: usize = 150;
/// Atoms.KD_TREE_THRESHOLD
pub const KD_TREE_THRESHOLD: usize = 15;

pub type P3 = [f64; 3];

#[inline(always)]
pub fn sqr_dist(a: &P3, b: &P3) -> f64 {
    let x = a[0] - b[0];
    let y = a[1] - b[1];
    let z = a[2] - b[2];
    x * x + y * y + z * z
}

#[inline(always)]
pub fn dist(a: &P3, b: &P3) -> f64 {
    sqr_dist(a, b).sqrt()
}

/// An ordered set of points with payload ids (e.g. atom indices), with a lazily built kd-tree.
pub struct PointSet {
    pub pts: Vec<P3>,
    pub ids: Vec<u32>,
    kd: OnceLock<KdTree3D>,
}

impl Clone for PointSet {
    fn clone(&self) -> Self {
        PointSet::new(self.pts.clone(), self.ids.clone())
    }
}

impl PointSet {
    pub fn new(pts: Vec<P3>, ids: Vec<u32>) -> PointSet {
        debug_assert_eq!(pts.len(), ids.len());
        PointSet { pts, ids, kd: OnceLock::new() }
    }

    pub fn empty() -> PointSet {
        PointSet::new(Vec::new(), Vec::new())
    }

    pub fn len(&self) -> usize {
        self.pts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    pub fn kd(&self) -> &KdTree3D {
        self.kd.get_or_init(|| KdTree3D::build(&self.pts))
    }

    /// Atoms.withKdTreeConditional()
    pub fn with_kd_tree_conditional(&self) {
        if self.len() > KD_TREE_THRESHOLD {
            self.kd();
        }
    }

    fn subset(&self, idx: impl Iterator<Item = usize>) -> PointSet {
        let mut pts = Vec::new();
        let mut ids = Vec::new();
        for i in idx {
            pts.push(self.pts[i]);
            ids.push(self.ids[i]);
        }
        PointSet::new(pts, ids)
    }

    /// Atoms.cutoutSphere
    pub fn cutout_sphere(&self, c: &P3, radius: f64) -> PointSet {
        let sqr = radius * radius;
        if self.len() >= USE_KDTREE_CUTOUT_SPHERE_THRESHOLD {
            let found = self.kd().find_within(c[0], c[1], c[2], sqr);
            self.subset(found.into_iter().map(|i| i as usize))
        } else {
            let idx: Vec<usize> = (0..self.len()).filter(|&i| sqr_dist(&self.pts[i], c) <= sqr).collect();
            self.subset(idx.into_iter())
        }
    }

    /// Atoms.cutoutSphere, returning positions in this set
    pub fn cutout_sphere_idx(&self, c: &P3, radius: f64, out: &mut Vec<u32>) {
        out.clear();
        let sqr = radius * radius;
        if self.len() >= USE_KDTREE_CUTOUT_SPHERE_THRESHOLD {
            self.kd().for_each_within(c[0], c[1], c[2], sqr, |i| out.push(i));
        } else {
            for i in 0..self.len() {
                if sqr_dist(&self.pts[i], c) <= sqr {
                    out.push(i as u32);
                }
            }
        }
    }

    /// Atoms.countWithinSphere
    pub fn count_within_sphere(&self, c: &P3, radius: f64) -> usize {
        let sqr = radius * radius;
        if self.len() >= USE_KDTREE_CUTOUT_SPHERE_THRESHOLD {
            self.kd().count_within(c[0], c[1], c[2], sqr)
        } else {
            (0..self.len()).filter(|&i| sqr_dist(&self.pts[i], c) <= sqr).count()
        }
    }

    /// Atoms.sqrDist(Atom): squared distance to the nearest point
    pub fn sqr_dist_to(&self, a: &P3) -> f64 {
        if let Some(kd) = self.kd.get() {
            kd.nearest_sqr_dist(a[0], a[1], a[2])
        } else {
            // PerfUtils.sqrDistL
            if self.is_empty() {
                return f64::MAX;
            }
            let mut min = f64::MAX;
            for p in &self.pts {
                let dx = a[0] - p[0];
                let dy = a[1] - p[1];
                let dz = a[2] - p[2];
                let next = dx * dx + dy * dy + dz * dz;
                if next < min {
                    min = next;
                }
            }
            min
        }
    }

    /// Atoms.areWithinDistance(Atom, dist)
    pub fn are_within_distance(&self, a: &P3, d: f64) -> bool {
        if let (Some(kd), true) = (self.kd.get(), self.len() > KD_TREE_THRESHOLD) {
            kd.nearest_sqr_dist(a[0], a[1], a[2]).sqrt() <= d
        } else {
            let sqr = d * d;
            self.pts.iter().any(|p| {
                let dx = a[0] - p[0];
                let dy = a[1] - p[1];
                let dz = a[2] - p[2];
                dx * dx + dy * dy + dz * dz <= sqr
            })
        }
    }

    /// Box.aroundAtoms(this).center and distance from center to max corner
    pub fn bounding_box(&self) -> (P3, P3) {
        if self.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let b = self.pts[0];
        let (mut min_x, mut min_y, mut min_z) = (b[0], b[1], b[2]);
        let (mut max_x, mut max_y, mut max_z) = (min_x, min_y, min_z);
        for p in &self.pts[1..] {
            min_x = jmin(min_x, p[0]);
            max_x = jmax(max_x, p[0]);
            min_y = jmin(min_y, p[1]);
            max_y = jmax(max_y, p[1]);
            min_z = jmin(min_z, p[2]);
            max_z = jmax(max_z, p[2]);
        }
        ([min_x, min_y, min_z], [max_x, max_y, max_z])
    }

    /// Atoms.cutoutShell(aroundAtoms, dist)
    pub fn cutout_shell(&self, around: &PointSet, d: f64) -> PointSet {
        if around.is_empty() {
            return PointSet::empty();
        }
        if around.len() == 1 {
            return self.cutout_sphere(&around.pts[0], d);
        }
        let (min, max) = around.bounding_box();
        let center = [(min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0, (min[2] + max[2]) / 2.0];
        let additional = dist(&center, &max);
        let of_atoms = self.cutout_sphere(&center, d + additional);
        cutout_shell_static(&of_atoms, around, d)
    }

    /// geometric centroid (PerfUtils.calculateCentroid)
    pub fn centroid(&self) -> P3 {
        let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
        for p in &self.pts {
            x += p[0];
            y += p[1];
            z += p[2];
        }
        let n = self.len() as f64;
        [x / n, y / n, z / n]
    }
}

/// static Atoms.cutoutShell(ofAtoms, aroundAtoms, dist)
pub fn cutout_shell_static(of_atoms: &PointSet, around: &PointSet, d: f64) -> PointSet {
    if around.is_empty() {
        return PointSet::empty();
    }
    around.with_kd_tree_conditional();
    let sqr = d * d;
    let idx: Vec<usize> = (0..of_atoms.len()).filter(|&i| around.sqr_dist_to(&of_atoms.pts[i]) <= sqr).collect();
    of_atoms.subset(idx.into_iter())
}

/// Atoms.cutoutShell with the per-atom distance filter run in parallel (same result and order).
pub fn cutout_shell_par(of: &PointSet, around: &PointSet, d: f64, threads: usize) -> PointSet {
    if around.is_empty() {
        return PointSet::empty();
    }
    if around.len() == 1 {
        return of.cutout_sphere(&around.pts[0], d);
    }
    let (min, max) = around.bounding_box();
    let center = [(min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0, (min[2] + max[2]) / 2.0];
    let additional = dist(&center, &max);
    let of_atoms = of.cutout_sphere(&center, d + additional);
    let sqr = d * d;
    let n = of_atoms.len();
    // "nearest around-point within d" == "any around-point within d": an exact grid test
    let grid = CountGrid::new(&around.pts, d);
    if grid.is_none() {
        around.with_kd_tree_conditional();
    }
    let test = |q: &P3| -> bool {
        match &grid {
            Some(g) => g.any_within(q, d),
            None => around.sqr_dist_to(q) <= sqr,
        }
    };
    let threads = threads.max(1).min(n.div_ceil(2048).max(1));
    let keep: Vec<bool> = if threads <= 1 {
        (0..n).map(|i| test(&of_atoms.pts[i])).collect()
    } else {
        let chunk = n.div_ceil(threads);
        std::thread::scope(|sc| {
            let hs: Vec<_> = (0..threads)
                .map(|t| {
                    let of_atoms = &of_atoms;
                    let test = &test;
                    sc.spawn(move || (t * chunk..((t + 1) * chunk).min(n)).map(|i| test(&of_atoms.pts[i])).collect::<Vec<bool>>())
                })
                .collect();
            hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
        })
    };
    of_atoms.subset((0..n).filter(|&i| keep[i]))
}

/// Exact sphere counting on a uniform grid: whole cells that are safely inside the sphere are
/// added without per-point tests, boundary cells use the same `d2 <= r*r` predicate as the kd-tree.
pub struct CountGrid {
    lo: [f64; 3],
    cell: f64,
    dims: [i64; 3],
    start: Vec<u32>,
    pts: Vec<P3>,
    /// original index of each point (grid order)
    ids: Vec<u32>,
}

impl CountGrid {
    pub fn new(points: &[P3], cell: f64) -> Option<CountGrid> {
        if points.is_empty() || !points.iter().all(|p| p.iter().all(|v| v.is_finite())) {
            return None;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in points {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let mut dims = [1i64; 3];
        for k in 0..3 {
            dims[k] = ((hi[k] - lo[k]) / cell).floor() as i64 + 1;
        }
        let ncell = dims[0].checked_mul(dims[1])?.checked_mul(dims[2])?;
        if ncell > (64 * points.len() as i64).max(1 << 20) {
            return None;
        }
        let idx = |p: &P3| -> usize {
            let c: Vec<i64> = (0..3).map(|k| (((p[k] - lo[k]) / cell).floor() as i64).clamp(0, dims[k] - 1)).collect();
            ((c[2] * dims[1] + c[1]) * dims[0] + c[0]) as usize
        };
        let mut start = vec![0u32; ncell as usize + 1];
        for p in points {
            start[idx(p) + 1] += 1;
        }
        for c in 0..ncell as usize {
            start[c + 1] += start[c];
        }
        let mut cursor = start.clone();
        let mut sorted = vec![[0.0; 3]; points.len()];
        let mut ids = vec![0u32; points.len()];
        for (i, p) in points.iter().enumerate() {
            let c = idx(p);
            sorted[cursor[c] as usize] = *p;
            ids[cursor[c] as usize] = i as u32;
            cursor[c] += 1;
        }
        Some(CountGrid { lo, cell, dims, start, pts: sorted, ids })
    }

    /// number of points with (q - p)^2 <= radius^2 (same arithmetic as KdTree3D.countWithinRadius)
    pub fn count_within(&self, q: &P3, radius: f64) -> usize {
        let sqr = radius * radius;
        let inner = radius - 1e-6 * (1.0 + radius);
        let inner2 = inner * inner;
        let mut count = 0usize;
        // per axis: index range and squared min / max distance from q to each cell slab
        let mut rng = [(0i64, 0i64); 3];
        let mut near = [[0.0f64; 64]; 3];
        let mut far = [[0.0f64; 64]; 3];
        for k in 0..3 {
            let a = (((q[k] - radius - self.lo[k]) / self.cell).floor() as i64).max(0);
            let b = (((q[k] + radius - self.lo[k]) / self.cell).floor() as i64).min(self.dims[k] - 1);
            if a > b || b - a >= 64 {
                return self.count_slow(q, radius);
            }
            rng[k] = (a, b);
            for i in a..=b {
                let lo = self.lo[k] + i as f64 * self.cell;
                let hi = lo + self.cell;
                let f = (q[k] - lo).abs().max((q[k] - hi).abs());
                let n = if q[k] < lo { lo - q[k] } else if q[k] > hi { q[k] - hi } else { 0.0 };
                far[k][(i - a) as usize] = f * f;
                near[k][(i - a) as usize] = n * n;
            }
        }
        let (x0, x1) = rng[0];
        let (y0, y1) = rng[1];
        let (z0, z1) = rng[2];
        for z in z0..=z1 {
            let (nz, fz) = (near[2][(z - z0) as usize], far[2][(z - z0) as usize]);
            for y in y0..=y1 {
                let (ny, fy) = (near[1][(y - y0) as usize], far[1][(y - y0) as usize]);
                if nz + ny > sqr * 1.000001 + 1e-9 {
                    continue;
                }
                let row = (z * self.dims[1] + y) * self.dims[0];
                for x in x0..=x1 {
                    let ix = (x - x0) as usize;
                    // cells entirely outside the sphere (with a safety margin) cannot contain hits
                    if nz + ny + near[0][ix] > sqr * 1.000001 + 1e-9 {
                        continue;
                    }
                    let c = (row + x) as usize;
                    let (s, e) = (self.start[c] as usize, self.start[c + 1] as usize);
                    if s == e {
                        continue;
                    }
                    if fz + fy + far[0][ix] < inner2 {
                        count += e - s;
                        continue;
                    }
                    for p in &self.pts[s..e] {
                        let dx = q[0] - p[0];
                        let dy = q[1] - p[1];
                        let dz = q[2] - p[2];
                        if dx * dx + dy * dy + dz * dz <= sqr {
                            count += 1;
                        }
                    }
                }
            }
        }
        count
    }

    pub fn new_indexed(points: &[P3], cell: f64) -> Option<CountGrid> {
        CountGrid::new(points, cell)
    }

    /// calls f(original index) for every point with (q - p)^2 <= radius^2 (grid order)
    pub fn for_each_within(&self, q: &P3, radius: f64, mut f: impl FnMut(u32)) {
        let sqr = radius * radius;
        let mut rng = [(0i64, 0i64); 3];
        for k in 0..3 {
            let a = (((q[k] - radius - self.lo[k]) / self.cell).floor() as i64).max(0);
            let b = (((q[k] + radius - self.lo[k]) / self.cell).floor() as i64).min(self.dims[k] - 1);
            if a > b {
                return;
            }
            rng[k] = (a, b);
        }
        for z in rng[2].0..=rng[2].1 {
            for y in rng[1].0..=rng[1].1 {
                let row = (z * self.dims[1] + y) * self.dims[0];
                let c0 = (row + rng[0].0) as usize;
                let c1 = (row + rng[0].1) as usize;
                let (s, e) = (self.start[c0] as usize, self.start[c1 + 1] as usize);
                for j in s..e {
                    let p = &self.pts[j];
                    let dx = q[0] - p[0];
                    let dy = q[1] - p[1];
                    let dz = q[2] - p[2];
                    if dx * dx + dy * dy + dz * dz <= sqr {
                        f(self.ids[j]);
                    }
                }
            }
        }
    }

    /// true iff some point has (q - p)^2 <= radius^2 (equivalent to `nearest_sqr_dist <= radius^2`)
    pub fn any_within(&self, q: &P3, radius: f64) -> bool {
        let sqr = radius * radius;
        let lim = sqr * 1.000001 + 1e-9;
        let mut rng = [(0i64, 0i64); 3];
        for k in 0..3 {
            let a = (((q[k] - radius - self.lo[k]) / self.cell).floor() as i64).max(0);
            let b = (((q[k] + radius - self.lo[k]) / self.cell).floor() as i64).min(self.dims[k] - 1);
            if a > b {
                return false;
            }
            rng[k] = (a, b);
        }
        let near = |k: usize, i: i64| -> f64 {
            let lo = self.lo[k] + i as f64 * self.cell;
            let hi = lo + self.cell;
            let n = if q[k] < lo { lo - q[k] } else if q[k] > hi { q[k] - hi } else { 0.0 };
            n * n
        };
        for z in rng[2].0..=rng[2].1 {
            let nz = near(2, z);
            for y in rng[1].0..=rng[1].1 {
                let ny = near(1, y);
                if nz + ny > lim {
                    continue;
                }
                let row = (z * self.dims[1] + y) * self.dims[0];
                for x in rng[0].0..=rng[0].1 {
                    if nz + ny + near(0, x) > lim {
                        continue;
                    }
                    let c = (row + x) as usize;
                    for p in &self.pts[self.start[c] as usize..self.start[c + 1] as usize] {
                        let dx = q[0] - p[0];
                        let dy = q[1] - p[1];
                        let dz = q[2] - p[2];
                        if dx * dx + dy * dy + dz * dz <= sqr {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    fn count_slow(&self, q: &P3, radius: f64) -> usize {
        let sqr = radius * radius;
        self.pts
            .iter()
            .filter(|p| {
                let dx = q[0] - p[0];
                let dy = q[1] - p[1];
                let dz = q[2] - p[2];
                dx * dx + dy * dy + dz * dz <= sqr
            })
            .count()
    }
}
