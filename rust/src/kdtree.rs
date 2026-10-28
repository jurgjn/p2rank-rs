//! Exact port of P2Rank's `KdTree3D` (immutable, median-split kd-tree with 24-point buckets).
//! The traversal order of range queries is part of the observable behaviour (it determines the
//! order of neighbour lists and therefore floating point summation order), so the build and query
//! algorithms are reproduced step by step.

const BUCKET_SIZE: usize = 24;

pub struct KdTree3D {
    xs: Vec<f64>,
    ys: Vec<f64>,
    zs: Vec<f64>,
    /// payload (index into the original point list) of each point in tree order
    ids: Vec<u32>,
    split_dims: Vec<u8>,
    split_vals: Vec<f64>,
    min_xs: Vec<f64>,
    min_ys: Vec<f64>,
    min_zs: Vec<f64>,
    max_xs: Vec<f64>,
    max_ys: Vec<f64>,
    max_zs: Vec<f64>,
    leaf_start: Vec<i32>,
    leaf_end: Vec<i32>,
    node_count: usize,
}

#[inline(always)]
fn sqr_dist(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> f64 {
    let dx = x1 - x2;
    let dy = y1 - y2;
    let dz = z1 - z2;
    dx * dx + dy * dy + dz * dz
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn sqr_dist_to_box(qx: f64, qy: f64, qz: f64, min_x: f64, min_y: f64, min_z: f64, max_x: f64, max_y: f64, max_z: f64) -> f64 {
    let mut d = 0.0;
    if qx < min_x {
        let dx = qx - min_x;
        d += dx * dx;
    } else if qx > max_x {
        let dx = qx - max_x;
        d += dx * dx;
    }
    if qy < min_y {
        let dy = qy - min_y;
        d += dy * dy;
    } else if qy > max_y {
        let dy = qy - max_y;
        d += dy * dy;
    }
    if qz < min_z {
        let dz = qz - min_z;
        d += dz * dz;
    } else if qz > max_z {
        let dz = qz - max_z;
        d += dz * dz;
    }
    d
}

/// Java Math.min / Math.max for doubles (NaN propagating, -0.0 < 0.0)
#[inline(always)]
pub fn jmin(a: f64, b: f64) -> f64 {
    if a != a {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a <= b {
        a
    } else {
        b
    }
}

#[inline(always)]
pub fn jmax(a: f64, b: f64) -> f64 {
    if a != a {
        return a;
    }
    if a == 0.0 && b == 0.0 && a.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a >= b {
        a
    } else {
        b
    }
}

struct Builder<'a> {
    xs: &'a mut [f64],
    ys: &'a mut [f64],
    zs: &'a mut [f64],
    ids: &'a mut [u32],
    split_dims: Vec<u8>,
    split_vals: Vec<f64>,
    min_xs: Vec<f64>,
    min_ys: Vec<f64>,
    min_zs: Vec<f64>,
    max_xs: Vec<f64>,
    max_ys: Vec<f64>,
    max_zs: Vec<f64>,
    leaf_start: Vec<i32>,
    leaf_end: Vec<i32>,
}

impl<'a> Builder<'a> {
    #[inline(always)]
    fn swap(&mut self, i: usize, j: usize) {
        self.xs.swap(i, j);
        self.ys.swap(i, j);
        self.zs.swap(i, j);
        self.ids.swap(i, j);
    }

    #[inline(always)]
    fn key(&self, dim: u8, i: usize) -> f64 {
        match dim {
            0 => self.xs[i],
            1 => self.ys[i],
            _ => self.zs[i],
        }
    }

    fn quickselect(&mut self, mut lo: usize, mut hi: usize, k: usize, dim: u8) {
        while lo < hi {
            if hi - lo == 1 {
                if self.key(dim, lo) > self.key(dim, hi) {
                    self.swap(lo, hi);
                }
                return;
            }
            let mid = lo + (hi - lo) / 2;
            if self.key(dim, lo) > self.key(dim, mid) {
                self.swap(lo, mid);
            }
            if self.key(dim, lo) > self.key(dim, hi) {
                self.swap(lo, hi);
            }
            if self.key(dim, mid) > self.key(dim, hi) {
                self.swap(mid, hi);
            }
            self.swap(mid, hi - 1);
            let pivot = self.key(dim, hi - 1);

            let mut i = lo as isize;
            let mut j = (hi - 1) as isize;
            loop {
                loop {
                    i += 1;
                    if !(self.key(dim, i as usize) < pivot) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if !(self.key(dim, j as usize) > pivot) {
                        break;
                    }
                }
                if i >= j {
                    break;
                }
                self.swap(i as usize, j as usize);
            }
            let i = i as usize;
            self.swap(i, hi - 1);

            let mut eq_lo = i;
            while eq_lo > lo && self.key(dim, eq_lo - 1) == pivot {
                eq_lo -= 1;
            }
            let mut eq_hi = i;
            while eq_hi < hi && self.key(dim, eq_hi + 1) == pivot {
                eq_hi += 1;
            }
            if k >= eq_lo && k <= eq_hi {
                return;
            }
            if k < eq_lo {
                hi = eq_lo - 1;
            } else {
                lo = eq_hi + 1;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_node(
        &mut self,
        node: usize,
        from: usize,
        to: usize,
        p_min_x: f64,
        p_min_y: f64,
        p_min_z: f64,
        p_max_x: f64,
        p_max_y: f64,
        p_max_z: f64,
    ) -> usize {
        let count = to - from;
        if count <= BUCKET_SIZE {
            let (mut x_min, mut y_min, mut z_min) = (self.xs[from], self.ys[from], self.zs[from]);
            let (mut x_max, mut y_max, mut z_max) = (x_min, y_min, z_min);
            for i in from + 1..to {
                if self.xs[i] < x_min {
                    x_min = self.xs[i];
                } else if self.xs[i] > x_max {
                    x_max = self.xs[i];
                }
                if self.ys[i] < y_min {
                    y_min = self.ys[i];
                } else if self.ys[i] > y_max {
                    y_max = self.ys[i];
                }
                if self.zs[i] < z_min {
                    z_min = self.zs[i];
                } else if self.zs[i] > z_max {
                    z_max = self.zs[i];
                }
            }
            self.min_xs[node] = x_min;
            self.min_ys[node] = y_min;
            self.min_zs[node] = z_min;
            self.max_xs[node] = x_max;
            self.max_ys[node] = y_max;
            self.max_zs[node] = z_max;
            self.leaf_start[node] = from as i32;
            self.leaf_end[node] = to as i32;
            return node + 1;
        }

        let x_width = p_max_x - p_min_x;
        let y_width = p_max_y - p_min_y;
        let z_width = p_max_z - p_min_z;
        let dim: u8 = if x_width >= y_width && x_width >= z_width {
            0
        } else if y_width >= z_width {
            1
        } else {
            2
        };

        let mid = from + count / 2;
        self.quickselect(from, to - 1, mid, dim);
        let split_val = self.key(dim, mid);
        self.split_dims[node] = dim;
        self.split_vals[node] = split_val;

        let (mut l_max_x, mut l_max_y, mut l_max_z) = (p_max_x, p_max_y, p_max_z);
        let (mut r_min_x, mut r_min_y, mut r_min_z) = (p_min_x, p_min_y, p_min_z);
        match dim {
            0 => {
                l_max_x = split_val;
                r_min_x = split_val;
            }
            1 => {
                l_max_y = split_val;
                r_min_y = split_val;
            }
            _ => {
                l_max_z = split_val;
                r_min_z = split_val;
            }
        }

        let left_count = self.build_node(2 * node + 1, from, mid, p_min_x, p_min_y, p_min_z, l_max_x, l_max_y, l_max_z);
        let right_count = self.build_node(2 * node + 2, mid, to, r_min_x, r_min_y, r_min_z, p_max_x, p_max_y, p_max_z);

        let l = 2 * node + 1;
        let r = 2 * node + 2;
        self.min_xs[node] = jmin(self.min_xs[l], self.min_xs[r]);
        self.min_ys[node] = jmin(self.min_ys[l], self.min_ys[r]);
        self.min_zs[node] = jmin(self.min_zs[l], self.min_zs[r]);
        self.max_xs[node] = jmax(self.max_xs[l], self.max_xs[r]);
        self.max_ys[node] = jmax(self.max_ys[l], self.max_ys[r]);
        self.max_zs[node] = jmax(self.max_zs[l], self.max_zs[r]);

        (node + 1).max(left_count.max(right_count))
    }
}

impl KdTree3D {
    /// Build from points; payload of each point is its index in `pts`.
    pub fn build(pts: &[[f64; 3]]) -> KdTree3D {
        let n = pts.len();
        let mut xs: Vec<f64> = pts.iter().map(|p| p[0]).collect();
        let mut ys: Vec<f64> = pts.iter().map(|p| p[1]).collect();
        let mut zs: Vec<f64> = pts.iter().map(|p| p[2]).collect();
        let mut ids: Vec<u32> = (0..n as u32).collect();
        if n == 0 {
            return KdTree3D {
                xs,
                ys,
                zs,
                ids,
                split_dims: vec![],
                split_vals: vec![],
                min_xs: vec![],
                min_ys: vec![],
                min_zs: vec![],
                max_xs: vec![],
                max_ys: vec![],
                max_zs: vec![],
                leaf_start: vec![],
                leaf_end: vec![],
                node_count: 0,
            };
        }
        let mut max_depth = 0;
        let mut tmp = n;
        while tmp > BUCKET_SIZE {
            max_depth += 1;
            tmp = (tmp + 1) / 2;
        }
        let mut max_nodes = (1usize << (max_depth + 1)) - 1;
        if max_nodes < 1 {
            max_nodes = 1;
        }

        let (mut x_min, mut y_min, mut z_min) = (xs[0], ys[0], zs[0]);
        let (mut x_max, mut y_max, mut z_max) = (xs[0], ys[0], zs[0]);
        for i in 1..n {
            if xs[i] < x_min {
                x_min = xs[i];
            } else if xs[i] > x_max {
                x_max = xs[i];
            }
            if ys[i] < y_min {
                y_min = ys[i];
            } else if ys[i] > y_max {
                y_max = ys[i];
            }
            if zs[i] < z_min {
                z_min = zs[i];
            } else if zs[i] > z_max {
                z_max = zs[i];
            }
        }

        let mut b = Builder {
            xs: &mut xs,
            ys: &mut ys,
            zs: &mut zs,
            ids: &mut ids,
            split_dims: vec![0; max_nodes],
            split_vals: vec![0.0; max_nodes],
            min_xs: vec![0.0; max_nodes],
            min_ys: vec![0.0; max_nodes],
            min_zs: vec![0.0; max_nodes],
            max_xs: vec![0.0; max_nodes],
            max_ys: vec![0.0; max_nodes],
            max_zs: vec![0.0; max_nodes],
            leaf_start: vec![-1; max_nodes],
            leaf_end: vec![0; max_nodes],
        };
        let node_count = b.build_node(0, 0, n, x_min, y_min, z_min, x_max, y_max, z_max);
        let Builder { split_dims, split_vals, min_xs, min_ys, min_zs, max_xs, max_ys, max_zs, leaf_start, leaf_end, .. } = b;
        KdTree3D {
            xs,
            ys,
            zs,
            ids,
            split_dims,
            split_vals,
            min_xs,
            min_ys,
            min_zs,
            max_xs,
            max_ys,
            max_zs,
            leaf_start,
            leaf_end,
            node_count,
        }
    }

    pub fn size(&self) -> usize {
        self.xs.len()
    }

    #[inline(always)]
    fn box_dist(&self, node: usize, qx: f64, qy: f64, qz: f64) -> f64 {
        sqr_dist_to_box(
            qx,
            qy,
            qz,
            self.min_xs[node],
            self.min_ys[node],
            self.min_zs[node],
            self.max_xs[node],
            self.max_ys[node],
            self.max_zs[node],
        )
    }

    /// Squared distance to the nearest point (+inf if empty)
    pub fn nearest_sqr_dist(&self, qx: f64, qy: f64, qz: f64) -> f64 {
        if self.size() == 0 {
            return f64::INFINITY;
        }
        let mut stack = [0usize; 64];
        let mut sp = 0;
        stack[sp] = 0;
        sp += 1;
        let mut best = f64::INFINITY;
        while sp > 0 {
            sp -= 1;
            let node = stack[sp];
            if node >= self.node_count {
                continue;
            }
            if self.box_dist(node, qx, qy, qz) >= best {
                continue;
            }
            if self.leaf_start[node] >= 0 {
                for i in self.leaf_start[node] as usize..self.leaf_end[node] as usize {
                    let d = sqr_dist(qx, qy, qz, self.xs[i], self.ys[i], self.zs[i]);
                    if d < best {
                        best = d;
                    }
                }
                continue;
            }
            let left = 2 * node + 1;
            let right = 2 * node + 2;
            let qv = match self.split_dims[node] {
                0 => qx,
                1 => qy,
                _ => qz,
            };
            if qv <= self.split_vals[node] {
                if right < self.node_count {
                    stack[sp] = right;
                    sp += 1;
                }
                if left < self.node_count {
                    stack[sp] = left;
                    sp += 1;
                }
            } else {
                if left < self.node_count {
                    stack[sp] = left;
                    sp += 1;
                }
                if right < self.node_count {
                    stack[sp] = right;
                    sp += 1;
                }
            }
        }
        best
    }

    /// Visit payloads of all points within squared radius, in KdTree3D.findWithinRadius order.
    #[inline]
    pub fn for_each_within(&self, qx: f64, qy: f64, qz: f64, sqr_radius: f64, mut f: impl FnMut(u32)) {
        if self.size() == 0 {
            return;
        }
        let mut stack = [0usize; 64];
        let mut sp = 0;
        stack[sp] = 0;
        sp += 1;
        while sp > 0 {
            sp -= 1;
            let node = stack[sp];
            if node >= self.node_count {
                continue;
            }
            if self.box_dist(node, qx, qy, qz) > sqr_radius {
                continue;
            }
            if self.leaf_start[node] >= 0 {
                for i in self.leaf_start[node] as usize..self.leaf_end[node] as usize {
                    let d = sqr_dist(qx, qy, qz, self.xs[i], self.ys[i], self.zs[i]);
                    if d <= sqr_radius {
                        f(self.ids[i]);
                    }
                }
                continue;
            }
            let left = 2 * node + 1;
            let right = 2 * node + 2;
            if right < self.node_count {
                stack[sp] = right;
                sp += 1;
            }
            if left < self.node_count {
                stack[sp] = left;
                sp += 1;
            }
        }
    }

    pub fn find_within(&self, qx: f64, qy: f64, qz: f64, sqr_radius: f64) -> Vec<u32> {
        let mut res = Vec::with_capacity(64);
        self.for_each_within(qx, qy, qz, sqr_radius, |i| res.push(i));
        res
    }

    pub fn count_within(&self, qx: f64, qy: f64, qz: f64, sqr_radius: f64) -> usize {
        let mut c = 0;
        self.for_each_within(qx, qy, qz, sqr_radius, |_| c += 1);
        c
    }
}
