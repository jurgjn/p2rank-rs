//! Random forest inference, port of FasterForest's `LegacyFlatBinaryForest`
//! (sum class probabilities over trees, then normalize).
//!
//! The forest can be loaded from P2Rank's `model.zst` (see `javaser`) or from the `model.flb`
//! cache file, which stores the in-memory node layout and is memory-mapped without copying.

use std::path::Path;

/// Internal node: split point, children (negative = leaf index) and attribute index.
#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct Node {
    split: f64,
    left: i32,
    right: i32,
    attr: u32,
    _pad: u32,
}

enum Storage {
    Owned { nodes: Vec<Node>, leaves: Vec<[f64; 2]> },
    Mapped { map: memmap2::Mmap, nodes_off: usize, n_nodes: usize, leaves_off: usize, n_leaves: usize },
}

pub struct Forest {
    pub num_trees: usize,
    pub num_attributes: usize,
    storage: Storage,
}

const MAGIC: &[u8; 8] = b"P2RFLB02";
const HEADER_LEN: usize = 32;

impl Forest {
    pub fn from_parts(
        num_trees: usize,
        num_attributes: usize,
        left: Vec<i32>,
        right: Vec<i32>,
        attr: Vec<i32>,
        split: Vec<f64>,
        leaves: Vec<[f64; 2]>,
    ) -> Result<Forest, String> {
        let nodes =
            (0..left.len()).map(|i| Node { split: split[i], left: left[i], right: right[i], attr: attr[i] as u32, _pad: 0 }).collect();
        let f = Forest { num_trees, num_attributes, storage: Storage::Owned { nodes, leaves } };
        f.validate()?;
        Ok(f)
    }

    #[inline(always)]
    fn nodes(&self) -> &[Node] {
        match &self.storage {
            Storage::Owned { nodes, .. } => nodes,
            Storage::Mapped { map, nodes_off, n_nodes, .. } => unsafe {
                std::slice::from_raw_parts(map.as_ptr().add(*nodes_off) as *const Node, *n_nodes)
            },
        }
    }

    #[inline(always)]
    fn leaves(&self) -> &[[f64; 2]] {
        match &self.storage {
            Storage::Owned { leaves, .. } => leaves,
            Storage::Mapped { map, leaves_off, n_leaves, .. } => unsafe {
                std::slice::from_raw_parts(map.as_ptr().add(*leaves_off) as *const [f64; 2], *n_leaves)
            },
        }
    }

    /// Checks that all child / leaf / attribute references are in range (so traversal can't go out of bounds).
    fn validate(&self) -> Result<(), String> {
        let nodes = self.nodes();
        let nl = self.leaves().len() as i64;
        let n = nodes.len() as i64;
        if self.num_trees as i64 > n {
            return Err("invalid forest: more trees than nodes".into());
        }
        for nd in nodes {
            for c in [nd.left as i64, nd.right as i64] {
                if (c >= 0 && c >= n) || (c < 0 && -c >= nl) {
                    return Err("invalid forest: child index out of range".into());
                }
            }
            if nd.attr as usize >= self.num_attributes.max(1) {
                return Err("invalid forest: attribute index out of range".into());
            }
        }
        Ok(())
    }

    /// Re-indexes the split attributes: `map[i]` is the feature vector index of model attribute `i`
    /// (for models trained on a subset of the computed features, e.g. `alphafold` without `bfactor`).
    pub fn remap_attributes(&mut self, map: &[usize], num_attributes: usize) -> Result<(), String> {
        let Storage::Owned { nodes, .. } = &mut self.storage else {
            return Err("cannot remap the attributes of a memory-mapped model".into());
        };
        for nd in nodes.iter_mut() {
            nd.attr = *map.get(nd.attr as usize).ok_or("invalid forest: attribute index out of range")? as u32;
        }
        self.num_attributes = num_attributes;
        self.validate()
    }

    /// Loads the `model.flb` cache (memory mapped).
    pub fn load_flb(path: &Path) -> Result<Forest, String> {
        let file = std::fs::File::open(path).map_err(|e| format!("cannot read model {}: {}", path.display(), e))?;
        let map = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| format!("cannot map model {}: {}", path.display(), e))?;
        let b = &map[..];
        if b.len() < HEADER_LEN || &b[0..8] != MAGIC || !cfg!(target_endian = "little") {
            return Err(format!("invalid model cache {}", path.display()));
        }
        let rd = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
        let (num_trees, num_attributes, n_nodes, n_leaves) = (rd(8), rd(12), rd(16), rd(20));
        let nodes_off = HEADER_LEN;
        let leaves_off = nodes_off + n_nodes * std::mem::size_of::<Node>();
        if b.len() != leaves_off + n_leaves * 16 {
            return Err(format!("invalid model cache size {}", path.display()));
        }
        // fault the pages in on a background thread, overlapping with structure parsing and surface
        // computation (concurrent page faults during inference are much slower)
        {
            let addr = map.as_ptr() as usize;
            let len = map.len();
            std::thread::spawn(move || {
                let mut sum = 0u8;
                let mut off = 0;
                while off < len {
                    // SAFETY: the mapping lives as long as the process uses the forest (never unmapped
                    // before exit in prank-rs); reads of a read-only mapping
                    sum = sum.wrapping_add(unsafe { std::ptr::read_volatile((addr + off) as *const u8) });
                    off += 4096;
                }
                std::hint::black_box(sum);
            });
        }
        // no full validation pass here: traversal uses checked indexing
        Ok(Forest { num_trees, num_attributes, storage: Storage::Mapped { map, nodes_off, n_nodes, leaves_off, n_leaves } })
    }

    /// Writes the `model.flb` cache atomically (temp file + rename).
    pub fn write_flb(&self, path: &Path) -> std::io::Result<()> {
        use std::io::Write;
        let nodes = self.nodes();
        let leaves = self.leaves();
        let tmp = path.with_extension(format!("flb.tmp{}", std::process::id()));
        let res = (|| {
            let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
            let mut h = [0u8; HEADER_LEN];
            h[0..8].copy_from_slice(MAGIC);
            h[8..12].copy_from_slice(&(self.num_trees as u32).to_le_bytes());
            h[12..16].copy_from_slice(&(self.num_attributes as u32).to_le_bytes());
            h[16..20].copy_from_slice(&(nodes.len() as u32).to_le_bytes());
            h[20..24].copy_from_slice(&(leaves.len() as u32).to_le_bytes());
            w.write_all(&h)?;
            for nd in nodes {
                w.write_all(&nd.split.to_le_bytes())?;
                w.write_all(&nd.left.to_le_bytes())?;
                w.write_all(&nd.right.to_le_bytes())?;
                w.write_all(&nd.attr.to_le_bytes())?;
                w.write_all(&0u32.to_le_bytes())?;
            }
            for l in leaves {
                w.write_all(&l[0].to_le_bytes())?;
                w.write_all(&l[1].to_le_bytes())?;
            }
            w.flush()?;
            drop(w);
            std::fs::rename(&tmp, path)
        })();
        if res.is_err() {
            let _ = std::fs::remove_file(&tmp);
        }
        res
    }

    /// Accumulates the class sums of trees `trees` for the instances `idx` (tree-major; each lane
    /// walks one instance, refilled with the next pending instance when it reaches a leaf).
    fn accumulate(&self, xs: &[&[f64]], idx: &[usize], trees: std::ops::Range<usize>, s0: &mut [f64], s1: &mut [f64]) {
        const LANES: usize = 16;
        let n = idx.len();
        let nodes = self.nodes();
        let leaves = self.leaves();
        for t in trees {
            let mut lane_pt = [0usize; LANES];
            let mut lane_node = [t as i32; LANES];
            let mut active = LANES.min(n);
            for l in 0..active {
                lane_pt[l] = idx[l];
            }
            let mut next_pt = active;
            while active > 0 {
                let mut l = 0;
                while l < active {
                    let nd = &nodes[lane_node[l] as usize];
                    let x = xs[lane_pt[l]][nd.attr as usize];
                    let children = [nd.left, nd.right];
                    let next = children[(!(x < nd.split)) as usize];
                    if next < 0 {
                        let p = &leaves[(-next) as usize];
                        let i = lane_pt[l];
                        s0[i] += p[0];
                        s1[i] += p[1];
                        if next_pt < n {
                            lane_pt[l] = idx[next_pt];
                            lane_node[l] = t as i32;
                            next_pt += 1;
                        } else {
                            active -= 1;
                            lane_pt[l] = lane_pt[active];
                            lane_node[l] = lane_node[active];
                            continue;
                        }
                    } else {
                        lane_node[l] = next;
                    }
                    l += 1;
                }
            }
        }
    }

    /// LegacyFlatBinaryForest.predictForBatch over a block of instances, optionally with early exit (fast mode): after `k` trees, instances whose
    /// partial score is below `cutoff` keep that partial score; `stages` = [(k, cutoff), ...].
    /// Instances that pass all stages get exactly the full-forest score.
    /// Returns, for each instance, the number of trees used (num_trees if not stopped early).
    pub fn predict_block_cascade(&self, xs: &[&[f64]], stages: &[(usize, f64)], out: &mut [f64]) -> Vec<u16> {
        self.predict_block_cascade_k(xs, stages, self.num_trees, out)
    }

    /// Like [`Self::predict_block_cascade`], with only the first `trees` trees making up the forest
    /// (random forest trees are i.i.d., so this is an unbiased smaller forest).
    pub fn predict_block_cascade_k(&self, xs: &[&[f64]], stages: &[(usize, f64)], trees: usize, out: &mut [f64]) -> Vec<u16> {
        let n = xs.len();
        let mut s0 = vec![0.0f64; n];
        let mut s1 = vec![0.0f64; n];
        let used = self.advance(xs, &mut s0, &mut s1, 0, stages, trees);
        for i in 0..n {
            out[i] = Self::normalized(s0[i], s1[i]);
        }
        used
    }

    /// NormalizationUtils.normalizedProb of the class sums
    pub fn normalized(s0: f64, s1: f64) -> f64 {
        if s0 + s1 == 0.0 {
            0.0
        } else {
            s1 / (s0 + s1)
        }
    }

    /// Resumable cascade: the instances have class sums `s0`/`s1` over the first `from` trees;
    /// adds trees up to `trees`, stopping instances early at the `stages` (k, cutoff) with k > from.
    /// Trees are added one by one in order, so resuming gives the same sums as one pass.
    /// Returns, for each instance, the number of trees in its sums.
    pub fn advance(&self, xs: &[&[f64]], s0: &mut [f64], s1: &mut [f64], from: usize, stages: &[(usize, f64)], trees: usize) -> Vec<u16> {
        let trees = trees.clamp(1, self.num_trees);
        let n = xs.len();
        for x in xs {
            assert!(x.len() + 1 >= self.num_attributes, "feature vector too short for model");
        }
        let mut alive: Vec<usize> = (0..n).collect();
        let mut used = vec![trees as u16; n];
        let mut done_trees = from.min(trees);
        for &(k, cutoff) in stages {
            let k = k.min(trees);
            if k <= done_trees {
                continue;
            }
            self.accumulate(xs, &alive, done_trees..k, s0, s1);
            done_trees = k;
            alive.retain(|&i| {
                if Self::normalized(s0[i], s1[i]) < cutoff {
                    used[i] = k as u16;
                    false
                } else {
                    true
                }
            });
        }
        self.accumulate(xs, &alive, done_trees..trees, s0, s1);
        used
    }
}
