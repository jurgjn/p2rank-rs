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

    /// LegacyFlatBinaryForest.predictForBatch over a block of instances: tree-major traversal,
    /// several instances advanced in lockstep. Per-instance sums are accumulated in tree order,
    /// exactly as in Java.
    pub fn predict_block(&self, xs: &[&[f64]], out: &mut [f64]) {
        const LANES: usize = 16;
        let n = xs.len();
        let mut s0 = vec![0.0f64; n];
        let mut s1 = vec![0.0f64; n];
        let nodes = self.nodes();
        let leaves = self.leaves();
        for x in xs {
            assert!(x.len() + 1 >= self.num_attributes, "feature vector too short for model");
        }
        for t in 0..self.num_trees {
            // each lane walks one instance; when it reaches a leaf it continues with the next
            // pending instance, so all lanes stay busy (independent loads overlap)
            let mut lane_pt = [0usize; LANES];
            let mut lane_node = [t as i32; LANES];
            let mut active = LANES.min(n);
            for l in 0..active {
                lane_pt[l] = l;
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
                            lane_pt[l] = next_pt;
                            lane_node[l] = t as i32;
                            next_pt += 1;
                        } else {
                            // retire the lane: move the last active lane here
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
        for i in 0..n {
            // NormalizationUtils.normalizedProb
            let sum = s0[i] + s1[i];
            out[i] = if sum == 0.0 { 0.0 } else { s1[i] / sum };
        }
    }
}
