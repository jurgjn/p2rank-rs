# p2rank-rs

> (Human-written text starts)
>
> This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (https://github.com/evocellnet/ai-hackathon/). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.
>
> The goal of this specific project was to port the P2Rank ligand-binding site prediction tool (`predict` command) from Java to Rust, with output byte-identical to the original, and then to add faster modes that trade small, documented amounts of accuracy for speed. It was executed by Claude Opus 5.5 with delegation to Opus 5.5 subagents.
>
> Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.
>
> This repository is a fork of P2Rank: the Java P2Rank code is the original, human-written project. The Rust port in `rust/` and this README (a benchmark and accuracy report) are entirely AI generated.
>
> (Human-written text ends)

`prank-rs` (in [`rust/`](rust/)) is a Rust port of the P2Rank `predict` command for the default
config and model. This page compares the speed of the unchanged Java P2Rank with the three
`prank-rs` modes and summarizes their accuracy. For the full documentation see
[`rust/README.md`](rust/README.md); the original P2Rank README is
[`README-P2Rank.md`](README-P2Rank.md).

| Mode | What it does |
|---|---|
| `exact` (default) | Output files byte-identical to Java P2Rank |
| `fast` | Random forest early exit for background points; pockets identical to `exact` on all 4429 coach420 + holo4k structures |
| `faster` | `fast` plus scoring only part of the protein surface; occasionally drops a low-ranked pocket |

## Install from GitHub

Requirements: a Rust toolchain (`cargo`, tested with 1.99.1) and a C compiler (zlib is built
statically). Java is **not** needed: the model files are part of the repository.

```bash
git clone --depth 1 https://github.com/jurgjn/p2rank-rs.git
cd p2rank-rs/rust
cargo build --release
target/release/prank-rs predict -f ../distro/test_data/1fbl.pdb -o out/
```

Run on a dataset file (`.ds`, as in P2Rank), selecting a mode:

```bash
target/release/prank-rs predict ../distro/test_data/fpocket.ds -mode faster -o out_ds/
```

> [!NOTE]
> `prank-rs` reads the model and score transformers from the `distro/` directory of the checkout.
> It finds it by looking next to and above the executable, so the binary works in place. To copy
> the binary elsewhere, set `P2RANK_HOME` to the checkout's `distro/` directory.

> [!TIP]
> The first run converts `distro/models/default/model.zst` into a memory-mapped cache
> (`model.flb`, loads in under 1 ms). This needs write access to `distro/models/default/`; run
> `prank-rs convert-model` once as the owner for a read-only installation.

## Speed

Measured on 2026-10-28 at commit `acc577ae`: median of 11 runs (single structures) or 5 runs
(datasets), output to local disk, default threads. Java P2Rank 2.6-alpha.8 on OpenJDK 25,
run offline (like `prank-rs -ccd reduced`, the default).

Wall time:

| | Java P2Rank | `exact` | `fast` | `faster` |
|---|---:|---:|---:|---:|
| `1fbl.pdb` (2.9k atoms) | 3.45 s | 0.079 s | 0.078 s | 0.077 s |
| `1AHP.pdb` (12.5k atoms) | 3.84 s | 0.141 s | 0.127 s | 0.116 s |
| coach420 (420 structures) | 8.9 s | 1.28 s | 0.89 s | 0.87 s |
| holo4k (4009 structures) | 30.0 s | 10.6 s | 6.2 s | 5.4 s |

CPU time (user + system, all threads):

| | Java P2Rank | `exact` | `fast` | `faster` |
|---|---:|---:|---:|---:|
| `1fbl.pdb` | 14.4 s | 0.4 s | 0.3 s | 0.2 s |
| `1AHP.pdb` | 16.6 s | 0.7 s | 0.4 s | 0.4 s |
| coach420 | 359 s | 39.8 s | 22.4 s | 19.6 s |
| holo4k | 1413 s | 548 s | 294 s | 244 s |

Speedup over Java P2Rank:

| | `exact` | `fast` | `faster` |
|---|---:|---:|---:|
| single structure, wall | 27-44x | 30-44x | 33-45x |
| coach420, wall / CPU | 7.0x / 9.0x | 10.0x / 16x | 10.3x / 18x |
| holo4k, wall / CPU | 2.8x / 2.6x | 4.9x / 4.8x | 5.5x / 5.8x |

> [!NOTE]
> For single structures most of the Java time is JVM start-up, Groovy config compilation and
> model loading, which `prank-rs` avoids; all three modes take under 0.15 s. On datasets the
> differences come from the per-structure work, where `fast` and `faster` skip most of the random
> forest evaluation for background surface points.

> [!WARNING]
> The machine was a shared node (2x AMD EPYC 9554, 64 of 128 cores allocated) with a load average
> between 7 and 65 during the runs. Expect run-to-run variation of about 5%, which is about the
> size of the `fast` to `faster` difference on coach420.

## Accuracy

Ligand-site benchmarks with P2Rank's own evaluation code (`rust/tools/eval_predictions.py`,
coach420(mlig) and holo4k(mlig)). DCA_4_0 / DCA_4_2: fraction of ligands whose closest atom is
within 4 A of the center of one of the top n / n+2 predicted pockets (n = number of ligands of
the protein).

| | coach420 DCA_4_0 | coach420 DCA_4_2 | holo4k DCA_4_0 | holo4k DCA_4_2 |
|---|---:|---:|---:|---:|
| Java P2Rank = `exact` | 0.6969 | 0.7332 | 0.7671 | 0.8346 |
| `fast` | 0.6969 | 0.7332 | 0.7671 | 0.8346 |
| `faster` | 0.6969 | 0.7332 | 0.7671 | 0.8346 |
| seed-to-seed SD (5 retrained models) | 0.0045 | 0.0028 | 0.0006 | 0.0007 |

Agreement with `exact`:

| | `fast` | `faster` |
|---|---:|---:|
| identical `predictions.csv`, coach420 / holo4k | 100% / 100% | 96.3% / 92.6% |
| top-1 pocket identical | 100% | 100% |
| residue score, mean abs. difference (holo4k) | 0.016 | 0.044 |
| residue pocket assignment agreement (holo4k) | 100% | 99.8% |

* **`exact`**: byte-identical output files to Java P2Rank in all validation runs (PDB, mmCIF,
  BinaryCIF and AlphaFold models; see the Validation table in `rust/README.md`).
* **`fast`**: identical pockets. Only residue scores and point scores of background points far from
  any pocket differ, by much less than between two retrained models (residue score mean abs.
  difference 0.113).
* **`faster`**: every pocket it finds is identical to `exact` (score, rank order, center,
  residues), because pockets are computed only from full random forest scores. It can miss a
  pocket when no point near it is among the initially scored ones: 11 pockets in 300 coach420
  structures and 290 in 3448 holo4k structures. None of them was ranked first; their median
  P2Rank probability is 0.004 (maximum 0.17). It never added or moved a pocket.

> [!TIP]
> If a low-ranked pocket can matter for your use, or you need output identical to P2Rank, use
> `exact` (identical files) or `fast` (identical pockets). `faster` is meant for large-scale runs
> where the top pockets are what counts.
