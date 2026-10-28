# prank-rs: P2Rank `predict` in Rust

`prank-rs` is a Rust port of the P2Rank `predict` command (default config, default model).
In the default `exact` mode it produces **byte-identical output files** to the Java P2Rank, and is
25-40x faster for a single structure and 3-10x faster on whole datasets. The optional `fast` mode
is another 1.5-1.7x faster on datasets, with documented differences far below P2Rank's own
seed-to-seed variation (see [Modes](#modes-exact-and-fast)).

```bash
cd rust
cargo build --release
P2RANK_HOME=../distro target/release/prank-rs predict -f ../distro/test_data/1fbl.pdb -o out/
```

> [!NOTE]
> `prank-rs` uses the model and score transformer files of a P2Rank installation (the `distro`
> directory). It finds it via `P2RANK_HOME`, or by looking next to / above the executable.

## Usage

```
prank-rs predict -f <structure> [-o <outdir>] [-mode exact|fast] [-threads <n>] [-visualizations 0|1] [-ccd reduced|full]
prank-rs predict <dataset.ds>  [-o <outdir>] [-mode exact|fast] [-threads <n>] [-visualizations 0|1] [-ccd reduced|full]
prank-rs convert-model [<model_dir>]
```

* Input formats: PDB (`.pdb`, `.ent`), mmCIF (`.cif`) and BinaryCIF (`.bcif`), optionally
  compressed with gzip (`.gz`) or zstd (`.zst`).
* Dataset files (`.ds`) are read like P2Rank does: `#` comments, `PARAM.` lines, `HEADER:` columns.
  The `chains` column is supported; other columns (ligands, predictions) do not affect `predict`.
* The default output directory is `<P2RANK_HOME>/test_output/predict_<name>`, as in P2Rank.
* `-c default` and `-m default` are accepted; any other config, model or parameter is rejected
  rather than silently giving different results.

> [!TIP]
> The first run converts `models/default/model.zst` (zstd-compressed Java serialization) to
> `models/default/model.flb`, a cache that is memory-mapped on later runs and loads in under 1 ms.
> Run `prank-rs convert-model` to create it explicitly (e.g. for a read-only installation, as
> the user owning it). The cache is ignored if it is older than `model.zst`.

## Modes: `exact` and `fast`

Following the mode vocabulary of
[uplifting-biomolecular-modeling](https://github.com/anthropics/uplifting-biomolecular-modeling):

* `exact` (default): output byte-identical to Java P2Rank.
* `fast`: faster, with small documented numeric differences that stay within the tool's own
  seed-to-seed variation.

`fast` changes three things:

1. **Random forest early exit, with pocket regions rescored exactly.** About 92% of SAS points
   score far below the 0.35 ligandability threshold. Points whose partial score after 10 trees is
   below 0.05, or after 20 trees below 0.1, keep that partial score, so most points need 10-20 instead of 100 trees.
   Then every early-stopped point within 4 A of a point scoring at least 0.2 is rescored with all
   100 trees, repeated until stable. Pockets are therefore computed from exact scores; only
   background points far from any pocket keep approximate scores.
2. **Residue scores and feature neighborhoods are summed in a different order** (grid instead of
   kd-tree traversal order). This only changes the last bits of floating point sums.
3. **The points file is compressed with gzip level 1 instead of 6** (identical content).

> [!IMPORTANT]
> In `fast` mode the `predictions.csv` files (pocket ranks, scores, probabilities, centers,
> residues, atoms) were byte-identical to `exact` mode for all 4429 coach420 and holo4k
> structures. Differences are limited to the residue scores in `residues.csv` and the point
> scores in the points file, from background points far from any pocket.

### Accuracy

The yardstick is P2Rank's own seed-to-seed variation: the default model retrained
(`traineval -t chen11-fpocket.ds -c config/train-default.groovy`) with different random seeds.

Ligand-site benchmarks (P2Rank's evaluation code, via `tools/eval_predictions.py`):

| | coach420 DCA_4_0 | coach420 DCA_4_2 | holo4k DCA_4_0 | holo4k DCA_4_2 |
|---|---:|---:|---:|---:|
| `exact` (default model) | 0.6969 | 0.7332 | 0.7671 | 0.8346 |
| `fast` (default model) | 0.6969 | 0.7332 | 0.7671 | 0.8346 |
| seed-to-seed SD (5 retrained models) | 0.0045 | 0.0028 | 0.0006 | 0.0007 |

Agreement on holo4k (4009 structures, `tools/compare_modes.py`), compared with the difference
between two retrained models (seeds 42 and 43), both run in `exact` mode:

| | `fast` vs `exact` | seed 42 vs seed 43 |
|---|---:|---:|
| identical `predictions.csv` | 100% | 0.02% |
| top-1 pocket center within 1 A | 100% | 85.9% |
| pocket score, mean relative difference | 0 | 8.4% |
| residue score, mean abs. difference | 0.016 | 0.113 |
| residue score, Pearson r | 0.99997 | 0.9984 |
| residue probability, mean abs. difference | 0.0014 | 0.0133 |
| residue pocket assignment agreement | 100% | 87.2% |

> [!NOTE]
> `fast` residue scores are on average slightly lower (by 0.015), because early-stopped
> background points tend to underestimate their full score. This is of the same size as the
> average difference between two seed-retrained models, and 7x smaller than their per-residue
> differences. The same `fast` settings applied to a retrained model (seed 42) also gave
> byte-identical `predictions.csv` for all of holo4k.

### `-ccd`: chemical component dictionary

How BioJava classifies residues depends on its chemical component dictionary, and P2Rank's
results depend on which one it gets:

* `reduced` (default): BioJava's bundled fallback set. This is what Java P2Rank uses **without
  network access** (and with an empty `/tmp/chemcomp` cache).
* `full`: the complete wwPDB CCD. This is what Java P2Rank uses when it **can download**
  definitions from RCSB (or has them cached in `/tmp/chemcomp`).

> [!IMPORTANT]
> The two modes give different results for structures with nonstandard residues: 20 of the 4009
> holo4k structures differ. This is a property of the Java P2Rank (see `misc/dev/backlog.md`);
> `prank-rs` reproduces either behaviour exactly.

## What is identical

For every structure, these files are byte-identical to the Java output:

* `<label>_predictions.csv`, `<label>_residues.csv`, `params.txt`
* `visualizations/<label>_pymol.pml`, `visualizations/<label>_chimerax.cxc`
* `visualizations/data/<label>_points.pdb.gz` (same zlib stream), `visualizations/data/<label>`
* `errors.csv`, `errors_aggregated.csv` for common failures (missing files, parse errors)

> [!NOTE]
> Not identical: `run.log` (similar content, different timings) and `errors_full.txt.gz`
> (contains Java stack traces). Like in P2Rank, the `git head:` line in `params.txt` comes from
> `git rev-parse` in the current directory.

## Validation

All output files compared against Java P2Rank 2.6-alpha.8 (JDK 25):

| Dataset | Structures | Files compared | Different |
|---|---:|---:|---:|
| `distro/test_data` (all PDB files) | 70 | 423 | 0 |
| coach420 | 420 | 2521 | 0 |
| holo4k | 4009 | 24055 | 0 |
| holo4k, Java with RCSB access vs `-ccd full` | 4009 | 24056 | 0 |
| chen11 | 251 | 1507 | 0 |
| fptrain | 222 | 1333 | 0 |
| ah4h (`.pdb.gz`) | 380 | 2281 | 0 |
| RCSB mmCIF (`.cif.gz`, random holo4k entries) | 299 | 1795 | 0 |
| RCSB BinaryCIF (`.bcif`) | 100 | 601 | 0 |
| AlphaFold DB models (`.pdb`, `.cif`, `.bcif`) | 72 | 433 | 0 |
| `test_data` CIF / BinaryCIF / zstd variants | 9 | 56 | 0 |

`cargo test` runs a regression subset (PDB, multi-chain, unusual elements, mmCIF, BinaryCIF)
against stored Java outputs in `tests/cases/`, checks `fast` mode against them within tolerance,
and checks number formatting against Java on 10k reference values. `tools/compare.sh` runs both implementations on a dataset and compares
all files.

## Performance

Wall time (median of several runs, local disk, shared 64-core node):

| | Java P2Rank | prank-rs `exact` | prank-rs `fast` |
|---|---:|---:|---:|
| `1fbl.pdb` (2.9k atoms) | 3.18 s | 0.083 s | 0.076 s |
| `1AHP.pdb` (12.5k atoms, 17k SAS points) | 3.71 s | 0.150 s | 0.124 s |
| coach420 (420 structures) | 9.2-11.8 s | 1.2 s | 0.82 s |
| holo4k (4009 structures) | 29-33 s | 9.5 s | 5.5 s |

CPU time on holo4k: Java 1324-1520 s, `exact` 503 s, `fast` 269 s. Peak memory is 2-5x lower
than Java. Set `P2RS_TIMING=1` to print a per-phase breakdown (`tools/phase_bench.py` reports the
minimum over several runs).

## How it works

The port follows the Java code path step by step, because floating point results depend on the
order of operations:

| Module | Ports |
|---|---|
| `pdb.rs`, `cif.rs`, `bcif.rs`, `chemcomp.rs` | BioJava `PDBFileParser`, ciftools-java tokenizer and BinaryCIF codecs, `CifStructureConsumerImpl`, entity assignment, alt locs, chem comp based group types |
| `protein.rs` | P2Rank residue chains, `AlternateChainReducer`, protein heavy atoms |
| `surface.rs` | faster-molecular-surface `DistinctPackedNumericalSurfaceV4` (tessellation and VdW radii tables dumped from CDK) |
| `kdtree.rs`, `atoms.rs` | P2Rank `KdTree3D` and `Atoms` (`cutoutSphere` / `cutoutShell`). The kd-tree traversal order defines neighbour order and therefore feature sums |
| `features.rs` | `chem`, `volsite`, `protrusion`, `bfactor`, `atom_table` features and their distance-weighted aggregation |
| `forest.rs`, `javaser.rs` | FasterForest `LegacyFlatBinaryForest`, read directly from the Java-serialized `model.zst` |
| `predict.rs` | clustering, pocket scoring, score transformers, residue scores |
| `output.rs`, `jfmt.rs`, `jdtoa.rs` | output files, including Java `DecimalFormat` (legacy `FloatingDecimal` digits), `Formatter` and Groovy `BigDecimal` division semantics |

Speed comes from avoiding JVM start-up, runtime Groovy compilation and Java deserialization,
plus optimizations that keep results exact: a memory-mapped model, tree-major random forest
inference with interleaved traversals, grid-based neighbor and shell tests where only a count or
an "any within distance" answer is needed, sparse feature aggregation (skipping zero terms
does not change IEEE sums here), and parallel surface, shell and output computation.

The constant tables (`src/data/`) were generated with `tools/gen_chemcomp.py` and the Java
programs in `tools/java/` (run against the P2Rank classpath).

## Limitations

* Only `predict` with the default config and model (`models/default`). No `rescore`,
  `eval-*`, training, conservation or AlphaFold-specific configs.
* PDB files with more than 500,000 atoms (where BioJava switches to CA-only parsing) and
  `.lzma` input are rejected.
* Number formatting assumes the `en_US` locale, which P2Rank forces at start-up.
