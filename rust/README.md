# prank-rs: P2Rank `predict` in Rust

`prank-rs` is a Rust port of the P2Rank `predict` command (default config, default model).
In the default `exact` mode it produces **byte-identical output files** to the Java P2Rank, and is
25-40x faster for a single structure and 3-10x faster on whole datasets. The optional `fast` mode
is another 1.5-1.7x faster on datasets, with documented differences far below P2Rank's own
seed-to-seed variation, and `faster` saves another 15-20% of CPU time by scoring only part of the
protein surface (see [Modes](#modes-exact-fast-and-faster)).

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
prank-rs predict -f <structure> [-o <outdir>] [-mode exact|fast|faster] [-threads <n>] [-visualizations 0|1] [-descriptors 0|1] [-ccd reduced|full]
prank-rs predict <dataset.ds>  [-o <outdir>] [-mode exact|fast|faster] [-threads <n>] [-visualizations 0|1] [-descriptors 0|1] [-ccd reduced|full]
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

## Modes: `exact`, `fast` and `faster`

Following the mode vocabulary of
[uplifting-biomolecular-modeling](https://github.com/anthropics/uplifting-biomolecular-modeling):

* `exact` (default): output byte-identical to Java P2Rank.
* `fast`: faster, with small documented numeric differences that stay within the tool's own
  seed-to-seed variation.
* `faster`: `fast` plus scoring only part of the surface. It occasionally drops a low-ranked
  pocket (see [`faster`](#faster-score-only-part-of-the-surface)).

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

### `faster`: score only part of the surface

`faster` builds on `fast` (same early exit, summation order and gzip level). Instead of computing
features and forest scores for every SAS point, it computes them on demand:

1. **Seeds.** Only points on 12 of the 42 tessellation directions are scored at first (the
   level-1 icosahedron directions are an exact subset of the level-2 ones used by P2Rank), about
   30% of all points, using the `fast` early exit.
2. **Expansion.** Every point within 4 A of a seed scoring at least 0.1, or of any point scoring
   at least 0.2, is scored with the same early exit, repeated until no new points are added.
3. **Completion.** Every point within 3.5 A (the extended pocket distance) of a point scoring at
   least 0.35 gets the full 100-tree score. The forest keeps per-point class sums, so completing a
   point adds the remaining trees instead of starting over, with bit-identical results.
4. **Background.** Early-stopped points with a partial score below 0.03 get 0.03, the RMS full
   score of such points on ah4h (a dataset used neither for training nor for this evaluation;
   their partial score is usually exactly 0). Points never scored get the distance-weighted mean
   of scored points within 3 A, or 0.03.

Clustering and pocket scores therefore only see full 100-tree scores, so every pocket `faster`
finds has the same score, rank order, center and residues as in `exact`. A pocket is lost only
if none of its points is near a seed scoring at least 0.1.

> [!WARNING]
> `faster` lost 11 of the pockets of 300 coach420 structures and 290 of 3448 holo4k structures.
> None was ranked first. The lost pockets have a median P2Rank probability of 0.004 (maximum
> 0.17). It never added or moved a pocket.

Same evaluation as above (coach420(mlig) and holo4k(mlig) proteins, compared with `exact`):

| | coach420 | holo4k |
|---|---:|---:|
| DCA_4_0 / DCA_4_2, `exact` | 0.6969 / 0.7332 | 0.7671 / 0.8346 |
| DCA_4_0 / DCA_4_2, `faster` | 0.6969 / 0.7332 | 0.7671 / 0.8346 |
| identical `predictions.csv` | 96.3% | 92.6% |
| top-1 pocket identical | 100% | 100% |
| residue score, mean abs. difference (bias) | 0.044 (-0.007) | 0.044 (-0.012) |
| residue pocket assignment agreement | 99.9% | 99.8% |

`faster` saves 17% of the CPU time of `fast` on holo4k and 13% on coach420 (see
[Performance](#performance)). Most of the remaining forest work is the full scoring of the
pocket regions.

> [!NOTE]
> Using only the first 70 of the 100 trees (dev override `P2RS_DEV_TREES=70`, i.e. an unbiased
> smaller random forest) saves another 10% of CPU time. On holo4k it keeps every agreement
> metric within seed-to-seed variation (pocket score mean relative difference 4.3%, top-1 center
> within 1 A 93.8%, residue score mean abs. difference 0.072), but it changes DCA_4_0 / DCA_4_2 to
> 0.7678 / 0.8381, more than the seed-to-seed SD (upwards, from more, smaller pockets). It is
> therefore not the default.

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

## Pocket descriptors (`-descriptors 1`)

With `-descriptors 1`, `prank-rs` also writes `<label>_pockets_descriptors.csv`: one row per
pocket (same names and ranks as `predictions.csv`) with cheap global properties of the kind used
by pocket fingerprints and druggability models. All other output files stay the same (byte-identical
to Java in `exact` mode), and the extra CPU time is under 1%.

| Column | Meaning | Inspired by |
|---|---|---|
| `residues` | residues lining the pocket | fpocket |
| `rgyr` | radius of gyration of the pocket SAS points (A) | DoGSiteScorer (size) |
| `axis1`, `axis2`, `axis3` | spread along the principal axes of the pocket points (A, largest first) | DoGSiteScorer (ellipsoid shape) |
| `shape` | `tunnel` (axis1/axis2 >= 1.8), `saucer` (axis2/axis3 >= 1.8) or `bowl` | |
| `buriedness` | mean number of buried protein atoms within 10 A of the pocket points (P2Rank's protrusion feature) | buriedness / depth |
| `aromatic`, `cation`, `anion`, `hydrophobic`, `acceptor`, `donor` | fraction of the pocket surface atoms with each VolSite pharmacophore property | FuzCav, VolSite |
| `polar_atoms` | fraction of N and O surface atoms | fpocket |
| `hydropathy` | mean Kyte-Doolittle hydropathy of the pocket residues | fpocket, DoGSiteScorer |
| `net_charge` | ARG + LYS - ASP - GLU over the pocket residues | fpocket |
| `bfactor` | mean B-factor of the surface atoms (flexibility; pLDDT for AlphaFold models) | |
| `score_max`, `score_mean` | max and mean ligandability score of the pocket points scoring >= 0.35 | |
| `flavour` | the pharmacophore property most enriched over an average pocket (z-score >= 1): `aromatic`, `basic`, `acidic`, `greasy`, `acceptor-rich`, `donor-rich`, else `balanced` | |

The `shape` thresholds (upper quartile of each axis ratio) and the means and SDs behind `flavour`
were computed on the 3575 pockets predicted on ah4h, a dataset used neither for training nor for
evaluation. There, about half the pockets are bowls and a quarter tunnels; a third are `balanced`.

> [!NOTE]
> These are descriptors, not predictions: none of them feeds back into pocket scores or ranks. On
> ah4h, the P2Rank pocket probability tracks size (Spearman 0.85-0.87 with `residues` and `rgyr`)
> and somewhat buriedness (0.26), but hardly the pharmacophore composition (|rho| <= 0.12).

Example (`1fbl`, a matrix metalloproteinase; pocket 1 is the S1' pocket next to the catalytic zinc):

```
name, rank, residues, rgyr, axis1, axis2, axis3, shape, buriedness, aromatic, cation, anion, hydrophobic, acceptor, donor, polar_atoms, hydropathy, net_charge, bfactor, score_max, score_mean, flavour
pocket1, 1, 18, 5.42, 4.79, 1.96, 1.60, tunnel, 136.4, 0.175, 0.000, 0.050, 0.325, 0.350, 0.300, 0.500, -1.13, 0, 21.3, 0.929, 0.613, donor-rich
```

References: fpocket ([Le Guilloux et al. 2009](https://doi.org/10.1186/1471-2105-10-168)),
DoGSiteScorer ([Volkamer et al. 2012](https://doi.org/10.1021/ci200454v)), FuzCav
([Weill and Rognan 2010](https://doi.org/10.1021/ci900349y)), PocketMatch
([Yeturu and Chandra 2008](https://doi.org/10.1186/1471-2105-9-543)), and the ProSPECCTs benchmark
of binding site comparison ([Ehrt et al. 2018](https://doi.org/10.1371/journal.pcbi.1006483)).

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
against stored Java outputs in `tests/cases/`, checks `fast` and `faster` modes against them within tolerance,
and checks number formatting against Java on 10k reference values. `tools/compare.sh` runs both implementations on a dataset and compares
all files.

## Performance

Wall time (median of several runs, local disk, shared 64-core node):

| | Java P2Rank | prank-rs `exact` | prank-rs `fast` | prank-rs `faster` |
|---|---:|---:|---:|---:|
| `1fbl.pdb` (2.9k atoms) | 3.18 s | 0.083 s | 0.076 s | 0.073 s |
| `1AHP.pdb` (12.5k atoms, 17k SAS points) | 3.71 s | 0.150 s | 0.124 s | 0.121 s |
| coach420 (420 structures) | 9.2-11.8 s | 1.2 s | 0.82 s | 0.80 s |
| holo4k (4009 structures) | 29-33 s | 9.5 s | 5.5 s | 5.0 s |

CPU time on holo4k: Java 1324-1520 s, `exact` 503 s, `fast` 269 s, `faster` 224-238 s. On
coach420: `exact` 36-38 s, `fast` 21.5 s, `faster` 18.5-18.8 s. Peak memory is 2-5x lower
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
