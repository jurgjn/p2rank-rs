#!/usr/bin/env python3
"""Evaluate prank-rs predictions with P2Rank's own (Java) evaluation code.

Runs `prank-rs predict` on the proteins of a ligand-annotated dataset (e.g. coach420(mlig).ds),
then scores the predicted pockets with `prank eval-rescore -rescorer IdentityRescorer` (which keeps
the pockets and their order). P2Rank's loader for its own *_predictions.csv reads pocket centers
from the wrong columns, so an adapter CSV with the layout it expects is written first. Only
center-based metrics (DCA, DCC) are meaningful this way.

usage: eval_predictions.py <dataset.ds> <workdir> [prank-rs predict args, e.g. -mode fast]
prints: the metrics line from the Java evaluation (stats.csv)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
P2 = os.path.abspath(os.path.join(HERE, "..", ".."))
METRICS = ["DCA_4_0", "DCA_4_2", "DCC_10_0", "DCC_10_2", "AVG_POCKETS", "LIGANDS"]


def main():
    ds, work, extra = sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3:]
    ds = os.path.abspath(ds)
    ds_dir = os.path.dirname(ds)
    header = ["protein"]
    items = []
    for line in open(ds):
        line = line.split("#")[0].strip()
        if not line or line.startswith("PARAM."):
            continue
        if line.startswith("HEADER:"):
            header = line.split()[1:]
            continue
        cols = line.split()
        items.append(dict(zip(header, cols)))
    os.makedirs(work, exist_ok=True)
    pred_dir = os.path.join(work, "pred")
    prot_ds = os.path.join(work, "proteins.ds")
    with open(prot_ds, "w") as f:
        for it in items:
            f.write(os.path.join(ds_dir, it["protein"]) + "\n")
    env = dict(os.environ, P2RANK_HOME=os.path.join(P2, "distro"))
    subprocess.run([os.path.join(P2, "rust", "target", "release", "prank-rs"), "predict", prot_ds, "-o", pred_dir,
                    "-visualizations", "0"] + extra, env=env, check=True, stdout=subprocess.DEVNULL)

    adapter_dir = os.path.join(work, "adapter")
    os.makedirs(adapter_dir, exist_ok=True)
    eval_ds = os.path.join(work, "eval.ds")
    with open(eval_ds, "w") as f:
        f.write("PARAM.PREDICTION_METHOD=p2rank\n\n")
        cols = ["prediction", "protein"] + [h for h in header if h not in ("protein", "prediction")]
        f.write("HEADER: " + " ".join(cols) + "\n")
        for it in items:
            label = os.path.basename(it["protein"])
            src = os.path.join(pred_dir, label + "_predictions.csv")
            dst = os.path.join(adapter_dir, label + "_predictions.csv")
            with open(src) as fi, open(dst, "w") as fo:
                rows = [r.split(",") for r in fi.read().splitlines()]
                fo.write("name,rank,score,x,x,center_x,center_y,center_z\n")
                for r in rows[1:]:
                    r = [c.strip() for c in r]
                    fo.write(",".join([r[0], r[1], r[2], "0", "0", r[6], r[7], r[8]]) + "\n")
            vals = [dst, os.path.join(ds_dir, it["protein"])] + [it[h] for h in cols[2:]]
            f.write("  ".join(vals) + "\n")

    out = os.path.join(work, "eval")
    subprocess.run([os.path.join(P2, "distro", "prank"), "eval-rescore", eval_ds, "-rescorer", "IdentityRescorer",
                    "-o", out, "-visualizations", "0"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    stats = {}
    for line in open(os.path.join(out, "stats.csv")):
        k, _, v = line.partition(",")
        stats[k.strip()] = v.strip()
    print(" ".join(f"{m}={stats.get(m, '?')}" for m in METRICS))


if __name__ == "__main__":
    main()
