#!/usr/bin/env python3
"""Compare prank-rs outputs of two modes (e.g. exact vs fast) for the same structures.

usage: compare_modes.py <exact_outdir> <fast_outdir>

Reports agreement of pockets (count, top-1 center shift, rank-wise center match within 1 A / 4 A)
and residues (score error, probability error, pocket assignment agreement).
"""
import glob
import math
import os
import sys


def read_csv(path):
    lines = open(path).read().splitlines()
    hdr = [h.strip() for h in lines[0].split(",")]
    return [dict(zip(hdr, [c.strip() for c in l.split(",")])) for l in lines[1:]]


def main():
    a_dir, b_dir = sys.argv[1], sys.argv[2]
    n = 0
    same_count = 0
    top1 = []
    match1 = match4 = total_pockets = 0
    res_n = 0
    res_abs = res_max = res_signed = 0.0
    prob_abs = prob_max = 0.0
    sxx = syy = sxy = sx = sy = 0.0
    pocket_agree = 0
    identical_pred = 0
    rel_err = []
    for pa in sorted(glob.glob(os.path.join(a_dir, "*_predictions.csv"))):
        label = os.path.basename(pa)[: -len("_predictions.csv")]
        pb = os.path.join(b_dir, label + "_predictions.csv")
        if not os.path.exists(pb):
            continue
        n += 1
        A, B = read_csv(pa), read_csv(pb)
        same_count += len(A) == len(B)
        cen = lambda p: (float(p["center_x"]), float(p["center_y"]), float(p["center_z"]))
        if A and B:
            top1.append(math.dist(cen(A[0]), cen(B[0])))
        identical_pred += open(pa).read() == open(pb).read()
        for i, p in enumerate(A):
            total_pockets += 1
            if i < len(B):
                d = math.dist(cen(p), cen(B[i]))
                match1 += d <= 1.0
                match4 += d <= 4.0
            # score error for the same pocket (matched by center, any rank)
            cands = [q for q in B if math.dist(cen(p), cen(q)) <= 1.0]
            if cands and float(p["score"]) > 0:
                rel_err.append((float(cands[0]["score"]) - float(p["score"])) / float(p["score"]))
        RA = read_csv(os.path.join(a_dir, label + "_residues.csv"))
        RB = read_csv(os.path.join(b_dir, label + "_residues.csv"))
        for ra, rb in zip(RA, RB):
            x, y = float(ra["score"]), float(rb["score"])
            res_n += 1
            res_abs += abs(x - y)
            res_signed += y - x
            res_max = max(res_max, abs(x - y))
            sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y
            dp = abs(float(ra["probability"]) - float(rb["probability"]))
            prob_abs += dp
            prob_max = max(prob_max, dp)
            pocket_agree += ra["pocket"] == rb["pocket"]
    cov = sxy / res_n - (sx / res_n) * (sy / res_n)
    r = cov / math.sqrt((sxx / res_n - (sx / res_n) ** 2) * (syy / res_n - (sy / res_n) ** 2))
    top1.sort()
    print(f"structures: {n}, identical predictions.csv: {identical_pred / n:.4f}")
    if rel_err:
        print(f"pocket score (same pocket): mean relative error {sum(abs(e) for e in rel_err) / len(rel_err):.5f}, mean bias {sum(rel_err) / len(rel_err):+.5f}")
    print(f"pockets: same pocket count {same_count / n:.4f}, rank-wise center match <=1A {match1 / total_pockets:.4f}, <=4A {match4 / total_pockets:.4f}")
    print(f"top-1 pocket center shift: median {top1[len(top1) // 2]:.3f} A, frac <=1A {sum(d <= 1 for d in top1) / len(top1):.4f}, frac <=4A {sum(d <= 4 for d in top1) / len(top1):.4f}")
    print(f"residue score: pearson {r:.6f}, mean |diff| {res_abs / res_n:.5f}, mean diff (bias) {res_signed / res_n:+.5f}, max |diff| {res_max:.4f}")
    print(f"residue probability: mean |diff| {prob_abs / res_n:.5f}, max |diff| {prob_max:.4f}; pocket label agreement {pocket_agree / res_n:.5f}")


if __name__ == "__main__":
    main()
