#!/usr/bin/env python3
"""Per-phase timing of prank-rs (minimum over N runs, from P2RS_TIMING output).
usage: phase_bench.py <structure> <runs> [prank-rs args]   (env vars are passed through)"""
import os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
P2 = os.path.abspath(os.path.join(HERE, "..", ".."))
f, runs, extra = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
best = {}
order = []
for _ in range(runs):
    env = dict(os.environ, P2RS_TIMING="1", P2RANK_HOME=os.path.join(P2, "distro"))
    out = subprocess.run([os.path.join(P2, "rust/target/release/prank-rs"), "predict", "-f", f, "-o", "/tmp/phase_bench_out"] + extra,
                         env=env, capture_output=True, text=True).stderr
    prev = 0.0
    for m in re.finditer(r"timing:\s+([\d.]+) ms\s+(.+)", out):
        t, name = float(m.group(1)), m.group(2).strip()
        if name not in best:
            order.append(name)
        best[name] = min(best.get(name, 1e9), t - prev)
        prev = t
print(" | ".join(f"{n} {best[n]:.1f}" for n in order), f"| total {sum(best.values()):.1f} ms")
