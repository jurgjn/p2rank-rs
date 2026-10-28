#!/bin/bash
# Runs Java P2Rank and prank-rs `predict` on the same dataset (or single structure) and compares
# all output files byte by byte (run.log excluded).
#
# usage: tools/compare.sh <dataset.ds | structure file> [extra prank-rs args...]
# env:   OUTDIR    output base dir (default: /tmp/prank-compare)
#        SKIP_JAVA=1  reuse existing Java output in $OUTDIR/java
#
# Requires a built distro (./gradlew assemble) and `cargo build --release`.
set -u
INPUT=$1; shift
HERE=$(cd "$(dirname "$0")/.." && pwd)
P2=$(cd "$HERE/.." && pwd)
OUT=${OUTDIR:-/tmp/prank-compare}
J=$OUT/java
R=$OUT/rust
if [ "${SKIP_JAVA:-0}" != "1" ]; then rm -rf "$J"; fi
rm -rf "$R"; mkdir -p "$J" "$R"

if [ -f "$INPUT" ] && [[ "$INPUT" != *.ds ]]; then ARGS=(-f "$INPUT"); else ARGS=("$INPUT"); fi

if [ "${SKIP_JAVA:-0}" != "1" ]; then
  /usr/bin/time -f "java: %e s" "$P2/distro/prank" predict "${ARGS[@]}" -o "$J" > "$OUT/java.stdout" 2>&1
  tail -1 "$OUT/java.stdout"
fi
/usr/bin/time -f "rust: %e s" env P2RANK_HOME="$P2/distro" "$HERE/target/release/prank-rs" predict "${ARGS[@]}" -o "$R" "$@" > "$OUT/rust.stdout" 2> "$OUT/rust.stderr"
tail -1 "$OUT/rust.stderr"

ok=0; bad=0
while IFS= read -r f; do
  if [ ! -f "$R/$f" ]; then echo "MISSING $f"; bad=$((bad+1)); continue; fi
  if cmp -s "$J/$f" "$R/$f"; then ok=$((ok+1)); else echo "DIFF $f"; bad=$((bad+1)); fi
done < <(cd "$J" && find . -type f ! -name run.log | sort)
while IFS= read -r f; do
  [ -f "$J/$f" ] || { echo "EXTRA $f"; bad=$((bad+1)); }
done < <(cd "$R" && find . -type f ! -name run.log | sort)
echo "identical=$ok different=$bad"
