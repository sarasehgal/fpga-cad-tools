#!/bin/bash
# usage: run1.sh <vpr> <outdir> <circuit> [vpr args...]
VPR=$1; OUT=$2; C=$3; shift 3
mkdir -p "$OUT"; cd "$OUT" || exit 1
[ -f vpr_stdout.log ] && grep -q "VPR succeeded" vpr_stdout.log && exit 0
S=$(date +%s.%N); "$VPR" /home/claude/e1/bench/k6_N10_mem32K_40nm.xml /home/claude/e1/bench/$C.blif "$@" > vpr_stdout.log 2>&1; echo "$(echo "$(date +%s.%N) - $S" | bc) s" > time.txt
rm -f *.net *.route *.post_routing *.blif *.rpt 2>/dev/null
