#!/bin/bash
# usage: mkjobs.sh <circuit...>  -> prints all Part A job lines
VPR=/home/claude/vtr9/build/vpr/vpr
BASE="--pack --place --place_algorithm bounding_box"
for c in "$@"; do
 for s in 1 2 3 4 5; do
  echo "$VPR runs/A2/$c/s$s $c $BASE --seed $s"
  echo "$VPR runs/A4/$c/s$s $c $BASE --seed $s --place_rlim_escape 1.0"
  for T in 1e-2 1e-3 1e-4 1e-5 1e-6; do
   echo "$VPR runs/A3/$c/T$T/s$s $c $BASE --seed $s --init_t $T --exit_t 1e-8 --alpha_t 0.9"
  done
 done
done
