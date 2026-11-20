#!/bin/bash
VPR=/home/claude/vtr9/build/vpr/vpr
declare -A W=([or1200]=72 [clma]=96 [pdc]=94 [riscv_core_lut6]=90)
SW=( "ipf0.1|--initial_pres_fac 0.1" "ipf2|--initial_pres_fac 2.0" "ipf10|--initial_pres_fac 10"
     "pfm1.1|--pres_fac_mult 1.1" "pfm2|--pres_fac_mult 2.0" "pfm3|--pres_fac_mult 3.0"
     "acc0.2|--acc_fac 0.2" "acc2|--acc_fac 2.0" "acc5|--acc_fac 5.0"
     "fipf0.5|--first_iter_pres_fac 0.5"
     "astar1.0|--astar_fac 1.0" "astar2|--astar_fac 2.0"
     "maxcrit0.9|--max_criticality 0.9" "maxcrit0.999|--max_criticality 0.999"
     "cexp2|--criticality_exp 2" "cexp4|--criticality_exp 4" "cexp8|--criticality_exp 8"
     "bbf1|--bb_factor 1" "bbf10|--bb_factor 10" )
for c in "$@"; do
 for s in 1 2 3 4 5; do
  echo "$VPR runs/B1/$c/s$s $c --route_chan_width ${W[$c]} --seed $s"
  echo "$VPR runs/B3a/$c/s$s $c --route_chan_width ${W[$c]} --seed $s --place_algorithm bounding_box"
  echo "$VPR runs/B3b/$c/s$s $c --route_chan_width ${W[$c]} --seed $s --place_algorithm bounding_box --max_criticality 0"
 done
 for e in "${SW[@]}"; do n=${e%%|*}; a=${e#*|}
  for s in 1 2 3 4 5; do echo "$VPR runs/B2/$c/$n/s$s $c --route_chan_width ${W[$c]} --seed $s $a"; done
 done
done
