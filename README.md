# FPGA CAD Tools in C++

Routing, placement and VPR experiments from ECE1387 CAD for Digital Circuit Synthesis and Layout, University of Toronto MEng, Fall 2026.

## maze-router/ (20/20)

Lee-Moore maze router for an island-style FPGA, one C++ file with an EZGL viewer.

- BFS over a routing-resource graph, wire sharing within multi-sink nets, rip-up and re-route that moves the failed net to the front (up to 100 attempts).
- Same code runs three architectures: baseline (Fc = W, Fs = 3), odd/even connection block, and an Fs = 4 switch block.
- An independent checker validates every route: sources reach sinks, only legal switches, no shared segments.

| Circuit | Min W (Fs=3) | Min W without rip-up | Min W (Fs=4) |
|---|---|---|---|
| cct1 | 2 | 5 | 2 |
| cct2 | 4 | 5 | 3 |
| cct3 | 5 | 6 | 4 |
| cct4 | 7 | 8 | 6 |

## analytical-placer/

Clique-model quadratic placement solved with UMFPACK (SuiteSparse 5.10.1), random-I/O sensitivity study over 20 trials, flow-based spreading after Darav et al. (FPGA 2019) with several psi schedules, and anchored re-solves. An EZGL viewer steps through every stage. Build UMFPACK with `build_umfpack.sh`, then `make`.

## vpr-experiments/

VTR 9 / VPR study on or1200, clma, pdc and riscv_core_lut6: anneal cost vs temperature, seed sensitivity (1 to 2% CV), starting-temperature sweep (too cold costs up to 17% HPWL), range-window removal, and router parameter sweeps. 660 seeded runs scripted across 132 configurations; `scripts/` holds the job generators and parsers, `all_results.csv` the aggregated results.

The EZGL graphics library (Vaughn Betz, U of T) and SuiteSparse are not included; see each README for setup.

Author: Sara Sehgal
