# The bistable run (decision 57 addendum 3, task 2; d57_bistable.py): P2_N0900_eta0.780000_G2_pin0, alpha 2 (M 1800), run 9 (seed 3661301716)

| quantity | run 9 | its nine siblings: median (range) |
|---|---|---|
| c (this run's own) | 26.692 | 24.160 (23.148 to 24.733) |
| nu_d | 0.19652 | 0.17786 (0.17043 to 0.18209) |
| tau_r (own fit) | 73.6 | 49.9 (41.1 to 98.8) |
| Q (own fit) | 45.4 | 28.1 (22.4 to 52.9) |
| spectrum share [%] | 86.4 | 82.3 (77.1 to 91.7) |
| sd / sd_E1 | 0.86 | 1.01 (0.86 to 1.11) |
| divider's mean displacement [sigma] | 0.0436 | 0.0015 (-0.0100 to 0.0150) |
| drift per 1000 [sigma] | -0.0013 | 0.0050 (-0.0161 to 0.0163) |
| Z | 20.52 | 20.54 (20.20 to 20.71) |
| psi6, record: left | 0.746 | 0.716 (0.653 to 0.765) |
| psi6, record: right | 0.750 | 0.689 (0.629 to 0.772) |
| orientation - pin [deg], record: left | 0.40 | -0.02 (-0.27 to 0.28) |
| orientation - pin [deg], record: right | -0.34 | -0.02 (-0.27 to 0.20) |
| e at the start [%]: left | -0.088 | -0.098 (-0.115 to -0.083) |
| e at the start [%]: right | -0.098 | -0.081 (-0.091 to -0.058) |
| e at the hold's end [%]: left | 2.301 | -2.374 (-2.928 to -1.857) |
| e at the hold's end [%]: right | -2.643 | -2.328 (-2.900 to -1.995) |
| e, the record's mean [%]: left | 2.056 | -2.335 (-2.528 to -2.052) |
| e, the record's mean [%]: right | -2.136 | -2.393 (-2.557 to -2.127) |
| e's sd over the record [%]: left | 0.436 | 0.326 (0.248 to 0.430) |
| e's sd over the record [%]: right | 0.359 | 0.343 (0.253 to 0.412) |
| e at the end [%]: left | 1.936 | -2.561 (-2.864 to -2.049) |
| e at the end [%]: right | -1.544 | -2.513 (-3.079 to -1.967) |
| ordered fraction at the end: left | 0.930 | 0.964 (0.934 to 0.997) |
| ordered fraction at the end: right | 0.990 | 0.957 (0.910 to 0.984) |
| line direction at the end (0 = x): left / right | 0 / 0 | 9 of 9 x / 9 of 9 x |
| divider - L0 at the end [sigma] | 0.0483 | 0.0095 (-0.0173 to 0.0312) |

- Run 9's left compartment: the first snapshot with e > +1 % at t = 128 sigma-time (snapshots every 64); e > 0 in 100 of its 101 snapshots (the hold and the record); its nine siblings' left compartments: e > 0 in 0 of their 900 snapshots after the first.

**The (cell, mass) E1 with and without run 9** (s1_verdict_261012.per_mass's steps; the cell's E1 over the 5 masses by s1_verdict_261012.e1, the other four masses unchanged):

| the runs | used | nu_d-bar | tau_r | Q | c (alpha 2) ± SE | the cell's E1 over the 5 masses ± SE (chi2 / dof) |
|---|---|---|---|---|---|---|
| all 10 | 10 | 0.17889 | 26.5 | 14.9 | 24.309 ± 0.310 | 23.889 ± 0.056 (7.21 / 4) |
| without run 9 | 9 | 0.17693 | 28.1 | 15.6 | 24.042 ± 0.176 | 23.890 ± 0.050 (6.13 / 4) |

**The scan of all 800 compartments** (the record-mean strain against the median of its cell and compartment over the cell's 50 runs): 1 lie more than 1 % away; the largest |deviation| in thermal sd: 13.3 (the next 1.3).

| cell | M | run | compartment | e, the record's mean [%] | the median of its cell and compartment [%] | deviation [thermal sd] |
|---|---|---|---|---|---|---|
| eta0.780000_G2_pin0 | 1800 | 9 | left | +2.056 | -2.306 | +13.3 |
