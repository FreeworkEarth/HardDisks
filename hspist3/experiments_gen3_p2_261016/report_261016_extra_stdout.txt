# The tables added to Report D57-2 (p2_report_extra_mac.py; decision 57 addendum 2, task 1)

- runs: 400; task logs with a CPU line: 400 of 400; runs whose node in run.log equals their task log's: 400; runs whose seed equals their task line's: 400

## A. The time per CPU model

Minutes per run: the worker's clock from the binary's start to the run's record (the binary and its file moves; run.log). Microseconds per executed event: those seconds over the run's executed events. The D56-1 timing model: 3.533 on E5-2670 v3, 4.242 on E5-4640 v2.

| CPU model | runs | minutes per run: median (range) | minutes at alpha 20: median (largest) | microseconds per event: median (range) | median against E5-2670 v3's |
|---|---|---|---|---|---|
| EPYC 9654 | 59 | 4.35 (2.90 to 6.38) | 5.69 (6.38) | 2.139 (1.916 to 2.600) | 0.602 |
| Gold 6342 | 41 | 5.10 (3.03 to 6.92) | 5.18 (6.92) | 2.209 (2.054 to 2.535) | 0.621 |
| EPYC 7543P | 32 | 5.04 (3.23 to 6.15) | 5.33 (5.40) | 2.284 (2.119 to 2.718) | 0.642 |
| Gold 6248 | 133 | 5.28 (3.23 to 8.00) | 7.05 (8.00) | 2.764 (2.157 to 3.113) | 0.777 |
| E5-2670 v3 | 125 | 6.25 (4.43 to 9.92) | 8.63 (9.92) | 3.555 (2.979 to 3.833) | 1.000 |
| E5-4640 v2 | 10 | 8.12 (5.17 to 8.93) | - (-) | 3.935 (3.493 to 4.607) | 1.107 |

- core-hours used (the runs' seconds): 36.85; the package's prediction 52.9 on E5-2670 v3, 63.5 on E5-4640 v2; the longest run 9.92 min (E5-2670 v3, eta0.780000_G1_pin0, alpha 20) against the array's limit of 25 min

## B. The executed events against the package's prediction, per cell

The prediction: gen_p2_koa_261016.plan() (r = 1.2004 (Z_upper - 1) + 0.0112 per disk and sigma-time with Z_upper = Liu's solid branch x 1.268, times N (T_eq + T_rec) and the rehearsal's factor per pin: 1.171 for pin 0, 1.002 for turned). Z measured: the hold's last quarter (mean over the runs).

| eta | G | pin | Z_upper (model) | Z measured | executed / predicted: median (range) | median per mass (alpha 1, 2, 5, 10, 20) |
|---|---|---|---|---|---|---|
| 0.72 | G1 | 0 | 12.72 | 10.92 | 0.734 (0.725 to 0.738) | 0.734, 0.735, 0.732, 0.731, 0.734 |
| 0.72 | G1 | turned | 12.72 | 10.84 | 0.848 (0.846 to 0.853) | 0.848, 0.847, 0.848, 0.848, 0.848 |
| 0.72 | G2 | 0 | 12.72 | 10.90 | 0.731 (0.721 to 0.737) | 0.732, 0.733, 0.733, 0.732, 0.728 |
| 0.72 | G2 | turned | 12.72 | 10.79 | 0.844 (0.842 to 0.848) | 0.844, 0.844, 0.843, 0.843, 0.844 |
| 0.78 | G1 | 0 | 18.13 | 20.34 | 0.981 (0.956 to 0.999) | 0.982, 0.983, 0.983, 0.976, 0.984 |
| 0.78 | G1 | turned | 18.13 | 18.07 | 1.006 (0.986 to 1.018) | 1.011, 1.004, 1.006, 1.004, 1.003 |
| 0.78 | G2 | 0 | 18.13 | 20.51 | 0.988 (0.970 to 0.998) | 0.986, 0.991, 0.990, 0.987, 0.988 |
| 0.78 | G2 | turned | 18.13 | 16.69 | 0.927 (0.924 to 0.933) | 0.928, 0.927, 0.925, 0.927, 0.927 |

## C. The orientation per cell

The close-packed line direction per compartment (x or y; p2_lattice.config_strain) in the initial state, the hold's last snapshot and the last snapshot, counted over the cell's 100 compartments (50 runs x 2); the record's psi6 orientation (psi6_t, the circular mean mod 60) within 7.5 deg of 0 (lines along x) or of 30 (along y); for the compartments whose direction changes, the first snapshot with the other direction (the direction from theta0: y when theta0 mod 60 lies within 15 deg of 30), its time from the start [sigma-time] and its phase (the hold lasts 5000).

| eta | G | pin | start: x / y | the hold's end: x / y | the end: x / y | record orientation: near 0 / near 30 | compartments that change | the first snapshot with the other direction: t, median (range) | in the hold / in the record |
|---|---|---|---|---|---|---|---|---|---|
| 0.72 | G1 | 0 | 100 / 0 | 0 / 100 | 0 / 100 | 0 / 100 | 100 of 100 | 117 (67 to 337) | 100 / 0 |
| 0.72 | G1 | turned | 0 / 100 | 0 / 100 | 0 / 100 | 0 / 100 | 0 of 100 | - | 0 / 0 |
| 0.72 | G2 | 0 | 100 / 0 | 0 / 100 | 0 / 100 | 0 / 100 | 100 of 100 | 135 (67 to 437) | 100 / 0 |
| 0.72 | G2 | turned | 0 / 100 | 0 / 100 | 0 / 100 | 0 / 100 | 0 of 100 | - | 0 / 0 |
| 0.78 | G1 | 0 | 100 / 0 | 100 / 0 | 100 / 0 | 100 / 0 | 0 of 100 | - | 0 / 0 |
| 0.78 | G1 | turned | 0 / 100 | 0 / 100 | 0 / 100 | 0 / 100 | 0 of 100 | - | 0 / 0 |
| 0.78 | G2 | 0 | 100 / 0 | 100 / 0 | 100 / 0 | 100 / 0 | 0 of 100 | - | 0 / 0 |
| 0.78 | G2 | turned | 0 / 100 | 0 / 100 | 0 / 100 | 0 / 100 | 0 of 100 | - | 0 / 0 |

## D. The entries of test (iv) with R > 2

Per (cell, mass, compartment) with R > 2: the run farthest from the median end strain of its 10, its end strain, its mean over the record's snapshots, its strain at the hold's end, and the end strains of the other 9 (strains in %).

| cell | alpha | compartment | R | run (seed) | end e | record mean e | e at the hold's end | the other 9: end e |
|---|---|---|---|---|---|---|---|---|
| eta0.780000_G2_pin0 | 2 | left | 4.32 | 9 (3661301716) | 1.936 | 2.056 | 2.301 | -2.864 to -2.049 |
