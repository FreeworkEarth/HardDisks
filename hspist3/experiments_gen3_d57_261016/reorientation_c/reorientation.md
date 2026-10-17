# The reorientation at eta 0.72 (decision 57 addendum 3, task 4; d57_reorientation.py)

## P2 (KOA): the two pin-0 cells, 50 runs x 2 compartments each

| cell | compartments | start x / y | end x / y | snapshots every [sigma] | psi6_t every [sigma] | first snapshot with y: median (range) | first y snapshot before 200 | first y snapshot before 300 | first y snapshot before 500 | first psi6_t within 7.5 deg of 30: median (range) | before 200 | before 300 | before 500 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G1, H 30 | 100 | 100 / 0 | 0 / 100 | 67-117 | 1.0-3.3 | 117 (67 to 337) | 90 of 100 | 99 of 100 | 100 of 100 | 78 (24 to 291) | 97 of 100 | 100 of 100 | 100 of 100 |
| G2, H 29.708 | 100 | 100 / 0 | 0 / 100 | 67-117 | 1.0-3.4 | 135 (67 to 437) | 75 of 100 | 94 of 100 | 100 of 100 | 100 (24 to 413) | 87 of 100 | 97 of 100 | 100 of 100 |

- The first snapshot with y is the first snapshot after the start (t = one snapshot interval) in G1, H 30: 55 of 100, G2, H 29.708: 44 of 100 compartments: the snapshots every 67-117 sigma-time do not show the rotation itself, only that it is over; psi6_t (every 1.0-3.4 sigma-time) shows the orientation's path but not the configurations.

## The Mac run (snapshots every 10 sigma-time; P2's line 21 with T_snap 10)

| t [sigma] | left: theta0 [deg], lines, ordered fraction, e [%] | right: theta0 [deg], lines, ordered fraction, e [%] |
|---|---|---|
| 0 | 0.0, x, 1.00, +2.14 | 60.0, x, 1.00, +2.15 |
| 10 | 0.3, x, 0.84, +1.77 | 0.2, x, 0.90, +1.33 |
| 20 | 0.8, x, 0.82, +0.62 | 0.0, x, 0.83, +0.79 |
| 30 | 0.1, x, 0.73, +1.06 | 3.3, x, 0.79, -0.45 |
| 40 | 59.1, x, 0.67, +0.97 | 28.0, y, 0.76, -0.87 |
| 50 | 59.7, x, 0.71, +0.42 | 29.5, y, 0.72, -0.04 |
| 60 | 59.6, x, 0.75, -0.86 | 29.2, y, 0.82, -0.58 |
| 70 | 31.0, y, 0.77, +0.17 | 29.5, y, 0.86, -0.02 |
| 80 | 30.3, y, 0.75, +1.11 | 30.3, y, 0.86, -1.04 |
| 90 | 30.2, y, 0.76, -0.98 | 30.4, y, 0.83, -0.43 |
| 100 | 30.0, y, 0.83, -0.66 | 30.4, y, 0.86, -1.03 |
| 110 | 30.1, y, 0.88, -1.38 | 30.0, y, 0.80, -2.12 |
| 120 | 29.8, y, 0.85, -1.62 | 30.6, y, 0.81, -0.93 |
| 130 | 30.2, y, 0.89, -0.33 | 29.9, y, 0.86, -0.63 |
| 140 | 29.9, y, 0.89, -0.49 | 29.9, y, 0.84, +0.05 |
| 150 | 29.9, y, 0.80, -0.17 | 30.5, y, 0.90, -0.57 |

- the left compartment: the first 10-sigma snapshot with lines along y at t = 70; the first psi6_t sample within 7.5 deg of 30 at t = 66.7; snapshots with lines along x after that, up to 600: 0

- the right compartment: the first 10-sigma snapshot with lines along y at t = 40; the first psi6_t sample within 7.5 deg of 30 at t = 34.2; snapshots with lines along x after that, up to 600: 0
