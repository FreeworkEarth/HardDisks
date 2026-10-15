# d56_b_lattice.py: the lattice start of the 14 Stage 2 cells (decision 56 task 2 b)

binary: experiments_gen3_d50_261014/bin_0544b7e8-dirty/00ALLINONE; per cell three starts (default, pin 0, pin turned), the start only (pa_orient.cmd), plus --gen3-dump-initial. Geometry: H = 10 sqrt(N/100), L0 = N pi / (8 H eta) per side, the divider at L0.

## Per start: the checks

| cell | start | exit, health | M4C (left / right) | left = right | smallest gaps in the dumped start: disks, walls, divider [sigma] | N_L / N_R | divider - L0 [sigma] |
|---|---|---|---|---|---|---|---|
| S2_N0900_eta0.690000 | default | 0, clean | - / - | yes | 0.0637, 0.0312, 0.0308 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.690000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0605, 0.0308, 0.0307 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.690000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0608, 0.0288, 0.0285 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.700000 | default | 0, clean | - / - | yes | 0.0605, 0.0302, 0.0306 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.700000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0609, 0.0289, 0.0327 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.700000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0601, 0.0284, 0.0294 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.720000 | default | 0, clean | - / - | yes | 0.0536, 0.0268, 0.0306 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.720000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0466, 0.0230, 0.0225 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.720000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0577, 0.0285, 0.0280 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.740000 | default | 0, clean | - / - | yes | 0.0422, 0.0220, 0.0211 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.740000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0317, 0.0152, 0.0156 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.740000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0428, 0.0219, 0.0261 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.760000 | default | 0, clean | - / - | yes | 0.0323, 0.0162, 0.0155 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.760000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0270, 0.0138, 0.0136 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.760000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0325, 0.0148, 0.0166 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.780000 | default | 0, clean | - / - | yes | 0.0222, 0.0109, 0.0109 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.780000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0223, 0.0118, 0.0104 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.780000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0204, 0.0101, 0.0100 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.800000 | default | 0, clean | - / - | yes | 0.0199, 0.0092, 0.0090 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.800000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0128, 0.0066, 0.0061 | 450 / 450 | +0.0e+00 |
| S2_N0900_eta0.800000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0186, 0.0087, 0.0097 | 450 / 450 | +0.0e+00 |
| S2_N1600_eta0.690000 | default | 0, clean | - / - | yes | 0.0695, 0.0332, 0.0367 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.690000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0677, 0.0337, 0.0341 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.690000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0706, 0.0354, 0.0327 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.700000 | default | 0, clean | - / - | yes | 0.0634, 0.0318, 0.0353 | 800 / 800 | +3.6e-15 |
| S2_N1600_eta0.700000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0626, 0.0311, 0.0316 | 800 / 800 | +3.6e-15 |
| S2_N1600_eta0.700000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0625, 0.0313, 0.0361 | 800 / 800 | +3.6e-15 |
| S2_N1600_eta0.720000 | default | 0, clean | - / - | yes | 0.0504, 0.0242, 0.0262 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.720000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0476, 0.0233, 0.0243 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.720000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0504, 0.0249, 0.0276 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.740000 | default | 0, clean | - / - | yes | 0.0433, 0.0216, 0.0221 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.740000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0399, 0.0193, 0.0232 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.740000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0435, 0.0213, 0.0207 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.760000 | default | 0, clean | - / - | yes | 0.0436, 0.0204, 0.0197 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.760000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0317, 0.0151, 0.0156 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.760000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0415, 0.0209, 0.0219 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.780000 | default | 0, clean | - / - | yes | 0.0305, 0.0145, 0.0159 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.780000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0254, 0.0127, 0.0128 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.780000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0303, 0.0150, 0.0146 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.800000 | default | 0, clean | - / - | yes | 0.0206, 0.0095, 0.0103 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.800000 | pin0 | 0, clean | lines along x / lines along x | yes | 0.0169, 0.0084, 0.0085 | 800 / 800 | +0.0e+00 |
| S2_N1600_eta0.800000 | pinturned | 0, clean | lines along y / lines along y | yes | 0.0198, 0.0103, 0.0097 | 800 / 800 | +0.0e+00 |

42 of 42 starts exit 0 with a clean health record and no overlap; the replica of g3_lattice_fit gives the binary's arrangement and spacings in 42 of 42

## Per cell: the arrangement, its strain and where it would be strain free

| cell | aspect W/H | default picks | pin 0: arrangement, vac, gap s, strain e, eta* | pin turned: arrangement, vac, gap s, strain e, eta* |
|---|---|---|---|---|
| S2_N0900_eta0.690000 | 0.5683 | 0 | para alt 29x16, 0, 0.118, -2.3 %, 0.7053 | para alt 17x27, 1, 0.111, +3.5 %, 0.7126 |
| S2_N0900_eta0.700000 | 0.5602 | 0 | para alt 29x16, 0, 0.114, -0.8 %, 0.7053 | para alt 17x27, 1, 0.111, +1.9 %, 0.7126 |
| S2_N0900_eta0.720000 | 0.5446 | turned | para alt 29x16, 0, 0.089, +2.2 %, 0.7053 | para alt 17x27, 1, 0.102, -1.1 %, 0.7126 |
| S2_N0900_eta0.740000 | 0.5298 | turned | para alt 29x16, 0, 0.060, +5.1 %, 0.7053 | para alt 17x27, 1, 0.080, -3.9 %, 0.7126 |
| S2_N0900_eta0.760000 | 0.5159 | turned | para alt 28x17, 12, 0.051, -2.2 %, 0.7763 | para alt 17x27, 1, 0.059, -6.5 %, 0.7126 |
| S2_N0900_eta0.780000 | 0.5026 | 0 | para alt 28x17, 12, 0.040, +0.5 %, 0.7763 | para alt 17x27, 1, 0.039, -9.0 %, 0.7126 |
| S2_N0900_eta0.800000 | 0.4900 | turned | para equ 27x17, 9, 0.024, -3.6 %, 0.8282 | para alt 16x29, 6, 0.034, +1.7 %, 0.8126 |
| S2_N1600_eta0.690000 | 0.5685 | turned | para alt 39x21, 0, 0.126, -1.4 %, 0.6990 | para equ 23x35, 5, 0.127, +0.7 %, 0.6944 |
| S2_N1600_eta0.700000 | 0.5604 | 0 | para alt 39x21, 0, 0.121, +0.1 %, 0.6990 | para equ 23x35, 5, 0.120, -0.8 %, 0.6944 |
| S2_N1600_eta0.720000 | 0.5448 | turned | para alt 39x21, 0, 0.090, +3.1 %, 0.6990 | para equ 23x35, 5, 0.097, -3.7 %, 0.6944 |
| S2_N1600_eta0.740000 | 0.5300 | turned | para alt 38x22, 17, 0.075, -1.6 %, 0.7517 | para alt 22x37, 3, 0.081, +2.3 %, 0.7564 |
| S2_N1600_eta0.760000 | 0.5161 | turned | para alt 38x22, 17, 0.059, +1.1 %, 0.7517 | para alt 22x37, 3, 0.077, -0.5 %, 0.7564 |
| S2_N1600_eta0.780000 | 0.5028 | turned | para equ 37x22, 14, 0.049, -1.2 %, 0.7893 | para alt 22x37, 3, 0.056, -3.1 %, 0.7564 |
| S2_N1600_eta0.800000 | 0.4902 | turned | para alt 36x23, 10, 0.033, -3.7 %, 0.8290 | para alt 22x37, 3, 0.037, -5.6 %, 0.7564 |

Entries: rows' orientation (para/perp), rows (equ/alt), n_s x n_r; crystal 0: close-packed lines along x; turned: along y (psi6 orientation 30 deg). vac: the vacancies (sites - N/2). gap s: the lattice's smallest surface gap before the jitter [sigma]. e: the spacing between neighbouring close-packed lines over (sqrt 3 / 2) x the spacing along them, minus 1 (negative: the lines pressed together, the crystal compressed across them and stretched along them). eta*: the packing fraction at which the same arrangement in the same height H is strain free (the compartment width W* that makes e = 0, by the binary's own spacing rules; eta* = N pi / (4 (2 W* + t) H)). W = L0 - t/2, t = 0.05.
