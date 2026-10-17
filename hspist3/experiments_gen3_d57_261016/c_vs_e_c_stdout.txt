# c against the relaxed strain e (decision 57 addendum 3, task 1; d57_c_vs_e.py)

e: the mean over the record's snapshots of both compartments' strain, then over the runs (SE = sd / sqrt(n)); strains in %. c: E1 over the 5 masses with its SE (Report D57-2). The pilot: alpha 5 only, 2 runs; its error bar is half the difference of its two runs' c.

| eta | G | pin (start) | H [sigma] | runs | e (record) ± SE | sd over runs | e left / right | e at the hold's end | e at the end | lines along y at the end | psi6 (record) | Z | c ± SE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.72 | G1 | turned | 30.0000 | 50 | -0.943 ± 0.010 | 0.070 | -0.953 / -0.932 | -0.940 | -0.953 | 100 of 100 | 0.647 | 10.84 | 14.189 ± 0.013 |
| 0.72 | G1 | 0 | 30.0000 | 50 | -0.822 ± 0.012 | 0.085 | -0.816 / -0.828 | -0.981 | -0.871 | 100 of 100 | 0.627 | 10.92 | 14.184 ± 0.011 |
| 0.72 | G2 | turned | 29.8333 | 50 | -0.133 ± 0.010 | 0.074 | -0.144 / -0.123 | -0.056 | -0.081 | 100 of 100 | 0.656 | 10.79 | 13.752 ± 0.013 |
| 0.72 | G2 | 0 | 29.7083 | 50 | +0.589 ± 0.016 | 0.112 | +0.587 / +0.590 | +0.659 | +0.580 | 100 of 100 | 0.638 | 10.90 | 13.512 ± 0.014 |
| 0.76 | pilot | turned | 30.0000 | 2 | +4.556 ± 0.007 | 0.010 | +4.469 / +4.642 | - | - | 4 of 4 | 0.707 | 15.73 | 15.958 ± 0.046 (indicative) |
| 0.78 | G2 | 0 | 29.9167 | 50 | -2.267 ± 0.047 | 0.330 | -2.228 / -2.307 | -2.230 | -2.254 | 0 of 100 | 0.702 | 20.51 | 23.889 ± 0.056 |
| 0.78 | G1 | 0 | 30.0000 | 50 | -1.867 ± 0.012 | 0.082 | -1.862 / -1.871 | -1.855 | -1.913 | 0 of 100 | 0.722 | 20.34 | 24.287 ± 0.038 |
| 0.78 | G2 | turned | 28.6667 | 50 | -0.336 ± 0.007 | 0.048 | -0.329 / -0.343 | -0.296 | -0.343 | 100 of 100 | 0.894 | 16.69 | 22.176 ± 0.011 |
| 0.78 | G1 | turned | 30.0000 | 50 | +2.459 ± 0.012 | 0.088 | +2.463 / +2.456 | +2.421 | +2.502 | 100 of 100 | 0.810 | 18.07 | 19.698 ± 0.027 |

The slopes d ln c / d e (the % change of c per % of e; a weighted straight line in ln c, effective variance; the scaled SE multiplies by sqrt(chi2/dof) when that exceeds 1):

| the points | n | e range [%] | d ln c / d e ± SE | SE scaled | chi2 / dof | c at e = 0 (the line) |
|---|---|---|---|---|---|---|
| eta 0.72, the three boxes (H 30 combined, 29.833, 29.708) | 3 | -0.882 to +0.589 | -3.454 ± 0.086 | 0.453 | 27.53 / 1 | 13.748 |
| eta 0.72, the four cells (check) | 4 | -0.943 to +0.589 | -3.454 ± 0.087 | 0.374 | 37.15 / 2 | 13.750 |
| eta 0.72, the segment H 30 combined -> 29.833 | 2 | -0.882 to -0.133 | -4.146 ± 0.163 | 0.163 | 0.00 / 0 | 13.677 |
| eta 0.72, the segment 29.833 -> 29.708 | 2 | -0.133 to +0.589 | -2.442 ± 0.201 | 0.201 | 0.00 / 0 | 13.708 |
| eta 0.72, the three boxes, point 1's SE_e = half the two G1 cells' e difference (sensitivity) | 3 | -0.882 to +0.589 | -3.076 ± 0.143 | 0.605 | 17.87 / 1 | 13.733 |
| eta 0.78, the two turned cells (H 30, 28.667) | 2 | -0.336 to +2.459 | -4.240 ± 0.056 | 0.056 | 0.00 / 0 | 21.863 |

- the two G1 cells at 0.72 (the same box, H 30): e differs by +0.121 ± 0.016 % (+7.7 SE) while c differs by -0.005 ± 0.017; the sensitivity row of the slopes takes half that e difference as point 1's error
- the pilot: its two runs' own c 16.006 and 15.913; record snapshots per run 34, 34; lines along y at its last record snapshot in 4 of 4 compartments
