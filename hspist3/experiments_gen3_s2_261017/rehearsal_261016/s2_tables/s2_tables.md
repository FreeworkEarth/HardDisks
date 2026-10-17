# The S2 tables (s2_tables_mac.py; registration sec. 4.9)

- data: /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d58/rehearsal_1/standin/hspist3/experiments_gen3_stage2_koa_261016/S2
- runs: 4; (cell, H): 4; cells: 4; flagged compartments: 0; (cell, H) with an end orientation off the pin: 0

## c0 per cell (the straight line of ln c on e through the heights; c0 at e = 0)

| N | eta | pin | points | c0 ± SE | d ln c / d e ± SE | chi2 / dof | e range [%] | psi6 | Z | sd / sd_E1 | nu_d / nu_pred | end off the pin | flagged compartments |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 900 | 0.78 | 0 | 0 | - | - | - / -2 | - to - | 0.937 | 18.13 | - | 1.239 | 0 | 0 |
| 900 | 0.76 | turned | 0 | - | - | - / -2 | - to - | 0.853 | 13.95 | - | 1.108 | 0 | 0 |
| 1600 | 0.72 | turned | 0 | - | - | - / -2 | - to - | 0.564 | 10.79 | - | 1.050 | 0 | 0 |
| 1600 | 0.80 | turned | 0 | - | - | - / -2 | - to - | 0.950 | 19.45 | - | 1.172 | 0 | 0 |

## c0(1600) - c0(900) per eta (turned)

| eta | c0 N 900 | c0 N 1600 | difference ± SE | z |
|---|---|---|---|---|

## VERDICT I (isotropy), N 900: z = (c0 pin 0 - c0 turned) / sqrt(SE_0^2 + SE_turned^2); both |z| < 3.5 -> ISOTROPY-CONSISTENT

| eta | c0 pin 0 | c0 turned | difference ± SE | z |
|---|---|---|---|---|
| 0.76 | - | - | - | - |
| 0.78 | - | - | - | - |

**VERDICT I: NO VERDICT (a c0 is missing)**

## Per (cell, H): E1 over the masses and the relaxed strain

| cell id | runs | masses in E1 | E1: c ± SE | chi2 | low Q | e [%] ± SE | end off the pin (compartments) | flagged compartments |
|---|---|---|---|---|---|---|---|---|
| S2_N0900_eta0.760000_k0697_pinturned | 1 | 0 | - | - | 0 | -0.389 (no SE) | 0 of 2 | 0 |
| S2_N0900_eta0.780000_k0685_pin0 | 1 | 0 | - | - | 0 | -0.216 (no SE) | 0 of 2 | 0 |
| S2_N1600_eta0.720000_k0866_pinturned | 1 | 0 | - | - | 0 | -0.015 (no SE) | 0 of 2 | 0 |
| S2_N1600_eta0.800000_k1024_pinturned | 1 | 0 | - | - | 0 | -0.165 (no SE) | 0 of 2 | 0 |

## Per (cell, H, mass): E1

| cell id | alpha | used / present | lost | out (> 2 lost) | Q | low Q | c ± SE |
|---|---|---|---|---|---|---|---|
| S2_N0900_eta0.760000_k0697_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N0900_eta0.780000_k0685_pin0 | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N1600_eta0.720000_k0866_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N1600_eta0.800000_k1024_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
