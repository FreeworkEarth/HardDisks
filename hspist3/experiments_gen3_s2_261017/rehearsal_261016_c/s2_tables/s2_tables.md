# The S2 tables (s2_tables_mac.py; registration sec. 4.9)

- data: /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d58/rehearsal_3/standin/hspist3/experiments_gen3_stage2_koa_261016/S2
- runs: 4 (with a p2_reduce error: 0); (cell, H): 4; cells: 4; m_<M> folders: 4, with a missing table: 0
- flagged compartments: 0; (cell, H) with an end orientation off the pin: 0; compartments with no strain in any snapshot: 0; runs without e: 0; snapshot compartments without a strain (skipped in the record means): 0
- conventions (not fixed by the registration; the docstring): a run left out for its flag is not 'lost'; c0 as registered needs exactly three heights (else NOT AS REGISTERED); Verdict I is INCONSISTENT when an eta with both c0 as registered gives |z| >= 3.5, else NO VERDICT (a c0 missing) or NO VERDICT AS REGISTERED (a c0 not as registered) before CONSISTENT; e per (cell, H) weighs each run equally (its record mean over the snapshots, both compartments; a snapshot without a strain is skipped; a run with a compartment with no strain in any snapshot has no e); the diagnostics are over all runs of the cell.

## c0 per cell (the straight line of ln c on e through the heights; c0 at e = 0)

| N | eta | pin | heights | c0 ± SE | d ln c / d e ± SE | chi2 / dof | e range [%] | psi6 | Z | sd / sd_E1 | nu_d / nu_pred | low-Q heights | end off the pin | flagged compartments |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 900 | 0.78 | 0 | 0 (NOT AS REGISTERED) | - | - | - / -2 | - to - | 0.937 | 18.13 | - | 1.239 | 0 | 0 | 0 |
| 900 | 0.76 | turned | 0 (NOT AS REGISTERED) | - | - | - / -2 | - to - | 0.853 | 13.95 | - | 1.108 | 0 | 0 | 0 |
| 1600 | 0.72 | turned | 0 (NOT AS REGISTERED) | - | - | - / -2 | - to - | 0.564 | 10.79 | - | 1.050 | 0 | 0 | 0 |
| 1600 | 0.80 | turned | 0 (NOT AS REGISTERED) | - | - | - / -2 | - to - | 0.950 | 19.45 | - | 1.172 | 0 | 0 | 0 |

## c0(1600) - c0(900) per eta (turned)

| eta | c0 N 900 | c0 N 1600 | difference ± SE | z |
|---|---|---|---|---|

## VERDICT I (isotropy), N 900: z = (c0 pin 0 - c0 turned) / sqrt(SE_0^2 + SE_turned^2); both |z| < 3.5 -> ISOTROPY-CONSISTENT

| eta | c0 pin 0 [heights] | c0 turned [heights] | difference ± SE | z |
|---|---|---|---|---|
| 0.76 | - | - | - | - |
| 0.78 | - | - | - | - |

**VERDICT I: NO VERDICT (a c0 or its SE is missing)**

## The literature elastic constants of the hard-disk crystal (descriptive, no verdict)

Not computed by this script: the comparison is written in the S2 report from the literature values.

## Per (cell, H): E1 over the masses and the relaxed strain

e per (cell, H): the mean over its runs of each run's record mean (both compartments; equal weight per run; a snapshot without a strain skipped).

| cell id | runs | runs with an error | masses in E1 | E1: c ± SE | chi2 | low Q | e [%] ± SE | runs without e | compartments with no strain in any snapshot | snapshot compartments without a strain | end off the pin (compartments) | flagged compartments |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S2_N0900_eta0.760000_k0697_pinturned | 1 | 0 | 0 (0 usable, fewer than 5: no E1) | - | - | 0 | -0.389 (no SE) | 0 | 0 | 0 | 0 of 2 | 0 |
| S2_N0900_eta0.780000_k0685_pin0 | 1 | 0 | 0 (0 usable, fewer than 5: no E1) | - | - | 0 | -0.216 (no SE) | 0 | 0 | 0 | 0 of 2 | 0 |
| S2_N1600_eta0.720000_k0866_pinturned | 1 | 0 | 0 (0 usable, fewer than 5: no E1) | - | - | 0 | -0.015 (no SE) | 0 | 0 | 0 | 0 of 2 | 0 |
| S2_N1600_eta0.800000_k1024_pinturned | 1 | 0 | 0 (0 usable, fewer than 5: no E1) | - | - | 0 | -0.165 (no SE) | 0 | 0 | 0 | 0 of 2 | 0 |

## Per (cell, H, mass): E1

| cell id | alpha | used / present | lost | out (> 2 lost) | Q | low Q | c ± SE |
|---|---|---|---|---|---|---|---|
| S2_N0900_eta0.760000_k0697_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N0900_eta0.780000_k0685_pin0 | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N1600_eta0.720000_k0866_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
| S2_N1600_eta0.800000_k1024_pinturned | 5 | 1 / 1 | 9 | OUT | - | 0 | - |
