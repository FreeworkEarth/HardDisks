# M2 harness, tests and scripts (261012 § 4.7.6), 2026-10-09 HST, this Mac (Apple M3 Max, macOS 26.6.2, Apple clang 17.0.0, arm64)

## Sources (engine-gen3, the commit that adds this folder)

- `hspist3/edmd_core/edmd_gen3.[ch]`: the generation-3 engine with dividers and pistons as bands (M2).
- `hspist3/edmd_core/tests/gen3_m2_harness.c`: the M2 harness (cells listed in its header).
- `hspist3/edmd_core/tests/gen3_body_rule_test.c`: white-box test of the spring divider's contact rule against an independent scan.
- `hspist3/edmd_core/tests/gen3_m1_harness.c`: the M1 harness, unchanged.
- `hspist3/validation/gen3_tolerances_261009.py` (amendment b) and `hspist3/validation/gen3_m2_summary_261009.py` (one row per cell): read the outputs below.

## Builds, from `hspist3/` (every build: 0 warnings under `-Wall -Wextra`)

```
cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o gen3_m2 edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
cc -O3 -ffp-contract=off -Wall -Wextra -o gen3_m1 edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o gen3_body_rule_test edmd_core/tests/gen3_body_rule_test.c -lm
cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -c edmd_core/edmd_gen3.c          (the engine alone)
```

## Runs, one process each, in this order, in a scratch folder (edmd.c reads `HD_CONTACT_AUDIT` once per process)

```
./gen3_m1 audit                     > m1_audit_on_m2.txt          (M1 cells on the M2 engine; compared with ../experiments_gen3_m1_261008/m1_audit_output.txt)
./gen3_body_rule_test 20000         > body_rule_test_output.txt
./gen3_m2 audit --quick             > m2_audit_quick_output.txt   (the same cells, T / 10 and the first 2000 events fully audited)
./gen3_m2 audit                     > m2_audit_output.txt
./gen3_m2 speed                     > m2_speed_output.txt         (no other run of mine; other applications were running: m2_speed_load.txt)
python3 validation/gen3_tolerances_261009.py  > gen3_tolerances_output.txt
python3 validation/gen3_m2_summary_261009.py  > m2_summary_output.txt
```

- `m1_audit_on_m2.txt` is not kept here: it is byte-identical to `experiments_gen3_m1_261008/m1_audit_output.txt` (SHA-256 in `m1_identity.txt`).
- The speed run writes gen2's event log `gen3_m2_speed_evlog.csv` into its working folder (the scratch folder; not kept).

## The runs whose outputs are here (2026-10-09 HST; engine and harness sources as committed in engine-gen3 ddae96c)

- 00:23–00:26: `gen3_m1 audit` (144.7 s) and `gen3_body_rule_test 20000` (about 10 s).
- 00:36–00:41: `gen3_m2 audit --quick` (30.3 s) and `gen3_m2 audit` (227.4 s), with the final harness.
- 00:35–00:36: `gen3_m2 speed` (32.0 s). It used the harness one step before the final one; that step changed only the information-only period estimator, which the speed mode does not call. Load averages before and after: `m2_speed_load.txt`.
- Then `python3 validation/gen3_m2_summary_261009.py > m2_summary_output.txt` and `python3 validation/gen3_tolerances_261009.py > gen3_tolerances_output.txt`.

## Changes after the harness's first full run (00:13; scratch only, every quoted output is from the final sources)

- **Engine (once):** the O(1) jump of the spring divider's contact search for slow approaches (261012 § 4.7.6, log).
- **Harness:**
  - added the late cradle cell (`cradle_round_late`) and the horizon of each class's max |dt|;
  - then two edits of the information-only period estimator (a detrended, windowed spectrum; "does not move" for a held divider).
  - The line diffs of the full outputs before and after those two edits showed differences only in the observables block.
