# M1 harness and E0 pre-check (261012 § 4.7.3), 2026-10-08 HST, this Mac (arm64, Apple clang 17)

## Harness (engine-gen3, the commit that adds this folder; source `hspist3/edmd_core/tests/gen3_m1_harness.c`)

Build, from `hspist3/`:

```
cc -O3 -ffp-contract=off -Wall -Wextra -o gen3_m1 edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
```

- The harness and `edmd_gen3.c` compile without warnings under `-Wall -Wextra`.
- `edmd.c` (the gen2 engine, unchanged) has its pre-existing warnings.

Runs, one process each, in this order (edmd.c reads `HD_CONTACT_AUDIT` once per process):

```
./gen3_m1 speed   > m1_speed_output.txt
./gen3_m1 audit   > m1_audit_output.txt
./gen3_m1 diverge > m1_diverge_output.txt
```

- **Speed ran first,** with no other run of mine on the machine.
- **Other applications were running.** `uptime` load averages were 8.04 8.56 7.98 before and 6.82 8.07 7.84 after. Hence each rate is timed 3 times (median and range).
- **The audit and divergence outputs are deterministic.** Two earlier full runs of the same engine code printed the same lines; one of them was before the z line and the timing repeats were added to the harness.

## E0 pre-check of the gen2 path (the gate's runner, on this Mac)

- **Binaries**, both built with `-O3 -ffp-contract=off` against the same SDL, kissfft and experiment_validation sources:
  - engine-gen3 at 2cdfe04: its 00ALLINONE.c differs from 7b08827 only by the item-0 guard;
  - 7b08827 (scratch build `bin_7b_O3`).
- **Runner**, from `hspist3/` on main:

```
python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin <engine-gen3 build> --ref-bin <7b08827 build> --out <out> --jobs 6 --cases ctrl_min,ctrl_leg
python3 cluster/resched_gate_261005/audit_runs_261007.py report --out <out> > e0_precheck_report.txt
```

The report's last verdict line needs the A-fixed and mode-1 cases, which this pre-check does not run. The per-case "plain vs ref" column is the E0 comparison.
