# 260909 — TASK A: the 504 `wall_overdue` flags, and temperature provenance (Claude Code)

Written 2026-09-09. Every number below was produced this session from the files named; commands are in the session log. No core physics was changed; no accepted trajectory was modified; nothing was committed.

## A.1 Verdict: **(a) — bookkeeping at initialization. Flag is harmless; the trajectories are kept.**

Reached with time stamps, not by inference. Summary of the evidence, then the detail.

| what | value |
|---|---|
| affected campaign / leaves | `campaign_r25_psi6_20260823/eta_0p019635, eta_0p026180, eta_0p039270` (L0 = 200, 150, 100) |
| flagged trajectories (ledger) | 162 + 180 + 162 = **504** of 675; total overdue events **662** |
| per-event time stamps in the leaf logs | **none** — the only lines mentioning "overdue" are the 162/180/162 `EDMD-HEALTH` summaries |
| trajectories reproduced with `--speed-sound-exact-seed` and a driver-side trace | **675 / 675** |
| ledger count reproduced exactly | **675 / 675** (0 mismatches) |
| time of every counter increment | **t = 0**, at initialization, before hold step 1 — `phase=` lines during `wall_hold` or `post_release`: **0** in all 675 runs |
| wall | **RIGHT** wall in all 1,350 face-particle records (never L/B/T, never the divider) |
| overdue by | **exactly 0 px** (`gap = 0`) — the particle centre sits precisely on the face `x = boxW − R` |
| particles involved | two per trajectory: 74 and 99 at L0 = 200 and 100; 66 and 83 at L0 = 150 |
| when the counter increments | iff the face particle's initial velocity points **into** the wall (`v_toward > 0`) |
| resulting count distribution | expected ¼ / ½ / ¼ for 0/1/2 of two particles → 56/112/56 of 225; observed 63/115/47 (L0 = 200), 45/116/64 (150), 63/115/47 (100). The L0 = 200 and L0 = 100 campaigns give **identical** counts because their per-run seeds are identical (same `--seed` base and the same (M, run) index), so the velocity draws are not independent across those densities and the ¼/½/¼ comparison is two samples, not three. Harmless for c_s(η). |
| "moving INTO wall" records | 109 + 216 + 100 + 237 = **662** = the ledger's total overdue events |

### The code path (`hspist3/edmd_core/edmd.c`)

```c
static int wall_time_from_gap(double gap, double speed, double* tcol){
    if (gap <= 0.0) { *tcol = 0.0; return 2; }   /* overdue */
    const double t = gap / speed;
    if (t <= 1e-12) return 0;
    *tcol = t; return 1;
}
...
static void schedule_walls(EDMD* S, int i){
    /* ##CHRIS: rc==2 means the wall collision was already overdue when scheduled. The old
       code discarded exactly those, which is how a particle escaped the box. Counted so the
       exposure is measurable rather than invisible. */
    double t; int rc;
    if((rc = collide_time_wall_L(S, &S->P[i], &t))) {
        if(rc==2) S->wall_overdue_count++;
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WL });
    }
    ... (R, B, T identical)
```
`schedule_walls` is called from `reschedule_all_internal` (line 811, initialization) and per particle after each of its collisions (line 789). The counter is incremented only for the four **outer** walls, only when a particle centre is at or past the face (`gap ≤ 0`) while moving toward it. The event is then scheduled at `t_col = 0` and resolved as an immediate elastic reflection — the physically correct outcome for a disk touching a wall with inward velocity. In one paragraph: **`wall_overdue` counts wall collisions that were already due at the moment they were scheduled.** Its purpose (2026-08-23) was to expose the class of event the pre-fix core silently dropped; a non-zero value means "a due-now wall bounce was handled", not "a collision was missed".

### Why the particles sit on the face, and why only L0 ≥ 100 (`hspist3/00ALLINONE.c`, one-wall initializer)

```c
const float diameter_px = 2.0f * (float)PARTICLE_RADIUS;                 // 5413
const float seed_pad_px = fmaxf(1e-4f, 1e-5f * diameter_px);             // 5414  -> 2.4e-4 px
const float margin_px   = (float)PARTICLE_RADIUS + seed_pad_px;          // 5415
...
float right_x_max = (float)XW2 - margin_px;                              // 5423
const float spacing_right_x = (num_cols_right > 1)
    ? (right_width / (float)(num_cols_right - 1)) : 0.0f;                // 5428
...
X[idx] = (num_cols_right > 1) ? right_x_min + col * spacing_right_x : ...   // placement loop
```
All of this is `float`. The last column is meant to land at `XW2 − R − 2.4e-4 px`. At L0 = 200 the box is ≈ 9800 px wide, where a float's resolution is ≈ 9.8e-4 px, so the 2.4e-4 px pad is rounded away and the column lands at exactly `XW2 − R`: `gap = 0`. At L0 = 75 (≈ 3800 px, resolution 2.4e-4 px) the pad survives, which is why the 2026-08-23 survey found 0 % exposure at L0 ≤ 75 and 66–79 % at L0 ≥ 100. The other three walls keep their pad (trace: min gap L/B/T = 2.44e-4 / 2.40e-4 / 2.44e-4 px in every run). Two particles are affected because the right compartment's lattice has two rows.

### Consequence for the data

- The 2000-step hold precedes release: 2000 × 0.4 px-time = 800 px-time = **33.3 σ-time** at dt = 0.4/24 = 0.01667 σ-time (unit confirmed by the leaf log, `samples=159840 T=2664 sigma-time`). At η = 0.02 that is only ≈ 3 collisions per particle, which is sufficient only because the initial velocities are already Maxwellian and rescaled to kT = 1 (A.2); the single reflection at t = 0 is irrelevant to the measurement window either way. *(Corrected 2026-09-10: an earlier version of this line said 800 σ-time.)*
- No overdue event occurs in any of the 675 trajectories after t = 0. The strict rule ("any health event ⇒ exclude"), applied literally, excluded 504 good trajectories in `analysis_paper1_20260908` because it does not distinguish an initialization bookkeeping event from a mid-run repair. **Recommendation:** keep the strict rule for `forced_advance`, `overlap_repair`, `clamp_repair`, and for `wall_overdue` at t > 0; treat `wall_overdue` at initialization as a seeder artefact to be documented, and retain these 675 trajectories for the sound-speed fits.
- **Open driver item (not implemented; not for existing data; not while any campaign runs):** the float seed pad `seed_pad_px = fmaxf(1e-4f, 1e-5f·d)` in `00ALLINONE.c` (line 5414) must be carried in double, or replaced by a fixed inset of 1e-3·d as `edmd_init_lattice_gas()` does, before any new c_s run at L0 ≥ 100. Nothing changes for the trajectories on disk.

### How the stamps were obtained

A driver-side (not core) diagnostic, gated by `HD_OVERDUE_TRACE=1`, in `00ALLINONE.c` (`##CHRIS`): at the start of the hold loop it scans all particles for `gap ≤ 0` at any outer wall and prints wall, gap and inward speed; after every hold and post-release sample it prints any change of `edmd_backend_wall_overdue_count()` with phase, step and σ-time. Nothing in it touches the dynamics. Reproductions ran with `--speed-sound-exact-seed=<seed>` and otherwise the leaf's `00_COMMAND.md` flags; the `wall_overdue` values reproduced the ledger for all 675 seeds, which also demonstrates that these trajectories are seed-deterministic.

## A.2 Temperature provenance

- `--kbt1` sets `cli_force_kbt_one = 1` (`00ALLINONE.c` line 3341). `kB_effective()` (line 1679) then returns `1/temperature_runtime`, so kB·T = 1 nominally.
- After the Maxwell–Boltzmann draw at `temperature_runtime` (lines 5578/5598), the initializer **rescales all velocities** so that KE = N·K_B·T exactly (lines 5612–5618: `scale = sqrt(target_ke/actual_ke); Vx[i] *= scale; Vy[i] *= scale`), then `equalize_temperature_per_segment(temperature_runtime)` (line 5641) rescales each compartment to `ke_target = n·kB_effective()·T_new` (line 1850). With `--kbt1` this makes **kBT = 1 exactly per compartment at t = 0, by construction**, not a finite-N draw (the ±5 % scatter seen in the pressure runner does not apply: that runner draws and never rescales).
- During the hold the divider is fixed and the outer walls are elastic, so the gas kinetic energy is conserved; T at release equals 1 up to float rounding.
- **Durable record:** none. Neither the wall trace (`Time,Wall_X,Displacement(σ),Left_Count,Right_Count,L0,eta,Center_X(σ),Seed,Target_Oscillations,Predicted_Frequency,Planned_Steps,Planned_Duration`) nor `speed_of_sound_psi6.csv` carries a temperature column. It is, however, exactly reproducible from the seed (see A.1) and exact by construction; if a per-trajectory record is wanted for the paper, the cheapest honest route is to add T at release to the trace as a column in a future run, not to manufacture one now.
- One wrinkle worth a sentence in methods: line 5613 uses the compile-time `K_B` while the per-segment step uses `kB_effective()`; with the default `T = 100, K_B = 0.01` both give kB·T = 1, so there is no inconsistency in the runs on disk, but the two would diverge if `temperature_runtime` were changed under `--kbt1`.
