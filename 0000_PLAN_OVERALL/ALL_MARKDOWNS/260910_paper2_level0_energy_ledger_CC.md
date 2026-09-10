# 260910 — Paper 2, Level 0: the energy ledger (Claude Code)

Replaces the 2026-09-10 morning version of this file, which was an inventory of what *could not* be computed. The ledger now closes. Definitions follow `260910_hand_calculations_cs_to_dissipation_COWORK.tex` §7–9. No core physics changed; no accepted trajectory modified; nothing committed.

---

## 1. What was added (logging only)

**Per-sample trace columns** (`00ALLINONE.c`, energy-transfer writer). Appended **after** `SegEtas`, so every existing column index is unchanged (30 → 34):

```
KE_gas_total, KE_gas_left, KE_gas_right, Px_gas
```
`KE_gas_*` from `segment_ke[]` (segment 0 = left, the rest = right — the same convention as `Left_Count`/`Right_Count`); `Px_gas = Σ m v_x`. `X/Vx/Vy` are synced from the EDMD backend and `recompute_segment_stats_counts_and_temperature()` is called earlier in the same block, so the values are current at write time.

**Gated per-event collision log** (`edmd_core/edmd.c`, print-only, inside `resolve_wall`). Enabled by `HD_PISTON_EVENTS=<path>` in energy-transfer mode; the driver passes `PIXELS_PER_SIGMA` so the log carries σ-time and the core stays free of driver constants. Fires at both piston branches (`EV_PL`/`EV_PR`) and both divider branches:

```
t_sigma,kind,u_wall,v_before,v_after,dE,dp        kind ∈ {PL, PR, D0…D4}
```
`dE` is the quantity the core books as work — the **particle's** KE change for a prescribed-velocity wall, the **wall's** KE change for a finite-mass wall. Those are different quantities; the log records which by `kind` and the caller must not mix them. `dp = m(v_after − v_before)` is always the particle's momentum change.

**Regression (item 2d).** Two checks on the final build:

| check | result |
|---|---|
| gate smoke test (the 2-step negative test in `experiments_energy_transfer/00_COMMAND.md`) | exit 2, `requested_particle_count_mismatch` — **correct**: it requests 20 particles with `--particles-boxes=0,19`, and the fail-closed validator catches it |
| accepted speed-of-sound cell re-run (η = 0.261799, L0 = 15, M = 50, run 0, seed 2894599681, `--edmd-acc=0`) | `cmp` vs the accepted trace: **BIT-IDENTICAL**, 10388 rows |

The campaign has no η = 0.30 cell; 0.261799 is the nearest accepted one.

**Units (item 2c), resolved from the writer.**

- `W*_x_sigma` divides by `PIXELS_PER_SIGMA`; `W*_v` does **not** — and does not need to. With `x_σ = x_px/PPS` and `t_σ = t_px/PPS` the scale cancels, so **velocities are numerically identical in both systems**. `PARTICLE_MASS = 1.0` matches the core's hard-coded `m = 1.0`, so energies are in kT under `--kbt1`.
- **`W0` is the spring wall.** `prm.divider_mass[w] = cli_wall_masses[w]`, so `--wall-mass-factors=200;16000` gives W0 = 200 and W1 = 16000; `primary_wall_index = leftmost_wall_index` and the spring acts on `all_wall_positions[leftmost_wall_index]`, i.e. W0. **The `M2_<x>` in the mass-sweep trace filenames is W1's mass, not W0's.**

> **Correction to the morning version of this file.** Its ΔKE_div column paired `W0_v` with the filename's mass. W0's mass is 200 in every row of that sweep; the filename number is W1's. That table's ΔKE_div values are therefore wrong and it also omitted W1 entirely. It is superseded by §3 below, which is measured on runs where the ledger is complete.

---

## 2. Two driver bugs — found, fixed (go given 2026-09-10), gated

**(a) `--wall-positions` used the *default* L0.** In `initialize_simulation()` the requested position became a fraction `frac = pos/(2·L0_UNITS)` while `L0_UNITS` still held its default 20 — `L0_UNITS = cli_override_L0_units` happens ~70 lines later (line 4248). `--l0=39.27 --wall-positions=39.27` therefore landed the wall at 39.27/40 × 78.54 = **77.07 σ** and the fail-closed validator aborted (`initial_wall_position_mismatch`, tolerance 1e-5 σ). Every reference run used `--l0=20`, where the bug is invisible. **Fix:** use `cli_override_L0_units` when it is set. Wall now lands at the centre for any L0.

**(a′) A second, independent limit remains (not a bug, a quantisation).** The box width is an integer pixel count: 2·39.27·24 = 1884.96 → 1884 px = 78.5 σ, so 39.27 σ is not representable and the validator still (correctly) refuses it. **L0 must lie on the pixel grid, i.e. 24·L0 integer.** Nearest to 39.27 is **L0 = 39.25** (η = 0.100051), which is what the pilot uses.

**(b) Every energy-transfer run rewrote `experiments_energy_transfer/00_COMMAND.md`**, because `write_command_md_to_dir()` always targeted the fixed `g_energy_transfer_dir`. My first gate test clobbered the 2026-08-20 provenance file; it was restored verbatim. **Fix:** write it into the directory of `--energy-transfer-trace` when that flag is given; unchanged otherwise.

### Gates (pass/fail)

| gate | result |
|---|---|
| energy-transfer trace, `--l0=20 --wall-positions=20`, seed 7001, before vs after both fixes | **PASS** — byte-identical (`e02bf2c4c1fb35443ad60c711ad98bc7`) |
| its event log, same comparison | **PASS** — byte-identical (`edec0d576a4c2296f4859645f2a9eb6b`) |
| speed-of-sound η = 0.261799, seed 2894599681, vs the accepted trace | **PASS** — BIT-IDENTICAL |
| `00_COMMAND.md` written next to the trace; repo file untouched | **PASS** — repo copy still reads `Timestamp: 2026-08-20 17:14:08 CEST` |
| three regression suites (16 / 27 / 8 assertions) | **PASS** |

Both traces were also confirmed deterministic (two runs of the same build, byte-identical) before the fixes were applied, so the gate compares like with like.

## 2b. The T_i deficit — located, and it is a formula choice

`HD_KE_TRACE=1` prints gas KE and momentum at four points (print-only, `##CHRIS`):

| point | KE_tot | KE_left | KE_right | Px | kT_mean |
|---|---|---|---|---|---|
| 1 after velocity draw | 112.333 | 51.588 | 60.745 | −9.293 | 1.1233 |
| 2a after global kT rescale | **100.000** | 45.924 | 54.076 | −8.768 | **1.0000** |
| 2b after per-segment equalize | **97.678** | 49.430 | 48.248 | **7.4e−15** | **0.97678** |
| 3 end of hold (at release) | 97.678 | 49.430 | 48.248 | 0.0844 | 0.97678 |
| 4 first trace row | 97.678 | 49.430 | 48.248 | 0.0844 | — |

The global rescale (`initialize_simulation`, ≈ line 5620) hits kT = 1 **exactly**:
```c
double actual_ke = kinetic_energy();
double target_ke = particles_active * K_B * temperature_runtime;   /* N, compile-time K_B */
double scale = sqrt(target_ke / actual_ke);
```
All 2.32 % is then lost in `equalize_temperature_per_segment()`, which per segment rescales to `T_target` using `kB_effective()` (**per compartment**, not total) and *then* calls `remove_drift_segment(seg)`:
```c
double T_seg = ke_seg / (c * kB_effective());
double scale = sqrt((double)T_target / T_seg);
... Vx[i] *= scale; Vy[i] *= scale; ...
remove_drift_segment(seg);          /* <- after the normalisation */
```
Subtracting a segment's centre-of-mass velocity removes ⟨½ N_s m v_cm²⟩ = kT per segment on average (2 d.o.f. × ½kT); with two compartments that is ≈ 2 kT of 100, and Px drops from −8.768 to 7×10⁻¹⁵ in exactly that step. The hold conserves energy exactly (97.678111 → 97.678111 under the EDMD integrator, `--mode=edmd`), so nothing is lost afterwards.

**This is an ordering choice — drift removal after temperature normalisation rather than before — not a defect. Reported, not changed.** It means kT is 1 − d/(2N_s) ≈ 0.98 per compartment, and every quantity in §3 uses the measured T_i rather than the nominal 1.

## 3. Level 0 pilot at the requested geometry — everything closes

L0 = **39.25** (pixel grid, see 2a′), H = 10, N = 50/side, r = 0.5, `--edmd-acc=0`, elastic walls, divider mass factor 10⁹, prescribed right piston u = 0.05, gas boundary travels 3.93 σ (78.50 → 74.57; the raw `PistonR_x` moves 4.18 because the piston plane carries a 0.25 σ offset), then hold to t = 283.3 σ-time. 5 seeds, 17000 samples each.

**"Divider held" was approximated, not exact.** The `M ≤ 0` hold-phase branch is not reachable from the CLI in energy-transfer mode (`--wall-mass-factors=0` falls back to the default, and holding via `wall_hold_steps` for the whole run is impossible because `if (!wall_is_released) continue;` gates both the trace and the piston trigger). Mass 10⁹ is the practical equivalent: the divider moved 1.5×10⁻⁵ σ, ΔKE_div ≤ 1.2×10⁻⁵, and KE_gas_left was constant to 5 decimals (49.43047 → 49.43050), so the right compartment is isolated.

η: 0.100051 → 0.111183. KR adiabat over that range gives T_f/T_i = 1.14117, i.e. **W_qs = 7.059·T_i**. The prompt's 6.326·T_i is for η 0.10 → 0.11 exactly; a travel of 3.93 σ is 10 % of the *length*, which is an 11.1 % rise in density, not 10 %. Both are quoted below.

| seed | T_i | T_f | W | ΔKE_right | ΔKE_div | R | \|R\|/W | W − ΣdE | mom (gas) | mom (div) | W_qs | W_diss |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7001 | 0.96495 | 1.12226 | 7.86563 | 7.86560 | 8.0e−06 | −8.1e−12 | 1.0e−12 | −4.9e−11 | −2.4e−11 | −3.1e−05 | 6.8111 | +1.0545 |
| 7002 | 0.99422 | 1.14489 | 7.53363 | 7.53359 | 1.2e−05 | +4.3e−10 | 5.6e−11 | +4.3e−10 | −4.8e−11 | −3.4e−06 | 7.0177 | +0.5159 |
| 7003 | 0.97948 | 1.12109 | 7.08043 | 7.08040 | 5.9e−06 | −5.1e−10 | 7.2e−11 | −4.8e−10 | −5.0e−11 | −1.9e−05 | 6.9137 | +0.1667 |
| 7004 | 0.99115 | 1.13088 | 6.98629 | 6.98626 | 6.2e−06 | +2.9e−10 | 4.1e−11 | +3.0e−10 | −1.5e−11 | +8.2e−06 | 6.9961 | −0.0098 |
| 7005 | 0.95507 | 1.09496 | 6.99462 | 6.99461 | 3.7e−06 | +5.4e−11 | 7.7e−12 | +7.0e−11 | +3.0e−13 | +4.1e−08 | 6.7414 | +0.2532 |

- **Ledger residual R = W − (ΔKE_gas + ΔKE_div): |R|/W ≤ 7×10⁻¹¹.** Round-off.
- **Σ per-event dE vs logged W: ≤ 5×10⁻¹⁰ absolute.**
- **Momentum, gas:** Px_gas(final) − Σdp over all x-events (PL/PR/WL/WR/D*, from t = 0, where Px = 7×10⁻¹⁵) closes to **≤ 5×10⁻¹¹**.
- **Momentum, divider:** ΔPx_div + Σdp over post-release D events closes to 3×10⁻⁵ absolute. Not round-off, but it is *print precision*, not dynamics: the trace writes `W0_v` with `%.6e` and the comparison multiplies by M = 10⁹.
- Both momentum sums require the windows to match: the event log starts at t = 0 while the trace starts at release, and the divider is immobile before release. Aligning them changed the gas closure from 1.4×10⁻² to 1.9×10⁻¹¹.

**⟨W⟩ = 7.2921 ± 0.1752, ⟨W_qs⟩ = 6.8960 ± 0.0530, W_diss = +0.3961 ± 0.1851 kT = +5.7 % of W_qs.**

Within the predicted 0.3–0.9 kT band. With c_s(η = 0.1001) = 1.7443, u/c_s = 0.0287 and the collisionless upper-side estimate is 4.5·u/c_s = **12.9 %** — the measurement sits at 5.7 %, i.e. **below** it, which is what §9.2 of the tex predicts for the collisional regime (collisions let the gas re-equilibrate ahead of the wall). At the earlier, denser trial geometry (η = 0.196) the same comparison came out at 13.7 % against 10–13 %, so the ordering only became clean once the gas was actually dilute.

If instead the prompt's W_qs = 6.326·T_i (η 0.10 → 0.11 exactly) is used, W_qs = 6.181 ± 0.048 and W_diss = +1.111 ± 0.181 kT. The 7.059 coefficient is the one matching the compression actually performed.

## 4. What Level 0 still needs

1. Fix blocker (a), then repeat this pilot at the specified geometry (η = 0.10, L0 = 39.27) — expected W_qs = 6.3 kT, W_diss ≈ 0.3–0.8 kT.
2. Outer-wall impulses in the event log, to close the momentum balance (d).
3. A velocity-reversal switch for (e).
4. Explain the 2 % T_i deficit.
5. The existing `mass_sweep_eta02_N600` reference set was run with `--edmd-acc=1`, the backend that failed validation on 2026-08-26 (181 invalid of 450). Nothing from it is physics until repeated with `--edmd-acc=0`.
