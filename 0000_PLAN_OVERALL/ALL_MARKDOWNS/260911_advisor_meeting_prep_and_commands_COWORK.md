# 260911 — Advisor meeting prep (speed of sound) + the terminal commands

For Chris, for the meeting on 2026-09-12. Read the first two sections on paper; the commands are for the terminal in `hspist3/`. Everything here is what is on disk today; where I could not verify a flag by running it (no Mac shell from here), it says so.

---

## 1. What we do, how, and why — in one page

**The measurement.** Two gas compartments of 50 hard disks each (radius 0.5, box height 10), separated by a movable wall (the divider, mass M). Release the divider: the two gas columns act as springs and it oscillates. The oscillation frequency ν depends on the stiffness of the gas, which is set by the speed of sound. Román et al. (2002) derived the relation for exactly this geometry:

    ν = c_s · K / (2π L_eff),      cot K = (M / 2Nm) · K,      L_eff = L0 − 2r

K is a number between 0 and π/2 fixed by the mass ratio (heavy divider → small K, light → K → π/2). So one run gives one ν; nine divider masses give nine (ν, K) points that must lie on one straight line through the origin, and its slope is c_s. That is the instrument.

**The two "25"s (your question).** They are different things. `--repeats=25`: each (η, M) cell is run 25 times with independent seeds, so every ν has an error bar (the standard error over the 25). `--target-oscillations=25`: every single run is made long enough that the divider completes ≥ 25 oscillations, so the FFT peak is sharp enough to read ν. So the campaign is 32 densities × 9 masses × 25 repeats = 7200 trajectories, each ≥ 25 periods long.

**Why the density sweep is done by changing L0 (route A).** With N = 100 fixed, η = 3.927/L0: η = 0.02 means a compartment 196 diameters long, η = 0.65 means 6 diameters. That gives 32 densities from the ideal-gas end to the dense fluid with one code path and one particle number. Its known weakness: box size and density change together, so it validates the *instrument* against the equation of state (EOS), not the bulk fluid. The bulk claim needs the other sweep (below).

**What "compressibility" and "pressure" mean here** (the χ² PDF §1 has the slow version). Pressure P is the momentum the disks deliver to a wall per unit time per unit wall length. We report it as Z = PA/(N k_BT): Z = 1 is the ideal gas, Z > 1 because disks take up room (Z ≈ 8.4 at η = 0.65). The name "compressibility factor" for Z is historical and confusing — Z is a pressure ratio. The sound speed is the *change* of pressure with density at constant entropy, c_s² = (∂P/∂ρ)_S, and for hard disks that is exactly

    c_s² = (k_BT/m) [ Z + η Z′(η) + Z² ],     ideal gas: c_s = √2 (γ = 2 in 2D),     first order: c_s = √2 (1 + 2η).

So the same EOS Z(η) (Kolafa–Rottner) predicts both the pressure campaign and the sound-speed campaign. That is why the two together validate the simulator: one equation, two independent measurements.

**Your three proposed sweeps — do they make sense?**

| sweep | what changes | what stays fixed | what it tests | status |
|---|---|---|---|---|
| A. L0 varies, N = 50/50, r = 0.5 | density and box size together | N, Kn ≈ 0.14 | the instrument vs the EOS over the whole fluid | done (32 η, 25 repeats) |
| B. N = 500/500 with smaller r at the same η and the same L0 | box size in units of σ (grows 3.2×), N | η, L0 in absolute units | is c_s intensive, i.e. independent of box size at fixed η — the finite-size question | partly: route B (r sweep at N = 100, L0 = 20) and the N-ladders exist; the clean version is famB (N = 100…1600 at fixed η and fixed aspect ratio, r = 0.5), launched today |
| C. N per side from 1 upward at fixed L0, r = 0.5 (then again with smaller r) | density AND collisions per crossing | box | at large N: a density sweep at fixed geometry; at small N (< ~10 per side at L0 = 20) the gas stops being a fluid, γ → 3 and c_s → √3 = 1.73 — the Knudsen crossover | not run; the small-N end is the positive control we planned (N_side = 25, 10 at L0 = 200) |

Sweep C is sound, with one warning to say out loud: below ~10 disks per side the divider no longer measures a fluid's sound speed, it measures single-particle bouncing (Lua–Grosberg piston), so those points belong to a different curve. Above that it is a legitimate fixed-box density sweep. The order for the paper is A (done) → famB (bulk claim, running) → C as the control that shows the method knows the difference between a fluid and a few bouncing particles.

**The finding of this week (say this to Susanne).** The seeding rescaled each compartment to k_BT = 1 and *then* removed its centre-of-mass velocity, which throws away on average one k_BT per compartment: every run actually started at T_i = 1 − X/50 with X exponentially distributed around 1 (mean 0.98, spread 0.02). Since c_s ∝ √T, every c_s was ≈ 1 % low with 1 % run-to-run scatter. Reconstructed from the seeds for all 7200 runs and corrected: the dilute end moves from −0.2 % to +0.8 % above the EOS, the mid range from +1.1 % to +2.2 %. My reading: route A holds the Knudsen number (mean free path / compartment) fixed at 0.14 at every η, and a gas at fixed finite Knudsen has a sound speed a fraction of a percent above the hydrodynamic value (effective γ creeps from 2 toward 3). That makes it a prediction: the offset must shrink as N grows at fixed η (famB) and grow as N shrinks (sweep C). New runs use `--seed-drift-order=drift-first` (T_i = 1.000000 exactly); old runs are corrected in analysis and stay bit-identical on disk. Pressure was never affected (its runner measures T directly).

**Numbers to have ready.** Pressure: validated η ≤ 0.65 (four box sizes 400–2500 disks, −0.35 ± 0.06 % at 0.65), no bulk value at 0.67–0.69 (the fluid begins to order; the finite-size law breaks). Sound: instrument validated at the ≤ 1 % level for η ≤ 0.4 after the temperature correction; bulk numbers come from famB (5 densities × 4 sizes, running). Paper 2: energy ledger of a piston push closes to 10⁻¹¹; dissipation 5–10 % of the reversible work at u/c_s = 0.01–0.06; one open offset being diagnosed.

**Before the meeting, do on paper:** hand-calculation PDF §4 (c_s from the EOS) and §5 (the Román relation, including the two limits), χ² PDF §1–2 (pressure and the box). One hour.

---

## 2. Where the plot code is

| file | what | edit live? |
|---|---|---|
| `0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots/plot_cs_final_style.py` | the two current figures (`260910_FINAL_cs_vs_eta.pdf`, `260910_FINAL_cs_idealgas_zoom.pdf`): shaded regimes, SPT/Henderson/KR, Román points, our points | **yes — simplest**: 110 lines, no dependencies beyond numpy/matplotlib, reads one CSV |
| `260909_plots/routeA_refit_cs_vs_eta_20260909.csv` | the data the figure reads (eta, c_s, c_s_err, …). CC will write a temperature-corrected version today; point the script at that file when it exists | — |
| `hspist3/plot_speed_of_sound_edmd.py` | the campaign's own FINAL plotter (regions, EOS curves, Román points, `k_root_bisect`, weighted fits); called by the analysis pipeline | possible, but 85 kB |
| `hspist3/replot_final_from_summary.py` | redraws the campaign FINAL figure from `final_plots/combined_speed_of_sound_summary.csv` without re-reading traces (seconds instead of minutes) | yes |
| `hspist3/analyze_speed_of_sound_by_eta.py` | the full pipeline: traces → FFT → ν per run → fits → summary CSV → figures | no (slow) |

Regenerate the current figures (seconds):

    cd "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots"
    python3 plot_cs_final_style.py routeA_refit_cs_vs_eta_20260909.csv 260910_FINAL

To change things live: the regime colours are in `REGIONS` (line ~37), the EOS curves in `eos_curves()`, the y-range and labels in the figure blocks. Change, re-run the command, open the PDF.

---

## 3. Commands

All from `cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3`. Build first if in doubt: `make clean && make release` (the debug build is 5× slower).

**3a. Watch one simulation in the SDL window (GUI).** One divider, 50/50 disks, L0 = 20 per compartment (η = 0.196), divider mass 200× particle mass:

    ./00ALLINONE --no-experiments --mode=edmd --particles=100 --particles-boxes=50,50 \
      --particle-radius=0.5 --l0=20 --height=10 --num-walls=1 --wall-mass-factors=200 --kbt1

(`--no-experiments` forces the interactive window. If a flag is rejected, `./00ALLINONE --help | grep -i wall` lists the current names; CC can confirm the exact GUI line in one minute.) **About the window size:** the physics scale is 24 px per diameter and the window is sized from the box, then *clamped to the display* (`00ALLINONE.c` ~4576). The view scale in `gui_compute_scene_view_scale()` only enlarges (1.0–2.6×), it never shrinks, so a very long box (η = 0.02 → 2 × 196 σ = 9400 px) runs out of the window. For the demo stay at L0 ≤ ~35 (η ≥ 0.11); the right-hand panel (sliders, live log) and the histograms below are laid out from the remaining space.

**3b. One headless cell, then plot the divider position vs time.** This is exactly one accepted campaign cell (η = 0.262, L0 = 15) with the lightest divider, 3 repeats, into a scratch folder:

    mkdir -p /tmp/demo_cs
    ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
      --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 \
      --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=15.0 --wall-masses=50 \
      --repeats=3 --seed=2026082309 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 \
      --seed-drift-order=drift-first \
      --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 \
      --oscillation-max-steps=10000000 --speed-sound-log-stride=auto --speed-sound-run-dir=/tmp/demo_cs

Output: `/tmp/demo_cs/wall_x_positions_L0_150_wallmassfactor_50_run{0,1,2}.csv` with columns `Time, Wall_X, Displacement(σ), Left_Count, Right_Count, L0, eta, Center_X(σ), Seed, …`. Plot the three repeats:

    python3 - <<'EOF'
    import glob, csv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10,4))
    for f in sorted(glob.glob("/tmp/demo_cs/wall_x_positions_*_run*.csv")):
        rows = list(csv.DictReader(open(f)))
        t = [float(r["Time"]) for r in rows]; x = [float(r["Displacement(σ)"]) for r in rows]
        ax.plot(t, x, lw=0.8, label=f.split("_")[-1].replace(".csv",""))
    ax.set_xlabel("time (σ-time)"); ax.set_ylabel("divider displacement (σ)")
    ax.set_title("η = 0.262, L0 = 15, M = 50 m — three repeats"); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig("/tmp/demo_cs/divider_x_vs_time.png", dpi=150); print("wrote /tmp/demo_cs/divider_x_vs_time.png")
    EOF

(If the column header differs slightly, `head -1` the CSV and adjust the two names. The FFT of exactly this signal is what gives ν; `wall_x_FFT.py` in `hspist3/` is the old standalone version.)

**3c. The whole route-A pipeline yourself.** The campaign script does everything (one folder per η, 9 masses, seeds, health checks) and prints the analysis command at the end:

    ./run_eta_scan_campaign.sh OUTDIR REPEATS JOBS ETA_SET NPART
    # example — quick clean run, 5 repeats, all 32 densities (0.0196 … 0.76), 50/50 disks:
    ./run_eta_scan_campaign.sh experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/campaign_r5_demo_$(date +%Y%m%d) 5 6 all 100

`ETA_SET`: `roman` (12 dilute/mid densities), `dense` (0.55–0.76), `all`, `transition`, `ladder`. REPEATS 5 is ~1/5 of the 25-repeat campaign; that one took the machine the better part of a day, so start with `roman` if you want a result during the meeting. Two things to know before a *new* run: (i) the script does not yet pass `--seed-drift-order=drift-first` — ask CC to add it to `run_one()` (one line) so new runs start at T = 1 exactly; (ii) η = 0.005 needs L0 = 785 σ, which hits the float seed-pad issue at L0 ≥ 100 — stay at η ≥ 0.0196 until the pad is fixed. When the script finishes it prints the `analyze_speed_of_sound_by_eta.py … --write-final --roman-ref --theory-cs adiabatic` command; run it, and the FINAL figure appears in `OUTDIR/analysis/final_plots/`. To redraw that figure quickly afterwards:

    python3 replot_final_from_summary.py --summary OUTDIR/analysis/final_plots/combined_speed_of_sound_summary.csv --roman-ref

**3d. The fixed-geometry family (bulk claim).** `./run_finitesize_aspect_campaign.sh` (famA/famB runner); CC launched famB today with η = 0.10/0.30/0.50/0.60/0.65 and N = 100/400/900/1600. Do not start another one while it runs.

---

## 4. If Susanne asks…
- "Why route A first?" — one particle number, one code path, 32 densities; it is the instrument test. The bulk numbers come from famB.
- "How good is the ideal-gas end?" — slope c_s/√2 − 1 = 0.035 vs 2η = 0.039 at η = 0.02; after the temperature correction the value sits +0.8 % above the EOS, which we interpret as the finite-Knudsen offset of a fixed-Kn sweep; famB tests that.
- "Error bars?" — 25 repeats per mass → SEM per ν; nine masses → weighted line through the origin; the quoted error is the fit error; per-mass scatter 0.5–0.7 % is FFT-bin resolution (3 % per run), being replaced by a damped-cosine fit.
- "And above 0.65?" — pressure: no bulk value (finite-size law breaks, χ² = 34); sound: 6-diameter compartments, structure forms, not a fluid-branch value. Both papers stop the fluid claim at 0.65.
