#!/usr/bin/env python3
"""##CHRIS 2026-10-13: mode-equivalence gate (261012 sec. 1.9 C2), run on the Mac. One test cell, the pi/8 anchor
(H = L_0 = 10, N_s = 50, t = 0.05), in both modes; every quantity READ FROM EACH MODE'S OWN OUTPUT where the mode
writes it. Rule: equal to 1e-6 (or to the print precision where that is coarser -- stated).

Inputs: speed-of-sound = the Mac pi/8 pilot, m_50 run0 trace (confinement_pilot_20261013/mac_pi8_H10_L10);
energy-transfer = the held-divider gate run (paper1_confinement_modegate_20261013/epi8_H_H10_L10/x_0, seed 9700).

Code lines that write or set each quantity (00ALLINONE.c unless stated):
  box        322-327  SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA); XW2 = XW1 + SIM_WIDTH;  (XW1 = 200 px, line 315)
  t          281-283  set_wall_thickness_sigma(): wall_thickness_runtime = thickness_in_sigma * PIXELS_PER_SIGMA;
             4787-4788 in parse_cli_options() (3481), called once at 20751 before either experiment runs
  eta (SoS)  15800-15802  eta_nominal_const = N pi r^2 / (2 L0_UNITS H)   -> trace column "eta"
  eta (ET)   17289-17290  eta_nominal = N pi r^2 / (2 L0 H)                -> summary column "eta_nominal"
  t (ET)     17295        t_sigma = fabs(wall_thickness_runtime / PIXELS_PER_SIGMA) -> summary "wall_thickness_sigma"
  segments   3214-3221  compute_segment_bounds(): segment s = [prev, wall - t/2], next starts at wall + t/2 -> SegEtas
  Center_X   15797      center_x_sigma_const = (XW1 + XW2) / (2 PIXELS_PER_SIGMA)  (SoS trace)
"""
import glob, math, os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
SOS = glob.glob(os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/"
                             "confinement_pilot_20261013/mac_pi8_H10_L10/m_50/wall_x_positions_*_run0.csv"))[0]
ET = os.path.join(HS, "experiments_energy_transfer/paper1_confinement_modegate_20261013/epi8_H_H10_L10/x_0")
T_IN, R_IN = 0.05, 0.5                                      # command-line inputs of both runs

def main():
    s = pd.read_csv(SOS, nrows=1).iloc[0]
    es = pd.read_csv(os.path.join(ET, "summary_9700.csv")).iloc[0]
    et = pd.read_csv(os.path.join(ET, "tr_9700.csv"), nrows=1, low_memory=False).iloc[0]
    L0s, L0e = float(s["L0"]), float(es["L0"])
    # The SoS trace has no pre-release row (its first row is one step after release, Displacement -0.000251 at M = 50),
    # so the HELD position is read from run.log, "Initial wall_x = %.3f" in px (line 15789), minus XW1 = 200 px (line 315).
    # A first version compared the post-release first row and reported a spurious 2.5e-4 FAIL; corrected, disclosed in sec. 1.10.
    import re
    w0 = float(re.search(r"Initial wall_x = ([\d.]+)", open(os.path.join(os.path.dirname(SOS), "run.log")).read()).group(1))
    xs = (w0 - 200.0) / 24.0                                            # divider centre from the left wall, +- 2.1e-5
    xe = float(et["W0_x_sigma"]); te = float(es["wall_thickness_sigma"])
    Ns, Ne = int(s["Left_Count"]) + int(s["Right_Count"]), int(es["particles_total"])
    Hs = Ns * math.pi * R_IN ** 2 / (2 * L0s * float(s["eta"]))          # SoS does not write H: from its eta
    segL = lambda x, t, L0: (x - t / 2, 2 * L0 - x - t / 2)
    Ls, Le = segL(xs, T_IN, L0s), segL(xe, te, L0e)
    cnt = [int(v) for v in str(et["SegCounts"]).split(";")]; seta = [float(v) for v in str(et["SegEtas"]).split(";")]
    Lseg = [n * math.pi * float(es["radius_sigma"]) ** 2 / (float(es["height"]) * e) for n, e in zip(cnt, seta)]
    rows = [("eta (nominal)", float(s["eta"]), float(es["eta_nominal"]), 1e-6, "both written, %.6f"),
            ("L_0", L0s, L0e, 1e-6, "both written"),
            ("N", Ns, Ne, 0, "SoS: Left+Right counts; ET: particles_total"),
            ("H", Hs, float(es["height"]), 2e-5, "SoS does not write H; inferred from its eta (6-decimal print -> ~1e-5)"),
            ("divider centre from left wall (held)", xs, xe, 2.1e-5, "SoS: run.log Initial wall_x (px, 3 dec. -> 2.1e-5); ET: W0_x_sigma"),
            ("t", T_IN, te, 1e-6, "SoS does NOT write t (input shown); ET summary"),
            ("left free length", Ls[0], Le[0], 2.1e-5, "x - t/2 (inherits the 2.1e-5 of x)"),
            ("right free length", Ls[1], Le[1], 2.1e-5, "2 L_0 - x - t/2"),
            ("L_eff = free length - 2r", Ls[0] - 2 * R_IN, Le[0] - 2 * float(es["radius_sigma"]), 2.1e-5, "SoS r not written (input)")]
    print("### Mode-equivalence gate, pi/8 anchor (H = L_0 = 10, N_s = 50, t = 0.05)\n")
    print("| quantity | speed-of-sound | energy-transfer | abs. difference | tolerance | verdict | note |")
    print("|---|---|---|---|---|---|---|")
    ok_all = True
    for q, a, b, tol, note in rows:
        d = abs(a - b); ok = d <= tol; ok_all &= ok
        print(f"| {q} | {a:.6f} | {b:.6f} | {d:.1e} | {tol:.0e} | {'PASS' if ok else 'FAIL'} | {note} |")
    print(f"\nSegEtas cross-check (ET honours t): free length from SegEtas = {Lseg[0]:.5f} / {Lseg[1]:.5f} "
          f"vs x - t/2 = {Le[0]:.5f} (SegEtas printed to 6 decimals -> ~1e-5 sigma); with t ignored it would be {L0e:.5f}.")
    print("\nGrid exactness of every campaign geometry (the (int) cast at line 323 truncates 2 L_0 x 24 px):")
    sys.path.insert(0, HERE)
    import contextlib, io, paper1_confinement_prereg_20261012 as PR
    with contextlib.redirect_stdout(io.StringIO()): C, _ = PR.main()
    bad = [(c["L0"], c["H"]) for c in C if abs(2 * c["L0"] * 24 - round(2 * c["L0"] * 24)) > 1e-9 or abs(c["H"] * 24 - round(c["H"] * 24)) > 1e-9]
    print(f"  {len(C)} cells; 2 L_0 x 24 and H x 24 integer in all: {'YES' if not bad else 'NO ' + str(bad)}")
    print("\nOPEN, outside this campaign: the same (int) cast on the canonical A1v2 cells with non-grid L_0 (260919 table):\n")
    print("| eta | L_0 (table) | 2 L_0 x 24 px | box after (int) | box shortened by [sigma] | relative |")
    print("|---|---|---|---|---|---|")
    import csv, tests_20260913 as TT
    for r in csv.DictReader(open(TT.plot_path('260919_A1v2_final_cs_vs_eta.csv'))):
        L0 = float(r['L0']); px = 2 * L0 * 24
        if abs(px - round(px)) > 1e-6:
            cut = (px - int(px)) / 24
            print(f"| {float(r['eta']):.6f} | {L0} | {px:.4f} | {int(px)} | {cut:.4f} | {cut/(2*L0):.1e} |")
    print(f"\n**GATE: {'PASS' if ok_all else 'FAIL'}** -- every quantity both modes write agrees to its tolerance. t and r are not "
          "written by speed-of-sound mode; they are equal by construction (one global, set in parse_cli_options before "
          "either experiment runs), and H is pinned by the equal eta at fixed N, r, L_0.")

if __name__ == "__main__":
    main()
