#!/usr/bin/env python3
"""##CHRIS 2026-09-18: the compartment is 0.5 sigma shorter than every analysis so far assumed.

Found while checking the ramp protocol's travel. Three facts, each measured, not argued:

  * wall_thickness_sigma = 1 in every energy-transfer run (summary.csv), and the wall position
    flag sets the divider CENTRE, so the gas-side face is at 39.25 + 0.5 = 39.75.
  * the stop snapshots confirm it: the closest disk centre to the divider over 60 runs is 40.2867,
    and a disk of radius 0.5 cannot come closer than 40.25 to a face at 39.75.
  * the gas's right boundary during the hold is the box wall at 78.50 (the travel flag is measured
    from it: target = 78.50 - 3.93 = 74.57, which is where the piston stops). The piston parks at
    78.75, so its first 0.25 sigma of travel happens outside the gas: the piston displacement is
    4.18 sigma but the gas is compressed by exactly 3.93 sigma.

So the true right compartment runs 39.75 -> 78.50, i.e. L = 38.75, not the nominal L0 = 39.25 used
in paper2_level0_level1_20260916.py. That file's zwall() uses L twice -- once to convert the wall
impulse into Z (Z = sum|dp| L / (T N_s)) and once to set eta -- so both move:

    Z_true = Z_reported * (L_true / L_nominal),      eta_true = eta_reported * (L_nominal / L_true).

This recomputes the Level 1 quasi-static work on the corrected geometry. Analysis only, no re-runs:
the same event logs are read, only the lengths change.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

ET = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_energy_transfer"
DT = 1.0 / 60.0
NS, H, R = 50, 10.0, 0.5
WALL_THICK = 1.0                      # summary.csv: wall_thickness_sigma
FACE = 0.5 * WALL_THICK               # divider centre -> gas-side face
AREA_DISKS = NS * math.pi * R * R


def zwall(files, L, hold_steps=12000):
    """Z from the wall impulse over the hold, with L the TRUE face-to-face compartment length."""
    T = hold_steps * DT
    z = []
    for ev in files:
        e = pd.read_csv(ev)
        h = e[e["t_sigma"] < T]
        zr = h[h["kind"] == "WR"]["dp"].abs().sum() * L / (T * NS)
        zd = h[(h["kind"] == "D0") & (h["dp"] > 0)]["dp"].sum() * L / (T * NS)
        z.append(0.5 * (zr + zd))
    z = np.array(z)
    return z.mean(), z.std(ddof=1) / math.sqrt(len(z)), len(z)


def eta_of(L):
    return AREA_DISKS / (L * H)


def w_qs(points, L_i, L_f):
    """W_qs = N_s kT (T_f/T_i - 1) with ln(T_f/T_i) = int Z dln(eta) along the measured path."""
    pts = sorted(points)                       # (L_true, Z, sem)
    lne = np.array([math.log(eta_of(p[0])) for p in pts])
    Z = np.array([p[1] for p in pts])
    sZ = np.array([p[2] for p in pts])
    o = np.argsort(lne); lne, Z, sZ = lne[o], Z[o], sZ[o]
    # straight line in ln eta, weighted; the three/four points are collinear to well under their errors
    w = 1.0 / sZ ** 2
    S, Sx, Sxx = w.sum(), (w * lne).sum(), (w * lne * lne).sum()
    Sy, Sxy = (w * Z).sum(), (w * lne * Z).sum()
    D = S * Sxx - Sx ** 2
    a, b = (Sxx * Sy - Sx * Sxy) / D, (S * Sxy - Sx * Sy) / D
    chi2 = float((w * (Z - a - b * lne) ** 2).sum()) / max(1, len(lne) - 2)
    x0, x1 = math.log(eta_of(L_i)), math.log(eta_of(L_f))
    integ = a * (x1 - x0) + 0.5 * b * (x1 ** 2 - x0 ** 2)
    # Error through the SAME integral, with the full covariance. Var(a) and Var(b) are separately
    # huge because the intercept sits at ln eta = 0, far outside the data; their covariance is
    # correspondingly negative and cancels almost all of it. Dropping it inflates the error ~100x.
    va, vb, cab = Sxx / D, S / D, -Sx / D
    ca, cb = (x1 - x0), 0.5 * (x1 ** 2 - x0 ** 2)
    dint = math.sqrt(max(0.0, ca * ca * va + cb * cb * vb + 2.0 * ca * cb * cab))
    W = NS * (math.exp(integ) - 1.0)
    dW = NS * math.exp(integ) * dint * math.sqrt(max(1.0, chi2))
    return W, dW, integ, chi2


def main():
    print("### The geometry correction\n")
    print(f"divider centre 39.25, thickness {WALL_THICK} σ  ->  gas-side face 39.75")
    print("box right wall 78.50 (the travel flag measures from it), piston parks at 78.75\n")
    print("| | nominal L0 used before | true face-to-face L | η before | η true |")
    print("|---|---|---|---|---|")
    for nom in (39.25, 38.270833, 37.291667, 36.3125, 35.3125):
        Lt = nom - FACE
        print(f"| hold point | {nom:.6f} | {Lt:.6f} | {eta_of(nom):.6f} | {eta_of(Lt):.6f} |")
    print(f"| push start | 39.25 | {78.50 - 39.75:.6f} | {eta_of(39.25):.6f} | {eta_of(78.50 - 39.75):.6f} |")
    print(f"| push end | 35.3125 | {74.57 - 39.75:.6f} | {eta_of(35.3125):.6f} | {eta_of(74.57 - 39.75):.6f} |")
    print(f"\ngas compression Δx = 78.50 − 74.57 = {78.50 - 74.57:.2f} σ (the flag value; the piston itself moves 4.18)")

    print("\n### Z along the path, on the true lengths\n")
    print("| true L | η | seeds | Z_wall | sem | Z_KR | Z_wall/Z_KR |")
    print("|---|---|---|---|---|---|---|")
    pts = []
    fam = [(39.25, sorted(glob.glob(os.path.join(ET, "level0_Wqs_20260911", "u*", "ev_*.csv"))))]
    for nom, tag in ((38.270833, "L38p270833"), (37.291667, "L37p291667"),
                     (36.3125, "L36p3125"), (35.3125, "L35p3125")):
        f = sorted(glob.glob(os.path.join(ET, "level1_Zwall_path_20260916", tag, "ev_*.csv")))
        if f:
            fam.append((nom, f))
    for nom, files in fam:
        Lt = nom - FACE
        Z, s, n = zwall(files, Lt)
        e = eta_of(Lt); ZK = sos.Z_kolafa_rottner_2006(e)
        pts.append((Lt, Z, s))
        print(f"| {Lt:.6f} | {e:.6f} | {n} | {Z:.4f} | {s:.4f} | {ZK:.4f} | {Z / ZK:.4f} |")

    L_i, L_f = 78.50 - 39.75, 74.57 - 39.75
    W, dW, integ, chi2 = w_qs(pts, L_i, L_f)
    print(f"\npath: L {L_i:.4f} → {L_f:.4f} σ, η {eta_of(L_i):.6f} → {eta_of(L_f):.6f}")
    print(f"∫Z dlnη = {integ:.6f}  →  T_f/T_i = {math.exp(integ):.6f}   (fit χ²/dof = {chi2:.2f})")
    print(f"W_qs^finite, corrected geometry = {W:.4f} ± {dW:.4f} kT")
    print(f"W_qs^finite, as published 2026-09-17 = 7.4879 ± 0.0066 kT")

    print("\n### What that does to Level 1 and Level 2\n")
    W0, dW0 = 7.4728, 0.0089      # u -> 0 intercept of the quadratic fit, 760 trajectories
    g = W0 - W; sg = math.hypot(dW0, dW)
    print(f"Level 1: W(0) − W_qs = {g:+.4f} ± {sg:.4f} kT  ({abs(g) / sg:.1f} σ)  "
          f"{'PASS' if abs(g) < 2 * sg else 'FAIL'}")
    print(f"  the published version had {7.4728 - 7.4879:+.4f} ± 0.0110, i.e. the measured work sat BELOW")
    print("  the quasi-static bound, which is not allowed for an adiabatic compression; on the corrected")
    print(f"  geometry the sign is {'physical (W ≥ W_qs)' if g > 0 else 'still inverted'}.")
    for nm, dx in (("gas Δx", 3.93), ("piston displacement", 4.18)):
        print(f"  ζ₀Δx with {nm} = {dx}: {2.615 * dx:.2f};  ζ_hyd Δx = {0.084 * dx:.3f}")


if __name__ == "__main__":
    main()
