#!/usr/bin/env python3
"""##CHRIS 2026-09-22: Level 3 v4 -- which hypothesis explains the +1.0 to +1.9 kT residual?

v3 left the parameter-free 1-DOF model low by that much at the quasi-static end (2.8-4.2 sigma).
Two candidates, and they make opposite predictions:

  A1  THE GAS. The model drives the wall with the quasi-static force F_ad(L - x_p + s). If the real
      pressure at the far wall during the push is not that, the model is under-driven. Test: hold
      the spring wall rigid, read the momentum delivered to it per collision out of the event log
      (kind D0, column dp), bin it into F_meas(t), and drive the SAME ODE with F_meas instead of
      F_ad. If the residual is gas dynamics, the measured-force ODE lands on the data.

  A2  THE WALL. The model treats the gas as a massless spring, but N m = 100 is comparable to M_s.
      Malek Mansour, Garcia & Baras 2006 Eq. 18 carries M-hat = M + N m / 3 for exactly this.
      Analysis only -- rerun the v3 ODE with M -> M-hat and see if that is enough.

  A3  A SLOW POINT at u = 0.02. If the residual is gas dynamics (A1) the ratio eps_meas/eps_ODE
      must fall toward 1 as the push slows; if it is the wall model (A2) it stays flat, because
      M-hat/M does not depend on u.

Assumption, stated: with the wall held at M = 1e9 it cannot move, so F_meas(t) does not depend on
M_s and one held run per speed serves both M_s of the v3 grid.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

ET = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/"
      "experiments_energy_transfer")
V3 = f"{ET}/level3_master_preload_20260921"
V4 = f"{ET}/level3_v4_20260922"
OUT = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/"
       "paper2_energytransfer/experiments/final")

N, L, H, R = 100, 78.5, 10.0, 0.5
DX, KSPR = 7.96, 0.5
ETA = N * math.pi * R * R / (L * H)
Z = sos.Z_kolafa_rottner_2006(ETA); DZ = sos.dZ_kolafa_rottner_2006(ETA)
CS = math.sqrt(Z + ETA * DZ + Z * Z)
K_GAS = N * CS ** 2 / L ** 2            # adiabatic, = N m c_s^2 / L^2
FGAS = N * Z / L
MHAT_ADD = N / 3.0                      # Mansour-Garcia-Baras: M-hat = M + N m / 3

_LG = np.linspace(L - DX - 10.0, L + 10.0, 2001)
_Zs = np.array([sos.Z_kolafa_rottner_2006(N * math.pi * R * R / (x * H)) for x in _LG])
_lnT = np.concatenate([[0.0], np.cumsum(-0.5 * (_Zs[1:] + _Zs[:-1]) * np.diff(np.log(_LG / L)))])
# ##CHRIS 2026-09-22 BUGFIX: the cumulative sum sets lnT = 0 at the LEFT EDGE of the grid, not at
# L. The adiabat must be normalised T(L) = 1, or every force carries a constant factor
# exp(-Zbar ln(L/L_grid0)) ~ 0.72 -- and since F0 carries it too, the DIFFERENCE F_ad(L') - F0 that
# drives the ODE was ~28 % too weak. This shifted the v3 quasi-static fixed point from 0.840 to
# 0.625 sigma and made the whole v3 model look low.
_lnT = _lnT - np.interp(L, _LG, _lnT)
_FG = N * np.exp(_lnT) * _Zs / _LG
F_ad = lambda Lx: np.interp(Lx, _LG, _FG)
F0 = float(F_ad(L))


def dE(s):
    return FGAS * s + 0.5 * KSPR * s * s


def ode(ms, u, force=None, dt=1e-3, t_end=None, k=KSPR):
    """M s'' = [drive] - k s.  drive = F_ad(L - x_p + s) - F_ad(L), or a measured F(t) - F_ad(L)
    with the gas stiffness restored linearly (the held-wall run cannot know about the wall moving)."""
    tau = DX / u
    t_end = t_end or (tau + 210.0)
    n = int(t_end / dt)
    s = v = 0.0; best = 0.0; tb = 0.0
    for i in range(n):
        t_ = i * dt
        def acc(tt, ss):
            if force is None:
                return (float(F_ad(L - min(u * tt, DX) + ss)) - F0 - k * ss) / ms
            # measured drive; the wall's own recoil is not in F_meas, so add it back as -k_gas*s
            return (float(np.interp(tt, force[0], force[1])) - F0 - K_GAS * ss - k * ss) / ms
        k1v = acc(t_, s); k1s = v
        k2v = acc(t_ + dt/2, s + dt/2*k1s); k2s = v + dt/2*k1v
        k3v = acc(t_ + dt/2, s + dt/2*k2s); k3s = v + dt/2*k2v
        k4v = acc(t_ + dt, s + dt*k3s);     k4s = v + dt*k3v
        s += dt/6*(k1s + 2*k2s + 2*k3s + k4s); v += dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        e = dE(s)
        if e > best: best, tb = e, t_
    return best, tb - tau


def measured_force(u, dtbin=0.5):
    """F_meas(t) on the held spring wall: sum of |dp| per bin / bin width, averaged over seeds."""
    d = f"{V4}/A1_held_u{u}"
    files = sorted(glob.glob(f"{d}/ev_*.csv"))
    tmax = 0.0
    per = []
    for f in files:
        e = pd.read_csv(f)
        e = e[e["kind"] == "D0"]
        if len(e) == 0:
            continue
        per.append((e["t_sigma"].to_numpy(float), np.abs(e["dp"].to_numpy(float))))
        tmax = max(tmax, per[-1][0].max())
    nb = int(tmax / dtbin) + 1
    acc = np.zeros(nb)
    for t, dp in per:
        idx = (t / dtbin).astype(int)
        np.add.at(acc, idx, dp)
    F = acc / (dtbin * len(per))
    tc = (np.arange(nb) + 0.5) * dtbin
    return tc, F, len(per)


def cell_meas(run, name, ms):
    """dE_coh from a v3/v4 cell, baseline and peak-bias both from the causally blind window."""
    tails, blinds = [], []
    tt = None
    for tr in sorted(glob.glob(f"{run}/{name}/tr_*.csv")):
        t = pd.read_csv(tr, low_memory=False)
        time = t["Time"].to_numpy(float); E = t["SpringE"].to_numpy(float)
        v = np.abs(t["PistonR_v"].to_numpy(float)); mv = np.nonzero(v > 1e-12)[0]
        if len(mv) < 2: continue
        i0, i1 = mv[0], mv[-1]
        bl = np.nonzero(time - time[i0] <= 0.4 * L / CS)[0]; bl = bl[bl >= i0]
        if len(bl) < 100: continue
        b = float(E[bl].mean())
        tails.append(E[i1:] - b); blinds.append(E[bl] - b)
        ti = time[i1:] - time[i1]
        if tt is None or len(ti) < len(tt): tt = ti
    if not tails: return None
    n = min(len(x) for x in tails); nb = min(len(x) for x in blinds)
    Eb = np.mean([x[:n] for x in tails], axis=0)
    Es = np.std([x[:n] for x in tails], axis=0, ddof=1) / math.sqrt(len(tails))
    bias = float(np.max(np.mean([x[:nb] for x in blinds], axis=0)))
    j = int(np.argmax(Eb))
    return dict(dE=float(Eb[j]) - bias, ddE=float(Es[j]), n=len(tails), bias=bias)


def main():
    print("## Level 3 v4 -- what is the residual?\n")
    print(f"k_gas^ad = N m c_s^2/L^2 = {K_GAS:.5f}; F = {FGAS:.4f} kT/sigma; "
          f"M-hat adds N m/3 = {MHAT_ADD:.2f}\n")

    print("### A1 -- drive the ODE with the MEASURED force on a held wall\n")
    fm = {}
    for u in (0.05, 0.1):
        tc, F, nseed = measured_force(u)
        fm[u] = (tc, F)
        tau = DX / u
        pre = F[tc < 0.4 * L / CS].mean()
        during = F[(tc > tau * 0.5) & (tc < tau)].mean()
        after = F[(tc > tau + 50) & (tc < tau + 200)].mean()
        print(f"  u = {u:g}: {nseed} seeds; F_meas before the signal arrives {pre:.4f} "
              f"(static F = {FGAS:.4f}), late in the push {during:.4f}, after the stop {after:.4f}")
    print()

    print("| u | M_s | measured ΔE_coh | ODE v3 (F_ad) | ODE A1 (F_meas) | ODE A2 (M-hat) | "
          "|meas−A1|/σ | |meas−A2|/σ | |meas−v3|/σ |")
    print("|---|---|---|---|---|---|---|---|---|")
    rows = []
    for u in (0.05, 0.1):
        for ms in (50, 200):
            c = cell_meas(V3, f"k0.5_M{ms}_u{u}", ms)
            if not c: continue
            e_v3, _ = ode(ms, u)
            e_a1, _ = ode(ms, u, force=fm[u])
            e_a2, _ = ode(ms + MHAT_ADD, u)
            s1 = abs(c["dE"] - e_a1) / c["ddE"]; s2 = abs(c["dE"] - e_a2) / c["ddE"]
            s3 = abs(c["dE"] - e_v3) / c["ddE"]
            rows.append((u, ms, c, e_v3, e_a1, e_a2, s1, s2, s3))
            print(f"| {u:g} | {ms} | {c['dE']:.3f} ± {c['ddE']:.3f} | {e_v3:.3f} | {e_a1:.3f} | "
                  f"{e_a2:.3f} | {s1:.1f} | {s2:.1f} | {s3:.1f} |")

    print("\n### A3 -- the slow point: does the ratio fall toward 1?\n")
    print("| u | τ_push/(L/c_s) | M_s | measured ΔE_coh | ODE v3 | ratio meas/ODE |")
    print("|---|---|---|---|---|---|")
    for u, run, pat in ((0.02, V4, "A3_M{}_u0.02"), (0.05, V3, "k0.5_M{}_u0.05"),
                        (0.1, V3, "k0.5_M{}_u0.1")):
        for ms in (50, 200):
            c = cell_meas(run, pat.format(ms), ms)
            if not c: continue
            e_v3, _ = ode(ms, u, t_end=DX / u + 210.0)
            print(f"| {u:g} | {(DX/u)/(L/CS):.2f} | {ms} | {c['dE']:.3f} ± {c['ddE']:.3f} | "
                  f"{e_v3:.3f} | {c['dE']/e_v3:.2f} |")

    print("\n### Criterion (iv): which hypothesis do the data select?\n")
    if rows:
        m1 = max(r[6] for r in rows); m2 = max(r[7] for r in rows); m3 = max(r[8] for r in rows)
        print(f"  worst |meas − ODE|/σ:  A1 measured force {m1:.1f}σ | A2 M-hat {m2:.1f}σ | "
              f"v3 baseline {m3:.1f}σ")
        best = min((m1, "A1 (gas dynamics)"), (m2, "A2 (gas inertia)"), (m3, "neither (v3)"))
        print(f"  A1 passes (<1σ at all four)" if m1 < 1 else f"  A1 does NOT pass (<1σ required)")
        print(f"  closest hypothesis: {best[1]} at {best[0]:.1f}σ")

    # figure: measured vs quasi-static force
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    for u, col in ((0.05, "tab:blue"), (0.1, "tab:red")):
        tc, F = fm[u]
        ax.plot(tc, F, lw=1.0, color=col, alpha=0.85, label=f"measured, u = {u:g}")
        tau = DX / u
        tq = np.linspace(0, tc.max(), 600)
        ax.plot(tq, [float(F_ad(L - min(u * x, DX))) for x in tq], "--", color=col, lw=1.4,
                label=f"quasi-static $F_{{ad}}$, u = {u:g}")
    ax.axhline(FGAS, color="0.5", ls=":", lw=1.0)
    ax.text(5, FGAS * 1.02, "standing force $F = NkTZ/L$", fontsize=8, color="0.4")
    ax.set_xlabel("time from the start of the push  [σ-time]")
    ax.set_ylabel("force on the held wall  [$k_BT/\\sigma$]")
    ax.set_xlim(0, 400); ax.grid(True, ls=":", alpha=0.6); ax.legend(fontsize=8.5)
    ax.set_title("A1: what the gas actually delivers to the far wall", fontsize=11)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/260922_level3_v4_measured_force.{ext}", dpi=200)
    print(f"\nfigure -> {OUT}/260922_level3_v4_measured_force.png/.pdf")


if __name__ == "__main__":
    main()
