#!/usr/bin/env python3
"""##CHRIS 2026-09-21: Level 3 v3 -- the one-degree-of-freedom model, no free parameter.

What the spring is doing is NOT impedance matching. Look at where the peaks sit: the maximum of
eps_coh tracks tau_push/T_w, not Z_w/Z_g. That is a driven oscillator: a pressure rise applied
faster than the wall's own period overshoots, up to twice the static displacement; one applied
slower is followed quasi-statically. The acoustic pulse (A u^2 ~ 1 kT at u = 0.2) is a small
correction to a 16 kT bulk pressure rise, so an impedance formula for the pulse was never the
controlling physics.

The model, with the piston ramp driving the gas spring against the load spring:

    M s'' = k_gas (x_p(t) - s) - k s ,    x_p(t) = min(u t, dx) ,    s(0) = s'(0) = 0

k_gas = 0.02457 is MEASURED (Kolafa-Rottner at the run's own eta), k and M are the run's flags, and
there is nothing else to set -- no damping, no fitted amplitude. The apparatus is a PRE-LOADED
spring, so the stored energy is first order in the wall displacement:

    dE(s) = F s + k s^2 / 2 ,    F = N kT Z / L = 1.5749 kT/sigma ,

and the prediction is E_max = max_t dE(s(t)). The pre-load cancels out of the DYNAMICS -- at s = 0
the spring's pull exactly balances the standing gas force, which is why the equation above has no
constant term -- but it does not cancel out of the ENERGY, and that is the whole point.

Pass is claimed only at u = 0.2, where the gas is quasi-static enough for a single gas coordinate
to stand in for the field (tau_push = 40 against L/c_s = 45, marginal). At u = 0.5 and 1.0 the push
is short compared with the sound traversal, the gas cannot be one spring, and the deviation from
the model is reported AS the acoustic contribution -- measured, not assumed -- with the impedance
estimate printed beside it for comparison and no pass claimed.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

ET = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/"
      "experiments_energy_transfer")
RUN = f"{ET}/level3_master_preload_20260921"
N, L, H, R = 100, 78.5, 10.0, 0.5
DX, KSPR = 7.96, 0.5
A_STEP, dA_STEP = 25.91, 4.28
MASSES = (2, 10, 50, 200, 1000)
SPEEDS = (0.05, 0.1, 0.2, 0.5, 1.0)
CELLS = [(m, u) for u in (0.2, 0.5, 1.0) for m in MASSES] + \
        [(m, u) for u in (0.05, 0.1) for m in (50, 200)]

ETA = N * math.pi * R * R / (L * H)
Z = sos.Z_kolafa_rottner_2006(ETA)
DZ = sos.dZ_kolafa_rottner_2006(ETA)
STIFF = Z + ETA * DZ
K_GAS = N * STIFF / L ** 2
CS = math.sqrt(STIFF + Z * Z)
ZG = (N / (L * H)) * CS * H
FGAS = N * Z / L
# ##CHRIS 2026-09-21: THE GAS SPRING IS ADIABATIC, NOT ISOTHERMAL. There is no heat bath in this
# box, so a compression heats the gas and the restoring force rises faster than Z(eta) alone says.
# For a 2D hard-disk gas U = N kT, so dU = -P dV gives dlnT/dlnL = -Z, and
#     k_gas^ad = -dF/dL|_S = (N T / L^2) (Z + eta Z' + Z^2) = (N T / L^2) c_s^2 ,
# exactly Paper 1's measured sound speed. It is 2.01x the isothermal value here, and using the
# isothermal one put the model a factor ~2.3 below the measured wall displacement. Cross-check of
# the adiabat itself: it predicts T_f/T_i = 1.1432 for this 10 % compression against Level 1's
# measured 1.149, agreeing to 0.5 %.
K_GAS_ISO = K_GAS
K_GAS = N * CS ** 2 / L ** 2
S_QS = DX * K_GAS / (KSPR + K_GAS)

# Exact adiabatic gas force F(L), tabulated once and interpolated: the ODE then carries the real
# equation of state instead of a stiffness linearised at the start.
_LG = np.linspace(L - DX - 8.0, L + 8.0, 2001)
def _T_of(Lt):
    lnL = np.log(_LG / L)
    Zs = np.array([sos.Z_kolafa_rottner_2006(N * math.pi * R * R / (x * H)) for x in _LG])
    lnT = np.concatenate([[0.0], np.cumsum(-0.5 * (Zs[1:] + Zs[:-1]) * np.diff(lnL))])
    return np.exp(lnT)
_TG = _T_of(_LG)
_FG = np.array([N * _TG[i] * sos.Z_kolafa_rottner_2006(N * math.pi * R * R / (_LG[i] * H)) / _LG[i]
                for i in range(len(_LG))])
def F_ad(Lx):
    return np.interp(Lx, _LG, _FG)
F0_AD = float(F_ad(L))


def dE(s):
    """Stored energy above the pre-loaded baseline: first order in s, because the spring is loaded."""
    return FGAS * s + 0.5 * KSPR * s * s


def model(ms, u, k=KSPR, dt=1e-3, t_end=None):
    """Integrate M s'' = k_gas (x_p - s) - k s with RK4. Returns (t, s, E_max, t_max)."""
    tau = DX / u
    if t_end is None:
        t_end = tau + 220.0
    n = int(t_end / dt)
    s = v = 0.0
    best, tbest = 0.0, 0.0
    ts, ss = [], []

    def acc(t_, s_):
        xp = min(u * t_, DX)
        # exact adiabatic gas force, not a stiffness linearised at the start
        return (float(F_ad(L - xp + s_)) - F0_AD - k * s_) / ms

    for i in range(n):
        t_ = i * dt
        k1v = acc(t_, s);                 k1s = v
        k2v = acc(t_ + dt/2, s + dt/2*k1s); k2s = v + dt/2*k1v
        k3v = acc(t_ + dt/2, s + dt/2*k2s); k3s = v + dt/2*k2v
        k4v = acc(t_ + dt, s + dt*k3s);     k4s = v + dt*k3v
        s += dt/6 * (k1s + 2*k2s + 2*k3s + k4s)
        v += dt/6 * (k1v + 2*k2v + 2*k3v + k4v)
        e = dE(s)
        if e > best:
            best, tbest = e, t_
        if i % 200 == 0:
            ts.append(t_); ss.append(s)
    return np.array(ts), np.array(ss), best, tbest - tau


def cell(ms, u):
    d = f"{RUN}/k0.5_M{ms}_u{u}"
    tails, W, base_means, ped, blinds = [], [], [], [], []
    tt = None
    for tr in sorted(glob.glob(f"{d}/tr_*.csv")):
        t = pd.read_csv(tr, low_memory=False)
        time = t["Time"].to_numpy(float)
        v = np.abs(t["PistonR_v"].to_numpy(float))
        mv = np.nonzero(v > 1e-12)[0]
        if len(mv) < 2:
            continue
        i0, i1 = mv[0], mv[-1]
        E = t["SpringE"].to_numpy(float)
        # Causally blind window: 0.4 L/c_s, not 0.8: a strong push outruns c_s, and at u = 1 the 0.8 window
        # was NOT blind -- its "bias" came out at 12 kT, i.e. signal. 0.4 needs a 2.5 c_s front to leak.
        # its spread is the pedestal.
        blind = np.nonzero(time - time[i0] <= 0.4 * L / CS)[0]
        blind = blind[blind >= i0]
        if len(blind) < 100:
            continue
        b = float(E[blind].mean())
        base_means.append(b); ped.append(float(E[blind].std(ddof=1)))
        blinds.append(E[blind] - b)
        W.append(float(t["PistonWork"].iloc[i1] - t["PistonWork"].iloc[i0]))
        tail = E[i1:] - b                      # baseline MEAN subtracted, per seed
        ti = time[i1:] - time[i1]
        if tt is None or len(ti) < len(tt):
            tt = ti
        tails.append(tail)
    if not tails:
        return None
    n = min(len(x) for x in tails)
    Ebar = np.mean([x[:n] for x in tails], axis=0)
    Esem = np.std([x[:n] for x in tails], axis=0, ddof=1) / math.sqrt(len(tails))
    j = int(np.argmax(Ebar))
    # The estimator is max-over-time of a noisy ensemble mean, which is biased HIGH: with ~18000
    # time points even pure noise has a positive maximum. Measure that bias with the same estimator
    # on the causally blind window, where the true signal is exactly zero, and subtract it. This is
    # self-calibrating -- no model, no assumption about the noise.
    nb = min(len(x) for x in blinds)
    bias = float(np.max(np.mean([x[:nb] for x in blinds], axis=0)))
    return dict(n=len(tails), t=tt[:n], Ebar=Ebar, Esem=Esem, bias=bias,
                dE=float(Ebar[j]) - bias, dE_raw=float(Ebar[j]),
                ddE=float(Esem[j]), tmax=float(tt[j]),
                W=float(np.mean(W)), dW=float(np.std(W, ddof=1) / math.sqrt(len(W))),
                base=float(np.mean(base_means)), ped=float(np.mean(ped)))


def main():
    print("## Level 3 v3 -- the one-degree-of-freedom model\n")
    print(f"eta = {ETA:.4f}, Z = {Z:.4f}, Z + etaZ' = {STIFF:.4f}, c_s = {CS:.4f}, "
          f"L/c_s = {L/CS:.1f} sigma-time")
    print(f"k_gas isothermal = {K_GAS_ISO:.5f};  ADIABATIC (N/L^2) c_s^2 = {K_GAS:.5f}  "
          f"(ratio {K_GAS/K_GAS_ISO:.2f})")
    print(f"k = {KSPR}, F = N kT Z/L = {FGAS:.4f} kT/sigma; the ODE uses the exact adiabatic F(L)")
    print(f"s_qs = dx k_gas/(k+k_gas) = {S_QS:.4f} sigma  ->  quasi-static dE = {dE(S_QS):.4f} kT")
    print(f"pre-load stored before anything happens: F^2/2k = {FGAS**2/(2*KSPR):.3f} kT; "
          f"equilibrium mean with thermal = {FGAS**2/(2*KSPR) + 0.5*KSPR/(KSPR+K_GAS):.3f} kT\n")

    print("### The two limits\n")
    print(f"  k -> 0    pushrod: the wall follows the piston, s -> dx, eps -> 1 trivially -- "
          f"nothing is captured FROM THE GAS, the piston does the work directly.")
    print(f"  k -> inf  rigid wall: s -> 0, eps -> 0.")
    print(f"  in between the quasi-static capture is F s_qs + k s_qs^2/2 = {dE(S_QS):.4f} kT, "
          f"which is FIRST order in s.\n")

    rows = {}
    print("### Measured against the 1-DOF prediction (no free parameter)\n")
    print("| u | τ_push/(L/c_s) | M_s | T_w [σ] | τ/T_w | ΔE raw | bias | ΔE_coh corrected [kT] | "
          "1-DOF predicted | diff | σ | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (ms, u) in CELLS:
        c = cell(ms, u)
        if not c:
            continue
        rows[(ms, u)] = c
        _, _, Epred, _ = model(ms, u)
        tw = 2 * math.pi / math.sqrt((KSPR + K_GAS) / ms)
        tau = DX / u
        diff = c["dE"] - Epred
        nsig = abs(diff) / c["ddE"] if c["ddE"] > 0 else float("nan")
        c["pred"] = Epred
        verdict = ("within 2σ" if nsig <= 2 else f"{nsig:.1f}σ")
        print(f"| {u:g} | {tau/(L/CS):.2f} | {ms} | {tw:.0f} | {tau/tw:.2f} | {c['dE_raw']:.2f} | {c['bias']:.2f} | "
              f"{c['dE']:.3f} ± {c['ddE']:.3f} | {Epred:.3f} | {diff:+.3f} | {nsig:.1f} | {verdict} |")
    print()

    print("### Criterion: the model where the gas really IS quasi-static\n")
    print(f"  tau_push/(L/c_s): u=0.05 -> {(DX/0.05)/(L/CS):.2f}, u=0.1 -> {(DX/0.1)/(L/CS):.2f}, "
          f"u=0.2 -> {(DX/0.2)/(L/CS):.2f}. Only the first two exceed 1, i.e. only there does the")
    print("  sound cross the box during the push and a single gas coordinate stand in for the field.")
    print("  u = 0.2 was originally chosen for this test; the data say it does not qualify.\n")
    q = [((m, u), rows[(m, u)]) for u in (0.05, 0.1) for m in (50, 200) if (m, u) in rows]
    bad = [(m, abs(r["dE"] - r["pred"]) / r["ddE"]) for m, r in q
           if abs(r["dE"] - r["pred"]) / r["ddE"] > 2]
    print(f"  {len(q)} quasi-static cells; within 2σ: {len(q) - len(bad)}/{len(q)}")
    if bad:
        print("  outside 2σ: " + ", ".join(f"M_s={m} ({s:.1f}σ)" for m, s in bad))
    print(f"  VERDICT: {'PASS' if not bad else 'FAIL'}\n")

    print("### The acoustic contribution, measured as the deviation from the model\n")
    print("| u | M_s | τ_push/(L/c_s) | ΔE measured | 1-DOF | deviation [kT] | A u² | T_imp·Au² |")
    print("|---|---|---|---|---|---|---|---|")
    for u in (0.5, 1.0):
        for ms in MASSES:
            if (ms, u) not in rows:
                continue
            r = rows[(ms, u)]
            zw = math.sqrt(ms * (KSPR + K_GAS))
            T = 4 * zw * ZG / (zw + ZG) ** 2
            print(f"| {u:g} | {ms} | {(DX/u)/(L/CS):.2f} | {r['dE']:.3f} ± {r['ddE']:.3f} | "
                  f"{r['pred']:.3f} | {r['dE']-r['pred']:+.3f} | {A_STEP*u*u:.2f} | {T*A_STEP*u*u:.2f} |")
    print("\n  No pass is claimed at these speeds: the push is shorter than the sound traversal, so a")
    print("  single gas coordinate cannot stand in for the field, and the deviation is the thing to")
    print("  explain, not a failure of the model outside its stated domain.\n")

    print("### Baseline, and one bias worth stating\n")
    print("| cell | baseline mean [kT] | equilibrium ⟨E⟩ | pedestal rms |")
    print("|---|---|---|---|")
    eq = FGAS ** 2 / (2 * KSPR) + 0.5 * KSPR / (KSPR + K_GAS)
    for ms in MASSES:
        if (ms, 0.2) in rows:
            r = rows[(ms, 0.2)]
            print(f"| M_s = {ms}, u = 0.2 | {r['base']:.3f} | {eq:.3f} | {r['ped']:.3f} |")
    print(f"\n  The wall is released from rest, so it carries no thermal energy at t = 0 and the")
    print(f"  blind window (45 σ-time, shorter than T_w for the heavy walls) samples a partly")
    print(f"  unthermalised state. The baselines above sit below the equilibrium {eq:.3f} kT, so")
    print(f"  ΔE_coh is biased HIGH by that difference -- a few tenths of a kT, stated rather than")
    print(f"  corrected, since correcting it would mean assuming the equilibrium the run has not"
          f" reached.")

    # ---- figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, axm) = plt.subplots(1, 2, figsize=(12.2, 4.7))
    for ms, col in zip(MASSES, ("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple")):
        c = rows.get((ms, 0.2))
        if not c:
            continue
        ax.plot(c["t"], c["Ebar"], lw=1.3, color=col, label=f"$M_s$ = {ms}")
        ax.fill_between(c["t"], c["Ebar"] - c["Esem"], c["Ebar"] + c["Esem"],
                        color=col, alpha=0.18, lw=0)
        ts, ss, _, _ = model(ms, 0.2)
        tau = DX / 0.2
        m = ts >= tau
        ax.plot(ts[m] - tau, dE(ss[m]), lw=1.1, ls="--", color=col, alpha=0.85)
    ax.axhline(dE(S_QS), color="0.4", lw=1.0, ls=":")
    ax.text(3, dE(S_QS) * 1.06, "quasi-static  $Fs_{qs}+ks_{qs}^2/2$", fontsize=8, color="0.35")
    ax.set_xlabel("time after the piston stops  [σ-time]")
    ax.set_ylabel(r"$\Delta E_{spring}$ above the baseline  [$k_BT$]")
    ax.set_xlim(0, 200); ax.grid(True, ls=":", alpha=0.6)
    ax.legend(fontsize=8, title="u = 0.2; dashed = 1-DOF model", title_fontsize=8)
    ax.set_title("Measured vs the parameter-free model", fontsize=11)

    for u, col in zip((0.2, 0.5, 1.0), ("tab:blue", "tab:purple", "tab:red")):
        M = [m for m in MASSES if (m, u) in rows]
        if not M:
            continue
        axm.errorbar(M, [rows[(m, u)]["dE"] for m in M],
                     yerr=[rows[(m, u)]["ddE"] for m in M],
                     fmt="o", color=col, ms=5, capsize=3, lw=1.2, label=f"u = {u:g} measured")
        axm.plot(M, [rows[(m, u)]["pred"] for m in M], "--", color=col, lw=1.2, alpha=0.8)
    axm.set_xscale("log"); axm.set_xlabel("wall mass  $M_s$")
    axm.set_ylabel(r"$\Delta E_{coh}$  [$k_BT$]")
    axm.grid(True, ls=":", alpha=0.6); axm.legend(fontsize=8)
    axm.set_title("Dashed = 1-DOF prediction, no free parameter", fontsize=11)
    fig.tight_layout()
    out = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/"
           "paper2_energytransfer/experiments/final/260921_level3_1dof")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    print(f"\nfigure -> {out}.png/.pdf")


if __name__ == "__main__":
    main()
