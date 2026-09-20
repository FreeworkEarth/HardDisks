#!/usr/bin/env python3
"""##CHRIS 2026-09-21: Level 3 in the master box -- does the spring catch the acoustic pulse?

Arm (b) only: M_s is scanned at fixed k = 0.5 kT/sigma^2 because that is the discriminator. The
series formula E_qs = k/2 (dx k_gas/(k+k_gas))^2 contains NO M_s, so it predicts a flat line across
the ladder; the impedance formula T_imp = 4 Z_w Z_g/(Z_w+Z_g)^2, Z_w = sqrt(M_s(k+k_gas)),
Z_g = rho c_s H, predicts an interior maximum at Z_w = Z_g. The u ladder is nested inside: the
acoustic term grows as u^2, the quasi-static one does not.

The primary observable is the ENSEMBLE-MEAN E_spring(t). The per-seed maximum is the wrong
estimator and the 2026-09-19 pilot showed why -- it samples the tail of the thermal pedestal and
reports it as capture. The pedestal is incoherent, so it survives averaging only as a constant,
which is subtracted, and its uncertainty falls as 1/sqrt(seeds).

One subtlety that sets the noise floor. The spring holds the standing gas force F = N kT Z/L, so
the wall sits d = F/k from the anchor and E = k/2 (d + delta)^2 has a term F*delta that is LINEAR
in the thermal displacement. Its rms is F sqrt(kT/(k+k_gas)) = 2.17 kT at k = 0.5, far above the
kT/2 one would guess from equipartition alone. It averages to zero over seeds, which is exactly why
the ensemble mean is the estimator and the per-seed peak is not.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

ET = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/"
      "experiments_energy_transfer")
RUN = f"{ET}/level3_master_20260920"
N, L, H, R = 100, 78.5, 10.0, 0.5
DX, KSPR = 7.96, 0.5
A_STEP, dA_STEP = 25.91, 4.28          # Level 2, stepped piston
MASSES = (2, 10, 50, 200, 1000)
SPEEDS = (0.2, 0.5, 1.0)

ETA = N * math.pi * R * R / (L * H)
Z = sos.Z_kolafa_rottner_2006(ETA)
DZ = sos.dZ_kolafa_rottner_2006(ETA)
STIFF = Z + ETA * DZ
K_GAS = N * STIFF / L ** 2
CS = math.sqrt(STIFF + Z * Z)
RHO = N / (L * H)
ZG = RHO * CS * H
FGAS = N * Z / L


def e_qs(k=KSPR):
    return 0.5 * k * (DX * K_GAS / (k + K_GAS)) ** 2


def t_imp(ms, k=KSPR):
    zw = math.sqrt(ms * (k + K_GAS))
    return 4 * zw * ZG / (zw + ZG) ** 2, zw


def cell(ms, u):
    """Seed-averaged spring energy after the stop, plus the per-seed work and ledger."""
    d = f"{RUN}/k0.5_M{ms}_u{u}"
    tails, W, ped, ledger = [], [], [], []
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
        # These traces begin AS the piston starts -- the hold is not logged -- so there is no
        # pre-push window to average. There is a better one anyway: the wall sits L = 78.5 sigma
        # from the piston, so NOTHING the piston does can reach it for L/c_s = 45 sigma-time.
        # The first 0.8 L/c_s of the push is therefore causally pre-signal at the wall, and is the
        # baseline and the pedestal estimate both.
        t_blind = 0.8 * L / CS
        blind = np.nonzero(time - time[i0] <= t_blind)[0]
        blind = blind[blind >= i0]
        if len(blind) < 100:
            continue
        base_slice = E[blind]
        base = float(base_slice.mean())
        ped.append(float(base_slice.std(ddof=1)))       # the pedestal's own spread, per seed
        W.append(float(t["PistonWork"].iloc[i1] - t["PistonWork"].iloc[i0]))
        tail = E[i1:] - base
        ti = time[i1:] - time[i1]
        if tt is None or len(ti) < len(tt):
            tt = ti
        tails.append(tail)
        # Ledger: W_in = dKE_gas + dKE_wall + dE_spring, every term referenced to the SAME index.
        # That index is NOT 0: the spring is armed a step after the first trace row, so SpringE[0]
        # is a spurious exact zero, and using it inflated the residual to 2e-1 -- an analysis bug,
        # not a physics one. Reference to the first armed sample instead.
        armed = np.nonzero(E != 0.0)[0]
        ir = int(armed[0]) if len(armed) else i0
        g = float(t["KE_gas_total"].iloc[-1] - t["KE_gas_total"].iloc[ir])
        vw = float(t["W0_v"].iloc[-1])
        w_ke = 0.5 * ms * vw * vw - 0.5 * ms * float(t["W0_v"].iloc[ir]) ** 2
        spr = float(E[-1] - E[ir])
        tot = float(t["PistonWork"].iloc[-1] - t["PistonWork"].iloc[ir])
        ledger.append(abs(tot - (g + w_ke + spr)) / max(abs(tot), 1e-30))
    if not tails:
        return None
    n = min(len(x) for x in tails)
    Ebar = np.mean([x[:n] for x in tails], axis=0)
    Esem = np.std([x[:n] for x in tails], axis=0, ddof=1) / math.sqrt(len(tails))
    j = int(np.argmax(Ebar))
    return dict(n=len(tails), t=tt[:n], Ebar=Ebar, Esem=Esem,
                Ecoh=float(Ebar[j]), dEcoh=float(Esem[j]), tcoh=float(tt[j]),
                W=float(np.mean(W)), dW=float(np.std(W, ddof=1) / math.sqrt(len(W))),
                ped=float(np.mean(ped)), ledger=float(np.max(ledger)))


def main():
    print("## Level 3 -- the spring against the acoustic pulse (master box, geometry C)\n")
    print(f"eta = {ETA:.4f}, Z = {Z:.4f}, Z + etaZ' = {STIFF:.4f}, c_s = {CS:.4f}")
    print(f"k_gas = {K_GAS:.5f} kT/sigma^2, k = {KSPR} (static offset F/k = {FGAS/KSPR:.2f} sigma)")
    print(f"Z_g = rho c_s H = {ZG:.4f}; impedance match at M_s = {ZG**2/(KSPR+K_GAS):.2f}")
    print(f"E_qs (no M_s dependence, the same in every row below) = {e_qs():.4f} kT")
    print(f"pedestal, predicted linear term F sqrt(kT/(k+k_gas)) = "
          f"{FGAS*math.sqrt(1/(KSPR+K_GAS)):.3f} kT\n")

    rows = {}
    for u in SPEEDS:
        print(f"### u = {u:g}   (A u^2 = {A_STEP*u*u:.2f} kT launched beyond W_qs)\n")
        print("| M_s | Z_w | T_imp | pedestal rms [kT] | ⟨W_in⟩ [kT] | E_coh [kT] | t_coh [σ] | "
              "E_coh/σ | ε_coh | E_ac = T·Au² | E_qs |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        for ms in MASSES:
            c = cell(ms, u)
            if not c:
                print(f"| {ms} | — | — | — | — | (no data) | | | | | |")
                continue
            rows[(ms, u)] = c
            T, zw = t_imp(ms)
            sig = c["Ecoh"] / c["dEcoh"] if c["dEcoh"] > 0 else float("nan")
            print(f"| {ms} | {zw:.2f} | {T:.3f} | {c['ped']:.3f} | {c['W']:.3f} ± {c['dW']:.3f} | "
                  f"{c['Ecoh']:.3f} ± {c['dEcoh']:.3f} | {c['tcoh']:.1f} | {sig:.1f} | "
                  f"{c['Ecoh']/c['W']:.4f} | {T*A_STEP*u*u:.2f} | {e_qs():.4f} |")
        print()

    print("### Pass criteria\n")
    want = [(m, u) for u in SPEEDS for m in MASSES if u >= 0.2]
    have = [k for k in want if k in rows]
    if len(have) < len(want):
        print(f"(ii) CANNOT BE JUDGED: {len(have)}/{len(want)} cells have data. "
              f"An all() over an empty set is not a pass.")
    else:
        bad = [(m, u) for (m, u) in have if rows[(m, u)]["Ecoh"] <= 3 * rows[(m, u)]["dEcoh"]]
        print(f"(ii) E_coh > 3σ in every cell with u >= 0.2: "
              f"{'PASS' if not bad else 'FAIL in ' + str(bad)}  ({len(have)} cells)")
    for u in SPEEDS:
        got = [(m, rows[(m, u)]["Ecoh"] / rows[(m, u)]["W"]) for m in MASSES if (m, u) in rows]
        if not got:
            continue
        best = max(got, key=lambda p: p[1])[0]
        pred = max(MASSES, key=lambda m: t_imp(m)[0])
        print(f"(iii) u = {u:g}: eps_coh peaks at M_s = {best}, impedance predicts {pred}"
              f"  {'-> agrees' if best == pred else '-> DISAGREES'}")
    led = max((r["ledger"] for r in rows.values()), default=float("nan"))
    print(f"(i) ledger with the spring term, worst over all cells: {led:.3e}"
          f"  ({'PASS' if led <= 1e-8 else 'above 1e-8'})")

    # ---- figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, axe) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    for ms in MASSES:
        c = rows.get((ms, 0.5))
        if not c:
            continue
        ax.plot(c["t"], c["Ebar"], lw=1.3, label=f"$M_s$ = {ms}")
        ax.fill_between(c["t"], c["Ebar"] - c["Esem"], c["Ebar"] + c["Esem"], alpha=0.18, lw=0)
    ax.axhline(0, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("time after the piston stops  [σ-time]")
    ax.set_ylabel("seed-mean spring energy above the hold  [$k_BT$]")
    ax.set_xlim(0, 200); ax.grid(True, ls=":", alpha=0.6)
    ax.legend(fontsize=8.5, title="u = 0.5", title_fontsize=8.5)
    ax.set_title("Ensemble-mean ring-up, band = s.e.m.", fontsize=11)

    for u, col in zip(SPEEDS, ("tab:blue", "tab:purple", "tab:red")):
        M = [m for m in MASSES if (m, u) in rows]
        if not M:
            continue
        e = [rows[(m, u)]["Ecoh"] / rows[(m, u)]["W"] for m in M]
        de = [rows[(m, u)]["dEcoh"] / rows[(m, u)]["W"] for m in M]
        axe.errorbar(M, e, yerr=de, fmt="o-", color=col, ms=5, capsize=3, lw=1.2, label=f"u = {u:g}")
    axe.axvline(ZG ** 2 / (KSPR + K_GAS), color="0.4", ls="--", lw=1.1)
    axe.text(ZG ** 2 / (KSPR + K_GAS) * 1.15, axe.get_ylim()[1] * 0.9,
             "impedance\nmatch", fontsize=8, color="0.35")
    axe.set_xscale("log"); axe.set_xlabel("wall mass  $M_s$")
    axe.set_ylabel(r"$\epsilon_{coh} = E_{coh}/\langle W_{in}\rangle$")
    axe.grid(True, ls=":", alpha=0.6); axe.legend(fontsize=8.5)
    axe.set_title("Capture vs wall mass at fixed $k$", fontsize=11)
    fig.tight_layout()
    out = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/"
           "paper2_energytransfer/experiments/final/260921_level3_impedance")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    print(f"\nfigure -> {out}.png/.pdf")


if __name__ == "__main__":
    main()
