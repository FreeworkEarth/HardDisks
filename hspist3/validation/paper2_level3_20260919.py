#!/usr/bin/env python3
"""##CHRIS 2026-09-19: Paper 2 Level 3 pilot -- one gas drives a spring-loaded wall (geometry C).

Box: rigid wall | spring (k = 5) | free wall (M_s = 200) | gas (N_s = 50, eta = 0.1013) | piston.
Five speeds, 25 seeds each, dx = 3.93 sigma of gas. The wall period T_w = 2 pi sqrt(M_s/(k + k_gas))
is 39.5 sigma-time, and the pushes run 786 .. 20 sigma-time, so the ladder crosses from
tau_push >> T_w (the wall follows quasi-statically) to tau_push << T_w (the wall is kicked).

Quantities per seed, all from the trace so the spring term is in the same units as the work:
  W_in            piston work at the moment the piston stops
  E_spring,max    peak spring energy after release, measured from its pre-push level
  t_max           when that peak occurs, relative to piston stop
  partition       gas KE, wall KE and spring energy at the end of the record
Analysis only; the runs are in experiments_energy_transfer/level3_spring_20260919.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

ET = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_energy_transfer"
RUN = f"{ET}/level3_spring_20260919"
NS, H, R, L0c = 50, 10.0, 0.5, 38.75
MS, KSPR = 200.0, 5.0
SPEEDS = (0.005, 0.02, 0.05, 0.1, 0.2)

ETA = NS * math.pi * R * R / (L0c * H)
Z = sos.Z_kolafa_rottner_2006(ETA); DZ = sos.dZ_kolafa_rottner_2006(ETA)
K_GAS = NS * (Z + ETA * DZ) / L0c ** 2
OMEGA = math.sqrt((KSPR + K_GAS) / MS)
T_W = 2 * math.pi / OMEGA
CS = math.sqrt(Z + ETA * DZ + Z * Z)


def cells(u):
    out = []
    for tr in sorted(glob.glob(f"{RUN}/u{u:g}/tr_*.csv")):
        d = pd.read_csv(tr, low_memory=False)
        t = d["Time"].to_numpy(float)
        v = np.abs(d["PistonR_v"].to_numpy(float))
        mv = np.nonzero(v > 1e-12)[0]
        if len(mv) == 0:
            continue
        i0, i1 = mv[0], mv[-1]                      # push start .. piston stop
        W = float(d["PistonWork"].iloc[i1] - d["PistonWork"].iloc[i0])
        Es = d["SpringE"].to_numpy(float)
        base = float(np.mean(Es[max(0, i0 - 200):i0 + 1]))   # spring level before the push
        after = Es[i1:] - base
        j = int(np.argmax(after))
        vw = d["W0_v"].to_numpy(float)
        out.append(dict(
            W=W, Emax=float(after[j]), tmax=float(t[i1 + j] - t[i1]),
            Egas=float(d["KE_gas_total"].iloc[-1] - d["KE_gas_total"].iloc[i0]),
            Ewall=float(0.5 * MS * vw[-1] ** 2),
            Espr=float(Es[-1] - base),
            tail=after, ttail=t[i1:] - t[i1]))
    return out


def ms(a):
    a = np.asarray(a, float)
    return a.mean(), a.std(ddof=1) / math.sqrt(len(a))


def main():
    print("### Level 3 pilot -- geometry C, one gas driving a spring-loaded wall\n")
    print(f"eta = {ETA:.5f}, Z = {Z:.4f}, c_s = {CS:.4f}, k_gas = {K_GAS:.5f}, k = {KSPR}")
    print(f"omega_w = {OMEGA:.5f}, T_w = {T_W:.2f} sigma-time; 2L/c_s = {2*L0c/CS:.2f}")
    print(f"k/(k + k_gas) = {KSPR/(KSPR+K_GAS):.4f}  (the share of the total stiffness the spring holds)\n")

    print("| u | tau_push [σ] | tau_push/T_w | seeds | ⟨W_in⟩ [kT] | ⟨E_spring,max⟩ [kT] | t_max after stop [σ] | ε = ⟨E⟩/⟨W⟩ | ⟨E/W⟩ |")
    print("|---|---|---|---|---|---|---|---|---|")
    rows = {}
    for u in SPEEDS:
        c = cells(u)
        if not c:
            continue
        W, dW = ms([x["W"] for x in c]); E, dE = ms([x["Emax"] for x in c])
        tm, dtm = ms([x["tmax"] for x in c])
        r, dr = ms([x["Emax"] / x["W"] for x in c])
        tau = 3.93 / u
        rows[u] = (c, W, dW, E, dE)
        print(f"| {u:g} | {tau:.0f} | {tau/T_W:.1f} | {len(c)} | {W:.3f} ± {dW:.3f} | {E:.4f} ± {dE:.4f} | "
              f"{tm:.1f} ± {dtm:.1f} | {E/W:.4f} ± {dE/W:.4f} | {r:.4f} ± {dr:.4f} |")

    print("\n### Where the energy is at the end of the record\n")
    print("| u | ΔKE_gas | KE_wall | E_spring | sum | ⟨W_in⟩ |")
    print("|---|---|---|---|---|---|")
    for u in SPEEDS:
        if u not in rows: continue
        c, W, _, _, _ = rows[u]
        g, _ = ms([x["Egas"] for x in c]); w, _ = ms([x["Ewall"] for x in c]); s, _ = ms([x["Espr"] for x in c])
        print(f"| {u:g} | {g:.3f} | {w:.4f} | {s:.4f} | {g+w+s:.3f} | {W:.3f} |")

    print("\n### The two limits\n")
    # quasi-static: the wall sits where the compressed gas balances the spring
    for u in SPEEDS:
        if u not in rows: continue
        c, W, _, E, dE = rows[u]
        L_f = L0c - 3.93
        eta_f = NS * math.pi * R * R / (L_f * H)
        Zf = sos.Z_kolafa_rottner_2006(eta_f)
        dP = (NS * Zf / (L_f * H)) - (NS * Z / (L0c * H))     # pressure rise, kT per sigma^2
        x_eq = dP * H / (KSPR + K_GAS)
        E_qs = 0.5 * KSPR * x_eq ** 2
        tau = 3.93 / u
        lab = "quasi-static" if tau > 3 * T_W else ("impulsive" if tau < T_W / 2 else "crossover")
        print(f"  u = {u:<6g} tau/T_w = {tau/T_W:5.1f}  {lab:12s} "
              f"E_qs = {E_qs:.4f}   measured {E:.4f} ± {dE:.4f}   ratio {E/E_qs:.2f}")

    # figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, axe) = plt.subplots(1, 2, figsize=(11.6, 4.6))
    for u in SPEEDS:
        if u not in rows: continue
        c = rows[u][0]
        n = min(len(x["tail"]) for x in c)
        tt = c[0]["ttail"][:n]
        mean = np.mean([x["tail"][:n] for x in c], axis=0)
        ax.plot(tt, mean, lw=1.2, label=f"u = {u:g}")
    ax.axvline(T_W, color="0.5", ls=":", lw=1.0)
    ax.text(T_W * 1.05, ax.get_ylim()[1] * 0.92, "$T_w$", fontsize=9, color="0.4")
    ax.set_xlabel("time after the piston stops  [σ-time]")
    ax.set_ylabel("spring energy above its pre-push level  [$k_BT$]")
    ax.set_xlim(0, 200); ax.grid(True, ls=":", alpha=0.6)
    ax.legend(fontsize=8.5); ax.set_title("Spring ring-down, seed-averaged", fontsize=11)
    U = np.array([u for u in SPEEDS if u in rows])
    EPS = np.array([rows[u][3] / rows[u][1] for u in U])
    dEPS = np.array([rows[u][4] / rows[u][1] for u in U])
    axe.errorbar(U, EPS, yerr=dEPS, fmt="o-", color="tab:purple", ms=5, capsize=3, lw=1.2)
    axe.set_xscale("log"); axe.set_xlabel("piston speed  u  [σ/τ]")
    axe.set_ylabel(r"$\epsilon = \langle E_{spring,max}\rangle / \langle W_{in}\rangle$")
    axe.grid(True, ls=":", alpha=0.6); axe.set_title("Capture efficiency", fontsize=11)
    fig.tight_layout()
    out = T.plot_path("260919_level3_spring", write=True).replace("/paper1_speedofsound/", "/paper2_energytransfer/")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    print(f"\nfigure -> {out}.png/.pdf")


if __name__ == "__main__":
    main()
