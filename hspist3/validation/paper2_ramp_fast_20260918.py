#!/usr/bin/env python3
"""##CHRIS 2026-09-18: Level 2 items 1 and 3 of the 2026-09-18 batch.

Item 1  ramp against step at the same travel: is the u^2 term protocol-launched acoustic energy?
        Plus the decomposition of the excess at piston stop into flow, compression and rest, from
        the env-gated stop snapshots (10 bins across the compressed compartment).
Item 3  fast end, u = 3, 5, 10 at travel 1 sigma, against the single-hit prediction (audit Eq. 7)
        and its variance.

All W_qs values are recomputed on the corrected geometry (paper2_geometry_fix_20260918.py): the
divider is 1 sigma thick and the gas-side face is at 39.75, so the compartment is 38.75 sigma, not
39.25, and the gas is compressed by 3.93 sigma while the piston moves 4.18.
Analysis only.
"""
import glob, math, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos
from paper2_geometry_fix_20260918 import zwall, eta_of, w_qs, ET, FACE, NS, H

RAMP = f"{ET}/level2_ramp_20260918"
FAST = f"{ET}/level2_fast1sigma_20260918"
FAST7 = f"{ET}/level2_fast_20260918"
L_I, DIV_FACE, WALL_R = 38.75, 39.75, 78.50


def path_points():
    pts = [(39.25 - FACE, *zwall(sorted(glob.glob(f"{ET}/level0_Wqs_20260911/u*/ev_*.csv")), 39.25 - FACE)[:2])]
    for nom, tag in ((38.270833, "L38p270833"), (37.291667, "L37p291667"),
                     (36.3125, "L36p3125"), (35.3125, "L35p3125")):
        f = sorted(glob.glob(f"{ET}/level1_Zwall_path_20260916/{tag}/ev_*.csv"))
        if f:
            m, s, _ = zwall(f, nom - FACE)
            pts.append((nom - FACE, m, s))
    return pts


def cell(d):
    w = np.array([float(r) for r in pd.read_csv(f"{d}/summary.csv")["W_in_max"]])
    return w.mean(), w.std(ddof=1) / math.sqrt(len(w)), w.std(ddof=1), len(w)


def snapshots(d, nbins=10):
    """Seed-averaged profiles at the stop instant: coherent flow and compression energy.

    The coherent field is the ENSEMBLE mean of each bin's velocity and density; the per-seed
    scatter about it is thermal and must not be counted as flow. Averaging first is what
    separates them -- squaring first would charge the shot noise to the flow term.
    """
    fs = sorted(glob.glob(f"{d}/snap_*.csv"))
    if not fs:
        return None
    vs, ns, pistons = [], [], []
    for f in fs:
        h = dict(re.findall(r"(\w+)=([-\d.]+)", open(f).readline()))
        piston = float(h["piston_x_sigma"]) - 8.3333333  # screen frame -> box frame
        d_ = pd.read_csv(f, skiprows=1)
        x = d_["x_sigma"].to_numpy(float); vx = d_["vx"].to_numpy(float)
        k = x > DIV_FACE
        x, vx = x[k], vx[k]
        edges = np.linspace(DIV_FACE, piston, nbins + 1)
        idx = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
        n = np.bincount(idx, minlength=nbins).astype(float)
        v = np.array([vx[idx == b].mean() if n[b] else 0.0 for b in range(nbins)])
        vs.append(v); ns.append(n); pistons.append(piston)
    V = np.array(vs); N = np.array(ns); piston = float(np.mean(pistons))
    L = piston - DIV_FACE
    Ab = L * H / nbins
    n_s = len(fs)
    vbar = V.mean(axis=0); nbar = N.mean(axis=0)
    # Unbiased squares. E[<v>^2] = mu^2 + Var(v)/n, so the sampling variance of the seed mean must
    # come off before squaring, otherwise a pure-equilibrium run reports a flow energy of
    # (1/2) nbins kT / n_seeds and a compression energy of c^2 nbins / (2 n_seeds) -- which is
    # exactly the size of the raw numbers at u = 0.05, i.e. all of them.
    v2 = vbar ** 2 - V.var(axis=0, ddof=1) / n_s
    E_flow = 0.5 * float((nbar * v2).sum())
    eta = eta_of(L); Zf = sos.Z_kolafa_rottner_2006(eta)
    dZ = sos.dZ_kolafa_rottner_2006(eta)
    T = 1.1494                      # adiabatic heating over the push (Level 1)
    cs2 = T * (Zf + eta * dZ + Zf * Zf)
    rho = NS / (L * H)
    dn2 = (nbar - NS / nbins) ** 2 - N.var(axis=0, ddof=1) / n_s
    E_comp = float((cs2 * (dn2 / Ab ** 2) * Ab / (2 * rho)).sum())
    return dict(n=n_s, L=L, E_flow=E_flow, E_comp=E_comp, vbar=vbar, nbar=nbar,
                v2=v2, dn2=dn2, Ab=Ab, cs2=cs2, rho=rho, cs=math.sqrt(cs2))


def boot(d, base, nboot=200, nbins=10):
    """Bootstrap the decomposition over seeds: these terms sit close to their own sampling floor."""
    import re as _re
    fs = sorted(glob.glob(f"{d}/snap_*.csv"))
    V, N, P = [], [], []
    for f in fs:
        h = dict(_re.findall(r"(\w+)=([-\d.]+)", open(f).readline()))
        piston = float(h["piston_x_sigma"]) - 8.3333333
        dd = pd.read_csv(f, skiprows=1)
        x = dd["x_sigma"].to_numpy(float); vx = dd["vx"].to_numpy(float)
        k = x > DIV_FACE; x, vx = x[k], vx[k]
        edges = np.linspace(DIV_FACE, piston, nbins + 1)
        idx = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
        n = np.bincount(idx, minlength=nbins).astype(float)
        V.append(np.array([vx[idx == b].mean() if n[b] else 0.0 for b in range(nbins)]))
        N.append(n); P.append(piston)
    V, N = np.array(V), np.array(N)
    L = float(np.mean(P)) - DIV_FACE; Ab = L * H / nbins
    eta = eta_of(L); Zf = sos.Z_kolafa_rottner_2006(eta); dZ = sos.dZ_kolafa_rottner_2006(eta)
    cs2 = 1.1494 * (Zf + eta * dZ + Zf * Zf); rho = NS / (L * H)
    rng = np.random.default_rng(7)
    fl, cm = [], []
    for _ in range(nboot):
        i = rng.integers(0, len(V), len(V))
        vb = V[i].mean(axis=0); nb = N[i].mean(axis=0); ns = len(i)
        v2 = vb ** 2 - V[i].var(axis=0, ddof=1) / ns - base["v2"]
        dn2 = (nb - NS / nbins) ** 2 - N[i].var(axis=0, ddof=1) / ns - base["dn2"]
        fl.append(0.5 * float((nb * v2).sum()))
        cm.append(float((cs2 * (dn2 / Ab ** 2) * Ab / (2 * rho)).sum()))
    return float(np.std(fl, ddof=1)), float(np.std(cm, ddof=1))


def decompose(d, base=None):
    """Flow and compression at stop, with the equilibrium profile at the same L subtracted."""
    s = snapshots(d)
    if s is None:
        return None
    if base is not None:
        # the equilibrium state at L_f is not uniform: the disks layer against the walls, and that
        # standing density structure is not compression energy. Subtract the measured baseline.
        s["E_comp"] = float((s["cs2"] * ((s["dn2"] - base["dn2"]) / s["Ab"] ** 2) * s["Ab"]
                             / (2 * s["rho"])).sum())
        s["E_flow"] = 0.5 * float((s["nbar"] * (s["v2"] - base["v2"])).sum())
    return s


def main():
    pts = path_points()
    WQS, DWQS, _, _ = w_qs(pts, L_I, 74.57 - DIV_FACE)
    print(f"### Level 2 items 1 and 3, on the corrected geometry")
    print(f"W_qs^finite (Δx = 3.93 σ of gas) = {WQS:.4f} ± {DWQS:.4f} kT;  N_s m/6 = {NS / 6:.1f}, N_s m/2 = {NS / 2:.1f}\n")

    print("## Item 1 — ramp against step\n")
    print("| protocol | u | seeds | ⟨W⟩ [kT] | W − W_qs [kT] | (W−W_qs)/u² |")
    print("|---|---|---|---|---|---|")
    res = {}
    for mode in ("step", "ramp"):
        for u in (0.05, 0.10, 0.20):
            d = f"{RAMP}/{'step_' if mode == 'step' else ''}u{u:.2f}"
            m, s, sd, n = cell(d)
            g = m - WQS; sg = math.hypot(s, DWQS)
            res[(mode, u)] = (g, sg)
            print(f"| {mode} | {u:g} | {n} | {m:.4f} ± {s:.4f} | {g:+.4f} ± {sg:.4f} | {g / u ** 2:6.2f} ± {sg / u ** 2:.2f} |")
    print()
    A = {}
    for mode in ("step", "ramp"):
        u = np.array([0.05, 0.10, 0.20])
        y = np.array([res[(mode, uu)][0] for uu in u]); e = np.array([res[(mode, uu)][1] for uu in u])
        w = 1 / e ** 2
        a = (w * u ** 2 * y).sum() / (w * u ** 4).sum(); da = 1 / math.sqrt((w * u ** 4).sum())
        c2 = (w * (y - a * u ** 2) ** 2).sum() / (len(u) - 1)
        A[mode] = (a, da)
        print(f"A_{mode:4s} = {a:6.2f} ± {da:4.2f}  (χ²/dof = {c2:.2f})")
    d_ = A["step"][0] - A["ramp"][0]; sd_ = math.hypot(A["step"][1], A["ramp"][1])
    print(f"difference {d_:.2f} ± {sd_:.2f} ({d_ / sd_:.1f} σ);  A_ramp against N_s m/6 = 8.3: "
          f"{(A['ramp'][0] - NS / 6) / A['ramp'][1]:+.1f} σ;  against N_s m/2 = 25: "
          f"{(A['ramp'][0] - NS / 2) / A['ramp'][1]:+.1f} σ")

    BASE = snapshots(f"{RAMP}/step_u0.005")
    print("\n### Decomposition at piston stop (10 bins, seed-averaged, sampling bias removed,")
    print(f"###  equilibrium profile at the same L subtracted from {BASE['n']} runs at u = 0.005)\n")
    print("| protocol | u | snaps | excess W−W_qs | flow ½mΣN_b v̄_b² | compression | rest | flow/(N_s m u²/6) |")
    print("|---|---|---|---|---|---|---|---|")
    for mode in ("step", "ramp"):
        for u in (0.05, 0.10, 0.20):
            d = f"{RAMP}/{'step_' if mode == 'step' else ''}u{u:.2f}"
            s = decompose(d, BASE)
            if not s: continue
            exc = res[(mode, u)][0]
            rest = exc - s["E_flow"] - s["E_comp"]
            pred = NS * u * u / 6
            ef, ec = boot(d, BASE)
            print(f"| {mode} | {u:g} | {s['n']} | {exc:+.4f} ± {res[(mode, u)][1]:.4f} | {s['E_flow']:.3f} ± {ef:.3f} | "
                  f"{s['E_comp']:+.3f} ± {ec:.3f} | {rest:+.3f} | {s['E_flow'] / pred:.2f} |")

    print("\n## Item 3 — the fast end\n")
    print("Audit Eq. 7: ⟨W⟩ ≃ 2 N_s m (Δx/L) u², and Var W ≃ N_hit[4mu²kT + (2mu²)²], N_hit = N_s Δx/L.")
    for tag, root, dx in (("travel 1.00 σ", FAST, 1.00), ("travel 7.2–7.3 σ", FAST7, None)):
        print(f"\n**{tag}**\n")
        print("| u | seeds | Δx [σ] | ⟨W⟩ [kT] | Eq. 7 | ratio | Var W | Eq. 7 Var | ratio | A = (W−W_qs)/u² |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for u in (3, 5, 10):
            d = f"{root}/u{u}"
            if not os.path.exists(f"{d}/summary.csv"): continue
            m, s, sd, n = cell(d)
            tr = sorted(glob.glob(f"{d}/tr_*.csv") + glob.glob(f"{d}/tr_*.csv.gz"))[0]
            t = pd.read_csv(tr, usecols=["PistonR_x_sigma", "PistonR_v"], low_memory=False)
            v = np.abs(t["PistonR_v"].to_numpy(float)); x = t["PistonR_x_sigma"].to_numpy(float)
            mv = np.nonzero(v > 1e-12)[0]
            travel = abs(x[mv[-1]] - x[mv[0]])
            gas_dx = min(travel, travel - 0.25) if travel > 1.5 else travel   # first 0.25 σ is outside the gas
            gas_dx = travel - 0.25 if travel > 4.0 else travel
            W7 = 2 * NS * (gas_dx / L_I) * u * u
            nhit = NS * gas_dx / L_I
            var7 = nhit * (4 * u * u + (2 * u * u) ** 2)
            wq, _, _, _ = w_qs(pts, L_I, L_I - gas_dx)
            print(f"| {u} | {n} | {gas_dx:.3f} | {m:.2f} ± {s:.2f} | {W7:.2f} | {m / W7:.3f} | "
                  f"{sd ** 2:.1f} | {var7:.1f} | {sd ** 2 / var7:.3f} | {(m - wq) / u ** 2:.2f} |")


if __name__ == "__main__":
    main()
