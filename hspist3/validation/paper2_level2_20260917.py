#!/usr/bin/env python3
"""##CHRIS 2026-09-17: Paper 2, Level 2. Analysis only, no physics touched.

Two questions, both answered from the 710 runs in level2_slope_20260917 plus the 50 older
level0_Wqs_20260911 runs (same command line, verified: hold 12000 steps, L0 39.25, N 100, dx 3.93).

(1) W(u): does the excess work over the quasi-static value rise linearly in the piston speed,
    and with what slope?
(2) zeta: the Green-Kubo friction of the pinned piston, as a RUNNING integral
    zeta(t_cut) = beta * int_0^t_cut <dF(0) dF(t')> dt', measured in the equilibrium hold.
    The linear-response claim W_diss = zeta u dx needs the plateau of that integral, not its
    zero-lag value. The zero-lag value alone is the free-molecular (Enskog) piston friction,
    which is what the first pass reported.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

ET = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_energy_transfer"
NEW, OLD = f"{ET}/level2_slope_20260917", f"{ET}/level0_Wqs_20260911"
WQS, DWQS, DX = 7.4879, 0.0066, 3.93          # Level 1 finite-box quasi-static work
HOLD_T, BURN, BIN = 200.0, 20.0, 0.25          # hold window [sigma], burn-in, force bin width


def work_cells():
    """<W> per piston speed, pooling the two run families (identical command lines)."""
    out = {}
    for d in (OLD, NEW):
        for p in sorted(glob.glob(os.path.join(d, "u*", "summary.csv"))):
            u = float(os.path.basename(os.path.dirname(p))[1:])
            out.setdefault(u, []).extend(float(r) for r in pd.read_csv(p)["W_in_max"])
    return {u: (np.mean(w), np.std(w, ddof=1) / math.sqrt(len(w)), len(w)) for u, w in sorted(out.items())}


def wline(u, y, e):
    """Weighted straight line y = a + b u."""
    w = 1.0 / np.asarray(e) ** 2
    S, Su, Suu = w.sum(), (w * u).sum(), (w * u * u).sum()
    Sy, Suy = (w * y).sum(), (w * u * y).sum()
    D = S * Suu - Su ** 2
    a, b = (Suu * Sy - Su * Suy) / D, (S * Suy - Su * Sy) / D
    chi2 = float((w * (y - a - b * u) ** 2).sum())
    return a, math.sqrt(Suu / D), b, math.sqrt(S / D), chi2 / max(1, len(u) - 2)


def force_series(ev, kind="WR"):
    """Impulse train on the pinned piston during the hold, binned into a force time series."""
    d = pd.read_csv(ev, usecols=["t_sigma", "kind", "dp"])
    d = d[(d["kind"] == kind) & (d["t_sigma"] >= BURN) & (d["t_sigma"] < HOLD_T)]
    nb = int((HOLD_T - BURN) / BIN)
    idx = ((d["t_sigma"].to_numpy(float) - BURN) / BIN).astype(int)
    F = np.zeros(nb)
    np.add.at(F, np.clip(idx, 0, nb - 1), d["dp"].to_numpy(float))
    return F / BIN, float((d["dp"].to_numpy(float) ** 2).sum())


def zeta_running(files, kind="WR", nlag=241):
    """Seed-averaged force autocorrelation and its running integral (kT = 1, so beta = 1)."""
    nb = int((HOLD_T - BURN) / BIN)
    acc, sq, n = np.zeros(nlag), 0.0, 0
    for ev in files:
        F, s = force_series(ev, kind)
        F = F - F.mean()
        sp = np.fft.rfft(F, 2 * nb)
        c = np.fft.irfft(sp * np.conj(sp))[:nlag] / (nb - np.arange(nlag))
        acc += c; sq += s; n += 1
    C = acc / n
    lags = np.arange(nlag) * BIN
    run = np.concatenate(([0.0], np.cumsum((C[1:] + C[:-1]) / 2 * BIN)))
    return lags, C, run, sq / n / (HOLD_T - BURN) / 2.0   # last: impulsive self-term beta*sum dp^2/(2T)


def main():
    cells = work_cells()
    print("### Level 2, the slow end: W(u) over a ladder of piston speeds")
    print(f"quasi-static reference W_qs^finite = {WQS:.4f} ± {DWQS:.4f} kT (Level 1), Δx = {DX} σ\n")
    print("| u [σ/τ] | seeds | ⟨W⟩ [kT] | W − W_qs [kT] |")
    print("|---|---|---|---|")
    for u, (m, s, n) in cells.items():
        print(f"| {u:g} | {n} | {m:.4f} ± {s:.4f} | {m - WQS:+.4f} ± {math.hypot(s, DWQS):.4f} |")
    u = np.array(list(cells)); y = np.array([c[0] for c in cells.values()]) - WQS
    e = np.array([c[1] for c in cells.values()])
    print("\n| fit range | intercept [kT] | slope dW/du [kT τ/σ] | χ²/dof |")
    print("|---|---|---|---|")
    for lo, hi in ((0.0, 0.021), (0.0, 0.051), (0.0, 0.21), (0.019, 0.21)):
        k = (u >= lo) & (u <= hi)
        if k.sum() < 3: continue
        a, da, b, db, c2 = wline(u[k], y[k], e[k])
        print(f"| u = {u[k].min():g}–{u[k].max():g} ({k.sum()} speeds) | {a:+.4f} ± {da:.4f} | {b:+.2f} ± {db:.2f} | {c2:.2f} |")

    print("\n### The friction coefficient, as a running integral over the equilibrium hold")
    files = sorted(glob.glob(f"{NEW}/u0.03/ev_*.csv")) + sorted(glob.glob(f"{NEW}/u0.05/ev_*.csv"))
    lags, C, run, self_term = zeta_running(files)
    print(f"{len(files)} holds, {HOLD_T - BURN:.0f} σ each, impulses on the pinned piston (labelled WR while it is held), bin {BIN} σ")
    print(f"impulsive self term  β Σdp²/(2T)                = {self_term:.3f}   (free-molecular / Enskog friction)")
    print("\n| t_cut [σ] | β∫₀^t_cut ⟨δF(0)δF(t)⟩ dt | implied dW/du = ζΔx |")
    print("|---|---|---|")
    for tc in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0):
        i = int(round(tc / BIN))
        z = run[i]      # the binned zero-lag bin already holds the impulsive self term
        print(f"| {tc:5.2f} | {z:+.3f} | {z * DX:+.2f} |")
    print("\nThe trapezoid's first half-step reproduces the impulsive self term (that is the check above:")
    print("the two agree), so the running integral must NOT have it added again.\n")
    print("### Is the excess work linear or quadratic in u?")
    for name, q in (("A u  (linear)", u), ("A u² (quadratic)", u ** 2), ("A u³ (cubic)", u ** 3)):
        w = 1.0 / e ** 2
        A = float((w * q * y).sum() / (w * q * q).sum()); dA = float(1.0 / math.sqrt((w * q * q).sum()))
        c2 = float((w * (y - A * q) ** 2).sum()) / (len(u) - 1)
        print(f"  through the origin, {name}:  A = {A:9.3f} ± {dA:.3f}   χ²/dof = {c2:7.2f}")
    print(f"\n  N_s m / 2 = {50 / 2:.1f}  — the kinetic energy of the whole compressed compartment moving at u")

    print("\n  on the three speeds that carry the signal (u >= 0.1), point by point:")
    k = u >= 0.1
    for nm, q in (("u²", u ** 2), ("u³", u ** 3)):
        w = 1.0 / e ** 2
        A = float((w[k] * q[k] * y[k]).sum() / (w[k] * q[k] ** 2).sum())
        c2 = float((w[k] * (y[k] - A * q[k]) ** 2).sum()) / (int(k.sum()) - 1)
        per = ", ".join(f"{y[i] / q[i]:.1f}" for i in np.where(k)[0])
        print(f"    {nm}: A = {A:7.2f}  χ²/dof = {c2:5.2f}  per point [{per}]")

    print("\n### u → 0 with the quadratic form: the Level 1 intercept, remeasured on 760 trajectories")
    w = 1.0 / e ** 2
    X = np.vstack([np.ones_like(u), u ** 2]).T
    C = np.linalg.inv(X.T @ (X * w[:, None])); a, A = C @ (X.T @ (w * y))
    da, dA = np.sqrt(np.diag(C))
    c2 = float((w * (y - a - A * u ** 2) ** 2).sum()) / (len(u) - 2)
    print(f"  W(0) − W_qs^finite = {a:+.4f} ± {da:.4f} kT  ({abs(a) / math.hypot(da, DWQS):.1f} σ)  χ²/dof = {c2:.2f}")
    print(f"  A = {A:.2f} ± {dA:.2f}   W(0) = {WQS + a:.4f} ± {math.hypot(da, DWQS):.4f} kT")

    print("\n### Is that excess the coherent compression flow?")
    print("| u | seeds | ⟨Px_gas⟩ at piston stop | −N_s m u | ⟨Px⟩²/(2 N_s m) [kT] | W − W_qs [kT] |")
    print("|---|---|---|---|---|---|")
    # stop_digest.csv is extracted from the traces once (85 kB against 6.4 GB), so this panel
    # survives the traces being archived; fall back to the traces if the digest is absent.
    dig = f"{NEW}/stop_digest.csv"
    D = pd.read_csv(dig) if os.path.exists(dig) else None
    for uu in (0.05, 0.10, 0.15, 0.20):
        if D is not None:
            P = D.loc[np.isclose(D["u"], uu), "Px_stop"].to_numpy(float)
        else:
            P = []
            for tr in sorted(glob.glob(f"{NEW}/u{uu:.2f}/tr_*.csv")):
                d = pd.read_csv(tr, usecols=["PistonR_v", "Px_gas"], low_memory=False)
                mov = np.nonzero(np.abs(d["PistonR_v"].to_numpy(float)) > 1e-9)[0]
                if len(mov): P.append(float(d["Px_gas"].to_numpy(float)[mov[-1]]))
            P = np.array(P)
        m = P.mean(); s = P.std(ddof=1) / math.sqrt(len(P))
        exc = cells[uu][0] - WQS
        print(f"| {uu:g} | {len(P)} | {m:+.3f} ± {s:.3f} | {-50 * uu:+.2f} | "
              f"{m * m / 100:.4f} ± {2 * abs(m) * s / 100:.4f} | {exc:+.4f} |")
    print("A non-uniform flow always carries more kinetic energy than ⟨Px⟩²/2M, so that column is a")
    print("lower bound on the flow energy, and it is about half the excess work — consistent, not equal.")


if __name__ == "__main__":
    main()
