#!/usr/bin/env python3
"""##CHRIS 2026-09-15: A2 (fixed eta, N ladder) -- slow mode and damping against N, and the finite-size form.

Per trace: resonance fit (same Whittle fit as damping_test_20260915) giving Gamma/f_0; sd of the 5-period running
mean of x (the slow mode), absolute and relative to L0; and the Paper 1 primary frequency (largest bin at
f >= nu_pred/2.5) so c_s per (eta, N) comes from the same pass.
A2 geometry has L0 and H both ~ sqrt(N), so the equipartition amplitude L0/sqrt(N [Z+eta Z']) is flat in sigma and
falls as N^-1/2 when measured in units of L0. Acoustic attenuation: Gamma/f_0 ~ 1/L0 ~ N^-1/2.
Finite size: c_s(N) fitted with a + b/sqrt(N), a + b/N, and a + b/sqrt(N) + c/N.

usage: A2_slow_mode_20260915.py ROOT1,ROOT2,... OUT_PREFIX"""
import os, sys, math, csv, glob, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos
from damping_test_20260915 import fit_one

X_EDGE, R = 2.5, 0.5


def x_of(M, N, L0):
    return T.k_root(M / float(N)) / (2 * math.pi * (L0 - 2 * R))


def cell(cdir):
    parts = cdir.rstrip("/").split("/")
    eta = float(parts[-3].split("_")[1].replace("p", ".")); N = int(parts[-2][1:]); M = int(parts[-1][2:])
    out = dict(eta=eta, N=N, M=M, nu=[], g=[], sx=[], L0=None, nfail=0, disc=0)
    for r, p, disc in T.cell_runs(cdir, M):
        if disc:
            out["disc"] += 1; continue
        import pandas as pd
        d = pd.read_csv(p, usecols=["Time", "Displacement(σ)", "Predicted_Frequency", "L0"])
        t = d["Time"].to_numpy(float); x = d["Displacement(σ)"].to_numpy(float)
        nup = float(d["Predicted_Frequency"].iloc[0]); out["L0"] = float(d["L0"].iloc[0])
        dt = (t[-1] - t[0]) / (len(t) - 1)
        P, df = T._spectrum(x, dt)
        ncyc = len(x) * dt * nup; k = max(1, int(round(ncyc / X_EDGE)))
        out["nu"].append((k + int(np.argmax(P[k:]))) * df)
        u0, g = fit_one(P, df, nup)
        if np.isfinite(u0): out["g"].append(g / u0)
        else: out["nfail"] += 1
        W = max(3, int(round(5.0 / (nup * dt))))
        c = np.cumsum(np.insert(x - x.mean(), 0, 0.0))
        out["sx"].append(float(np.std((c[W:] - c[:-W]) / W)))
    return out


def loglog_slope(xv, yv):
    xv, yv = np.log(np.asarray(xv, float)), np.log(np.asarray(yv, float))
    ok = np.isfinite(xv) & np.isfinite(yv)
    if ok.sum() < 3: return float("nan"), float("nan")
    A = np.vstack([xv[ok], np.ones(ok.sum())]).T
    sol, res, *_ = np.linalg.lstsq(A, yv[ok], rcond=None)
    dof = max(1, ok.sum() - 2)
    s2 = float(res[0]) / dof if len(res) else 0.0
    cov = s2 * np.linalg.inv(A.T @ A)
    return float(sol[0]), float(math.sqrt(max(cov[0, 0], 0.0)))


def fit_forms(Ns, cs, sc):
    forms = {}
    for name, cols in (("a + b/√N", [lambda n: 1 / np.sqrt(n)]), ("a + b/N", [lambda n: 1 / n]),
                       ("a + b/√N + c/N", [lambda n: 1 / np.sqrt(n), lambda n: 1 / n])):
        k = len(cols) + 1
        if len(Ns) < k:
            continue
        A = np.vstack([np.ones_like(Ns)] + [f(Ns) for f in cols]).T / sc[:, None]
        y = cs / sc
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = (A @ beta - y)
        chi2 = float((resid ** 2).sum()); dof = len(Ns) - k
        cov = np.linalg.inv(A.T @ A)
        forms[name] = dict(a=float(beta[0]), a_err=float(math.sqrt(cov[0, 0])), chi2=chi2, dof=dof,
                           params=[float(b) for b in beta[1:]])
    return forms


if __name__ == "__main__":
    roots = [r for r in sys.argv[1].split(",") if r]; prefix = sys.argv[2]
    dirs = [c for root in roots for c in sorted(glob.glob(f"{root}/eta_*/N*/m_*"))]
    with Pool(10) as pool:
        cells = pool.map(cell, dirs, chunksize=1)
    merged = {}
    for c in cells:
        k = (c["eta"], c["N"], c["M"])
        m = merged.setdefault(k, dict(c, nu=[], g=[], sx=[], nfail=0, disc=0))
        for f in ("nu", "g", "sx"): m[f] += c[f]
        m["nfail"] += c["nfail"]; m["disc"] += c["disc"]
        if m["L0"] is None: m["L0"] = c["L0"]
    per_mass, pts = [], {}
    for (eta, N) in sorted({(k[0], k[1]) for k in merged}):
        cs = [v for k, v in merged.items() if (k[0], k[1]) == (eta, N) and len(v["nu"]) >= 3]
        if len(cs) < 3: continue
        L0 = cs[0]["L0"]
        xs = [x_of(c["M"], N, c["L0"]) for c in cs]; ys = [float(np.mean(c["nu"])) for c in cs]
        c_s = float(np.sum(np.array(xs) * np.array(ys)) / np.sum(np.array(xs) ** 2))
        scat = float(np.std(np.array(ys) / np.array(xs), ddof=1))
        for c in cs:
            per_mass.append(dict(eta=f"{eta:.2f}", N=N, M=c["M"], n=len(c["nu"]),
                                 gamma_over_f0=float(np.median(c["g"])) if c["g"] else float("nan"),
                                 sd_x_slow=float(np.mean(c["sx"])), sd_x_slow_over_L0=float(np.mean(c["sx"])) / c["L0"]))
        g_all = [x for c in cs for x in c["g"]]; sx_all = [x for c in cs for x in c["sx"]]
        pts[(eta, N)] = dict(L0=L0, c_s=c_s, scatter=scat, n_masses=len(cs), n_runs=sum(len(c["nu"]) for c in cs),
                             disc=sum(c["disc"] for c in cs), nfail=sum(c["nfail"] for c in cs),
                             gamma=float(np.median(g_all)), sd_slow=float(np.median(sx_all)),
                             sd_slow_rel=float(np.median(sx_all)) / L0,
                             masses=sorted(c["M"] for c in cs))
    T.write_csv(prefix + "_slow_mode_per_mass.csv", per_mass)
    print("### Slow mode and damping per (η, N)   [x_slow = 5-period running mean of x; Γ/f_0 from the per-trajectory resonance fit]")
    print("| η | N | L0 | masses | trajectories | fits failed | median Γ/f_0 | median sd(x_slow) [σ] | sd(x_slow)/L0 | c_s (X = 2.5) | ± mass scatter |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for (eta, N), p in sorted(pts.items()):
        print(f"| {eta:.2f} | {N} | {p['L0']:.4f} | {','.join(str(m) for m in p['masses'])} | {p['n_runs']} | {p['nfail']} | "
              f"{p['gamma']:.4f} | {p['sd_slow']:.4f} | {p['sd_slow_rel']:.6f} | {p['c_s']:.4f} | {p['scatter']:.4f} |")
    print("\n### Exponents from log-log fits against N")
    print("| η | N values | sd(x_slow) [σ] exponent | sd(x_slow)/L0 exponent | Γ/f_0 exponent |")
    print("|---|---|---|---|---|")
    exps = {}
    for eta in sorted({e for e, _ in pts}):
        rs = sorted((N, v) for (e, N), v in pts.items() if e == eta)
        Ns = np.array([N for N, _ in rs], float)
        a1, e1 = loglog_slope(Ns, [v["sd_slow"] for _, v in rs])
        a2, e2 = loglog_slope(Ns, [v["sd_slow_rel"] for _, v in rs])
        a3, e3 = loglog_slope(Ns, [v["gamma"] for _, v in rs])
        exps[eta] = (rs, a1, e1, a2, e2, a3, e3)
        print(f"| {eta:.2f} | {' '.join(str(int(n)) for n in Ns)} | {a1:+.3f} ± {e1:.3f} | {a2:+.3f} ± {e2:.3f} | {a3:+.3f} ± {e3:.3f} |")
    print("\n### Finite-size form (c_s against N), weighted by the mass scatter")
    print("| η | N values | form | a (c_∞) | ± | χ² | dof | KR | (a − KR)/KR [%] |")
    print("|---|---|---|---|---|---|---|---|---|")
    form_rows = []
    for eta in sorted({e for e, _ in pts}):
        rs = sorted((N, v) for (e, N), v in pts.items() if e == eta)
        Ns = np.array([N for N, _ in rs], float); c = np.array([v["c_s"] for _, v in rs])
        sc = np.array([max(v["scatter"], 1e-12) for _, v in rs]); kr = T.kr_cs(eta)
        for name, f in fit_forms(Ns, c, sc).items():
            print(f"| {eta:.2f} | {' '.join(str(int(n)) for n in Ns)} | {name} | {f['a']:.4f} | {f['a_err']:.4f} | {f['chi2']:.3f} | {f['dof']} | {kr:.4f} | {100*(f['a']-kr)/kr:+.2f} |")
            form_rows.append(dict(eta=f"{eta:.2f}", form=name, a=f"{f['a']:.5f}", a_err=f"{f['a_err']:.5f}",
                                  chi2=f"{f['chi2']:.4f}", dof=f["dof"], KR=f"{kr:.5f}", dev_pct=f"{100*(f['a']-kr)/kr:+.3f}"))
    T.write_csv(prefix + "_finite_size_forms.csv", form_rows)
    for eta in sorted({e for e, _ in pts}):
        fr = [r for r in form_rows if r["eta"] == f"{eta:.2f}"]
        if len(fr) > 1:
            a = np.array([float(r["a"]) for r in fr]); kr = T.kr_cs(eta)
            print(f"η = {eta:.2f}: c_∞ across forms {a.min():.4f} … {a.max():.4f} (spread {a.max()-a.min():.4f} = "
                  f"{100*(a.max()-a.min())/np.mean(a):.2f} %), deviation from KR {100*(a.min()-kr)/kr:+.2f} … {100*(a.max()-kr)/kr:+.2f} %")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    etas = sorted(exps); fig, axs = plt.subplots(2, len(etas), figsize=(3.5 * len(etas), 7), squeeze=False)
    for j, eta in enumerate(etas):
        rs, a1, e1, a2, e2, a3, e3 = exps[eta]
        Ns = np.array([N for N, _ in rs], float)
        ax = axs[0][j]
        y = np.array([v["sd_slow_rel"] for _, v in rs])
        ax.loglog(Ns, y, "o-", color="#2a78d6", label=f"sd(x_slow)/L₀, slope {a2:+.2f} ± {e2:.2f}")
        ax.loglog(Ns, y[0] * (Ns / Ns[0]) ** -0.5, "--", color="#52514e", lw=1, label="N^(−1/2) reference")
        ax.set_title(f"η = {eta:.2f}"); ax.set_xlabel("N"); ax.set_ylabel("sd(x_slow)/L₀")
        ax.grid(True, which="both", ls=":", alpha=0.6); ax.legend(fontsize=7)
        ax = axs[1][j]
        g = np.array([v["gamma"] for _, v in rs])
        ax.loglog(Ns, g, "s-", color="#eb6834", label=f"Γ/f₀, slope {a3:+.2f} ± {e3:.2f}")
        ax.loglog(Ns, g[0] * (Ns / Ns[0]) ** -0.5, "--", color="#52514e", lw=1, label="N^(−1/2) reference")
        ax.set_xlabel("N"); ax.set_ylabel("Γ/f₀"); ax.grid(True, which="both", ls=":", alpha=0.6); ax.legend(fontsize=7)
    fig.suptitle("A2: slow-mode amplitude and resonance width against system size (5-period running mean; fits on the same traces)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}_slow_mode_vs_N.{ext}", dpi=170)
    print("\nhealth: trajectories used", sum(p["n_runs"] for p in pts.values()), "discarded", sum(p["disc"] for p in pts.values()),
          "| fits failed", sum(p["nfail"] for p in pts.values()))
    print("wrote", prefix + "_slow_mode_vs_N.png/.pdf, _slow_mode_per_mass.csv, _finite_size_forms.csv")
