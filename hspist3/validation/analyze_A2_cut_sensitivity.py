#!/usr/bin/env python3
"""##CHRIS 2026-09-12: A2 (fixed-eta N ladder) -- per-mass c_s table and the
sensitivity of c_s(N -> inf) to the frequency estimator and its quality cut.

A2 geometry: N particles, N/2 per side, so N_side = N/2 and alpha = M/(2 N_side) = M/N.
L0 and H both scale as sqrt(N), so the aspect ratio is fixed at each eta.

Variants, all on the same cached fits, all under the strict health contract:
  cut0p03   damped-cosine nu, sigma_nu/nu < 3e-4  (the baseline used so far)
  cut0p10   damped-cosine nu, sigma_nu/nu < 1e-3
  nocut     binned power-spectrum peak, NO cut of any kind
The damped variants also drop frequency aliases (nu/nu_predicted outside [1/3, 3]);
the nocut variant deliberately drops nothing, which is the point of the comparison.

Outputs
  <prefix>_cs_per_mass.csv        item 1(a): c_s and its fit error per (eta, N, M)
                                  plus the weighted slope of c_s against log10 M
  <prefix>_cut_sensitivity.csv    item 1(b): c_s(N -> inf) per variant and the shifts
  <prefix>_cs_vs_N_<variant>.pdf/png
"""
import csv, math, os, sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

CACHE = sys.argv[1]
PREFIX = sys.argv[2]
R = 0.5
# ##CHRIS 2026-09-12: "robust" is the PRIMARY estimator. sigma_nu from the fit covariance
# underestimates the true run-to-run scatter by a median factor of ~92 (measured over 72
# cells), because the residual is driven, correlated motion and least-squares assumes
# independent residuals. Selecting on sigma_nu is therefore selecting on a quantity that
# is not a valid error. The robust variant never consults sigma_nu: it takes the median
# of nu per mass over repeats, rejects only points beyond 5 robust sigma of that median,
# and weights by MAD/sqrt(n) -- the error we can actually measure. The three cut variants
# are kept as a robustness table only.
VARIANTS = [("robust", "damped", None), ("cut0p03", "damped", 3e-4),
            ("cut0p10", "damped", 1e-3), ("nocut", "binned", None)]
ROBUST = "robust"
ALIAS_LO, ALIAS_HI = 1.0 / 3.0, 3.0

rows = list(csv.DictReader(open(CACHE)))


def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


def accepted(r, kind, cut):
    if any(int(r[k]) for k in ("forced_advance", "clamp_repair", "overlap_repair", "wall_overdue")):
        return None                                   # strict health contract, always
    if kind == "binned":
        nu = f(r["nu_binned"])
        return nu if nu > 0 else None
    nu, sig, pred = f(r["nu_damped"]), f(r["sigma_nu"]), f(r["nu_predicted"])
    if not (nu > 0):
        return None
    # cut is None for the robust variant: no precision cut at all, alias guard only.
    if cut is not None and (not (sig >= 0) or sig / nu >= cut):
        return None
    if pred > 0 and not (ALIAS_LO < nu / pred < ALIAS_HI):
        return None
    return nu


def x_of(M, N, L0):
    return sos.k_root_bisect(M / float(N)) / (2 * math.pi * (L0 - 2 * R))


def collect(kind, cut):
    """-> {(eta,N): {"L0":.., "M": {M: [nu,...]}}}"""
    d = {}
    for r in rows:
        nu = accepted(r, kind, cut)
        if nu is None: continue
        key = (float(r["eta"]), int(r["N"]))
        e = d.setdefault(key, {"L0": f(r["L0"]), "M": defaultdict(list)})
        e["M"][int(r["M"])].append(nu)
    return d


def mad(a):
    m = np.median(a)
    return 1.4826 * float(np.median(np.abs(a - m)))


def per_point(d, robust=False):
    """-> {(eta,N): dict with c_s, err, scatter, per-mass rows}"""
    out = {}
    for (eta, N), e in sorted(d.items()):
        L0 = e["L0"]
        xs, ys, ws, mass_rows = [], [], [], []
        for M, nus in sorted(e["M"].items()):
            if len(nus) < 3: continue
            a = np.array(nus); x = x_of(M, N, L0)
            if robust:
                m0 = float(np.median(a)); s0 = mad(a)
                a = a[np.abs(a - m0) <= 5 * s0] if s0 > 0 else a
                if len(a) < 3: continue
                mu = float(np.median(a))
                sd = mad(a) if mad(a) > 0 else float(a.std(ddof=1))
                sem = sd / math.sqrt(len(a))
            else:
                mu = float(a.mean()); sem = float(a.std(ddof=1)) / math.sqrt(len(a))
            xs.append(x); ys.append(mu); ws.append(max(sem, 1e-12))
            mass_rows.append(dict(M=M, n=len(a), nu=mu, nu_sem=sem,
                                  c_s=mu / x, c_s_err=sem / x))
        if len(xs) < 3: continue
        cs, cs_err = sos.weighted_linreg(np.array(xs), np.array(ys), np.array(ws),
                                         force_zero_intercept=True)[:2]
        csm = np.array([m["c_s"] for m in mass_rows])
        out[(eta, N)] = dict(L0=L0, c_s=cs, c_s_err=cs_err,
                             scatter=float(csm.std(ddof=1)), n_masses=len(csm),
                             n_runs=sum(m["n"] for m in mass_rows), masses=mass_rows)
    return out


def kr(eta):
    a = np.array([eta]); h = 1e-5
    Z = sos.Z_kolafa_rottner_2006(a)
    dZ = (sos.Z_kolafa_rottner_2006(a + h) - sos.Z_kolafa_rottner_2006(a - h)) / (2 * h)
    return float(sos.cs_adiabatic_2d_monatomic(Z, dZ, a, kbt=1, m=1)[0]) if eta <= 0.69 else float("nan")


def extrapolate(pts):
    """c_s = a + b/sqrt(N) per eta, weighted by the mass scatter."""
    res = {}
    for eta in sorted({e for e, _ in pts}):
        rs = sorted((N, v) for (e, N), v in pts.items() if e == eta)
        if len(rs) < 3:
            res[eta] = None; continue
        Ns = np.array([N for N, _ in rs], float)
        cs = np.array([v["c_s"] for _, v in rs])
        sc = np.array([max(v["scatter"], 1e-12) for _, v in rs])
        sl, sl_err, ic, ic_err = sos.weighted_linreg(1.0 / np.sqrt(Ns), cs, sc,
                                                     force_zero_intercept=False)[:4]
        res[eta] = dict(slope=sl, slope_err=sl_err, inf=ic, inf_err=ic_err, n_N=len(rs))
    return res


# ---------------- run every variant ----------------
V = {}
for name, kind, cut in VARIANTS:
    pts = per_point(collect(kind, cut), robust=(name == ROBUST))
    V[name] = dict(pts=pts, ext=extrapolate(pts))
    tot = sum(p["n_runs"] for p in pts.values())
    print(f"variant {name:<9} kind={kind:<7} "
          f"points {len(pts):>2}  trajectories in fits {tot:>4}")

BASE = ROBUST

# ---------------- item 1(a): per-mass table + trend in M ----------------
with open(f"{PREFIX}_cs_per_mass.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["variant", "eta", "N", "L0", "M", "n_runs", "nu_mean", "nu_sem",
                "c_s", "c_s_err", "c_s_point", "c_s_point_err",
                "slope_dcs_dlog10M", "slope_err", "t_stat"])
    trend_summary = []
    for name in (BASE, "cut0p03", "cut0p10", "nocut"):
        for (eta, N), p in sorted(V[name]["pts"].items()):
            lm = np.log10([m["M"] for m in p["masses"]])
            cm = np.array([m["c_s"] for m in p["masses"]])
            em = np.array([max(m["c_s_err"], 1e-12) for m in p["masses"]])
            if len(lm) >= 3:
                sl, se = sos.weighted_linreg(lm, cm, em, force_zero_intercept=False)[:2]
            else:
                sl = se = float("nan")
            t = sl / se if se and se == se and se > 0 else float("nan")
            if name == BASE:
                trend_summary.append((eta, N, sl, se, t, len(lm)))
            for m in p["masses"]:
                w.writerow([name, f"{eta:.2f}", N, f"{p['L0']:.6f}", m["M"], m["n"],
                            f"{m['nu']:.9g}", f"{m['nu_sem']:.4g}",
                            f"{m['c_s']:.5f}", f"{m['c_s_err']:.5f}",
                            f"{p['c_s']:.5f}", f"{p['c_s_err']:.5f}",
                            f"{sl:.5f}", f"{se:.5f}", f"{t:.2f}"])

print(f"\n=== 1(a) c_s per divider mass, variant {BASE} ===")
for (eta, N), p in sorted(V[BASE]["pts"].items()):
    print(f"\n eta={eta:.2f}  N={N:<5} L0={p['L0']:.3f}  "
          f"point c_s={p['c_s']:.5f} ± {p['c_s_err']:.5f} (fit), scatter {p['scatter']:.5f}")
    print(f"   {'M':>6} {'runs':>5} {'nu':>13} {'nu sem':>11} {'c_s':>10} {'± fit':>9}")
    for m in p["masses"]:
        print(f"   {m['M']:>6} {m['n']:>5} {m['nu']:>13.8g} {m['nu_sem']:>11.4g} "
              f"{m['c_s']:>10.5f} {m['c_s_err']:>9.5f}")

print(f"\n=== 1(a) does c_s trend with M?  (weighted slope of c_s vs log10 M, {BASE}) ===")
print(f"{'eta':>6} {'N':>6} {'slope':>12} {'± err':>10} {'t':>7}  verdict")
nsig = 0
for eta, N, sl, se, t, nm in trend_summary:
    v = "no trend" if abs(t) < 2 else ("TREND" if abs(t) >= 3 else "marginal")
    if abs(t) >= 3: nsig += 1
    print(f"{eta:>6.2f} {N:>6} {sl:>12.5f} {se:>10.5f} {t:>7.2f}  {v}")
print(f"{nsig} of {len(trend_summary)} points show |t| >= 3")

# ---------------- item 1(b): where does c_s(N->inf) move? ----------------
print("\n=== 1(b) c_s(N -> inf) by variant ===")
hdr = f"{'eta':>6} " + " ".join(f"{n:>22}" for n, _, _ in VARIANTS) + f" {'KR':>9}"
print(hdr)
etas = sorted({e for e, _ in V[BASE]["pts"]})
with open(f"{PREFIX}_cut_sensitivity.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["eta", "KR"] + [c for n, _, _ in VARIANTS
                                for c in (f"{n}_cs_inf", f"{n}_err", f"{n}_slope", f"{n}_n_N")]
               + [c for n in ("cut0p03", "cut0p10", "nocut")
                  for c in (f"shift_{n}_minus_robust", f"shift_{n}_sigma")])
    for eta in etas:
        cells, vals = [], {}
        for n, _, _ in VARIANTS:
            e = V[n]["ext"].get(eta)
            vals[n] = e
            cells.append(f"{e['inf']:>10.5f} ± {e['inf_err']:<9.5f}" if e else f"{'not fitted':>22}")
        print(f"{eta:>6.2f} " + " ".join(cells) + f" {kr(eta):>9.5f}")
        b = vals[BASE]
        out = [f"{eta:.2f}", f"{kr(eta):.5f}"]
        for n, _, _ in VARIANTS:
            e = vals[n]
            out += ([f"{e['inf']:.5f}", f"{e['inf_err']:.5f}", f"{e['slope']:.5f}", e["n_N"]]
                    if e else ["", "", "", ""])
        for n in ("cut0p03", "cut0p10", "nocut"):
            e = vals[n]
            if b and e:
                d = e["inf"] - b["inf"]; s = math.hypot(b["inf_err"], e["inf_err"])
                out += [f"{d:+.5f}", f"{d/s:+.2f}" if s > 0 else ""]
            else:
                out += ["", ""]
        w.writerow(out)

print("\n=== 1(b) shifts vs the 0.03 % baseline (flagged when |shift| > combined error) ===")
print(f"{'eta':>6} {'variant':>9} {'shift':>11} {'comb err':>10} {'ratio':>7}  flag")
for eta in etas:
    b = V[BASE]["ext"].get(eta)
    for n in ("cut0p03", "cut0p10", "nocut"):
        e = V[n]["ext"].get(eta)
        if not (b and e):
            print(f"{eta:>6.2f} {n:>9} {'--':>11} {'--':>10} {'--':>7}  not fitted in one variant"); continue
        d = e["inf"] - b["inf"]; s = math.hypot(b["inf_err"], e["inf_err"])
        print(f"{eta:>6.2f} {n:>9} {d:>+11.5f} {s:>10.5f} {d/s:>+7.2f}  "
              f"{'MOVED' if abs(d) > s else 'within error'}")

# ---------------- figures ----------------
for name, _, _ in VARIANTS:
    pts = V[name]["pts"]
    es = sorted({e for e, _ in pts})
    if not es: continue
    fig, axes = plt.subplots(1, len(es), figsize=(3.5 * len(es), 4.0), squeeze=False)
    for j, eta in enumerate(es):
        ax = axes[0][j]
        rs = sorted((N, v) for (e, N), v in pts.items() if e == eta)
        Ns = np.array([N for N, _ in rs], float)
        cs = np.array([v["c_s"] for _, v in rs])
        sc = np.array([v["scatter"] for _, v in rs])
        ax.errorbar(Ns, cs, yerr=sc, fmt="o", ms=5, capsize=3, color="#1f4e79",
                    label="A2 (bar = 1σ over masses)")
        k = kr(eta)
        if k == k: ax.axhline(k, color="red", lw=1.4, label=f"KR 2006 = {k:.4f}")
        ax.set_xscale("log"); ax.set_xticks(Ns); ax.set_xticks([], minor=True)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlim(Ns.min() / 1.6, Ns.max() * 1.6)
        ax.set_xlabel("N"); ax.set_title(f"η = {eta:.2f}", fontsize=10)
        if j == 0: ax.set_ylabel(r"$c_s$  ($\sigma/t$, $k_BT=m=\sigma=1$)")
        ax.grid(alpha=0.3, which="major"); ax.legend(fontsize=6.5, loc="best")
        ax.text(0.02, 0.03, "masses: " + ", ".join(str(v["n_masses"]) for _, v in rs),
                transform=ax.transAxes, fontsize=6, color="0.35",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    sub = {"robust": "damped-cosine ν · NO σ_ν cut · robust per-mass median, "
                     "weighted by measured scatter (MAD/√n)",
           "cut0p03": "damped-cosine ν, σ_ν/ν < 0.03 %",
           "cut0p10": "damped-cosine ν, σ_ν/ν < 0.10 %",
           "nocut": "binned power-spectrum peak, no quality cut"}[name]
    fig.suptitle("A2: speed of sound vs system size at fixed η and fixed aspect ratio\n"
                 f"{sub} · strict health contract · bar = 1σ scatter over piston masses",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0.02, 1, 0.90])
    for ext in ("pdf", "png"):
        fig.savefig(f"{PREFIX}_cs_vs_N_{name}.{ext}", dpi=170)
    plt.close(fig)
    print(f"wrote {PREFIX}_cs_vs_N_{name}.pdf / .png")
print(f"wrote {PREFIX}_cs_per_mass.csv and {PREFIX}_cut_sensitivity.csv")
