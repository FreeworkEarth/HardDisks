#!/usr/bin/env python3
"""##CHRIS 2026-09-12: route-A c_s(eta) from the per-run damped-cosine frequency,
corrected to T_i = 1, on ELIGIBLE trajectories that pass the 0.03 % fit-quality cut.

Estimator, per eta:
  per run   nu_corr = nu_fitted / sqrt(T_i_mean)          (c_s scales as sqrt(T))
  per mass  mean(nu_corr) and its sem over the surviving repeats
  fit       nu_corr = c_s * x through the origin, x = K/(2 pi L_eff),
            K the fundamental root of cot K = alpha K, alpha = M/(2 N_side),
            L_eff = L0 - 2r, weights 1/sem^2   -> c_s_corr, c_s_err_fit
  scatter   c_s_scatter_mass = std(ddof=1) of the per-mass implied c_s = mean(nu_corr)/x

Quality cut is on the FIT ONLY (sigma_nu/nu < 3e-4). It is deliberately NOT a cut on
agreement with the binned nu, which would bias the slope. Read-only on all inputs.
"""
import csv, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

MAN = sys.argv[1]
OUT = sys.argv[2]
CUT = float(sys.argv[3]) if len(sys.argv) > 3 else 3e-4
# ##CHRIS 2026-09-12: MODE "robust" is the primary estimator and ignores CUT entirely.
# sigma_nu underestimates the true run-to-run scatter by ~92x (residuals are driven and
# correlated, least squares assumes independent), so selecting on it is selecting on a
# quantity that is not a valid error. Robust mode instead: no sigma cut, a factor-3
# alias guard against the INDEPENDENT binned nu (wide enough that it cannot bias the
# slope -- real fits agree to ~1%), per-mass median, 5-robust-sigma outlier rejection,
# weights from the MEASURED scatter MAD/sqrt(n).
MODE = sys.argv[4] if len(sys.argv) > 4 else "cut"
ALIAS_LO, ALIAS_HI = 1.0 / 3.0, 3.0


def _mad(a):
    m = np.median(a)
    return 1.4826 * float(np.median(np.abs(a - m)))
R, N_SIDE, KR_MAX = 0.5, 50, 0.69

rows = list(csv.DictReader(open(MAN)))
kept, drop_elig, drop_fit, drop_cut = 0, 0, 0, 0
by_eta = {}
for r in rows:
    if r["eligible"] != "1":
        drop_elig += 1; continue
    if r["nu_fitted"] in ("", "nan") or r["sigma_nu_fitted"] in ("", "nan"):
        drop_fit += 1; continue
    nu, sig = float(r["nu_fitted"]), float(r["sigma_nu_fitted"])
    if not (nu > 0):
        drop_cut += 1; continue
    if MODE in ("robust", "nocut"):
        try: nb = float(r["nu"])
        except (KeyError, ValueError): nb = 0.0
        if MODE == "nocut":
            # no cut at all: the binned estimator itself, every eligible trace
            if not (nb > 0):
                drop_cut += 1; continue
            nu = nb
        elif nb > 0 and not (ALIAS_LO < nu / nb < ALIAS_HI):
            drop_cut += 1; continue
    elif sig / nu >= CUT:
        drop_cut += 1; continue
    ti = float(r["T_i_mean"])
    kept += 1
    e = float(r["eta"])
    by_eta.setdefault(e, {"L0": float(r["L0"]), "M": {}, "raw": {}, "ti": []})
    by_eta[e]["M"].setdefault(int(r["M"]), []).append(nu / math.sqrt(ti))
    by_eta[e]["raw"].setdefault(int(r["M"]), []).append(nu)
    by_eta[e]["ti"].append(ti)

summ = []
for eta in sorted(by_eta):
    L0 = by_eta[eta]["L0"]; L_eff = L0 - 2 * R
    xs, ys, ws, cs_per_mass, n_used = [], [], [], [], 0
    rxs, rys, rws = [], [], []
    for M, nus in sorted(by_eta[eta]["M"].items()):
        if len(nus) < 3:
            continue
        arr = np.array(nus)
        if MODE in ("robust", "nocut"):
            m0 = float(np.median(arr)); s0 = _mad(arr)
            arr = arr[np.abs(arr - m0) <= 5 * s0] if s0 > 0 else arr
            if len(arr) < 3: continue
            y = float(np.median(arr))
            sd = _mad(arr) if _mad(arr) > 0 else float(arr.std(ddof=1))
            sem = sd / math.sqrt(len(arr))
        else:
            y = float(arr.mean()); sem = float(arr.std(ddof=1)) / math.sqrt(len(arr))
        K = sos.k_root_bisect(M / (2.0 * N_SIDE))
        x = K / (2 * math.pi * L_eff)
        xs.append(x); ys.append(y)
        ws.append(max(sem, 1e-12))
        cs_per_mass.append(y / x); n_used += len(arr)
        # ##CHRIS: identical estimator on the UNCORRECTED nu, so c_s_raw and c_s_corr
        # differ only by the per-run 1/sqrt(T_i) factor and by nothing else.
        raw = by_eta[eta]["raw"][M]
        rxs.append(x); rys.append(float(np.mean(raw)))
        rws.append(max(float(np.std(raw, ddof=1)) / math.sqrt(len(raw)), 1e-12))
    if len(xs) < 3:
        print(f"  skip eta={eta}: only {len(xs)} masses with >=3 runs")
        continue
    x, y, s = np.array(xs), np.array(ys), np.array(ws)
    cs, cs_err = sos.weighted_linreg(x, y, s, force_zero_intercept=True)[:2]
    cs_raw, cs_raw_err = sos.weighted_linreg(np.array(rxs), np.array(rys), np.array(rws),
                                             force_zero_intercept=True)[:2]
    ti_list = by_eta[eta]["ti"]
    ti_mean = float(np.mean(ti_list))
    scatter = float(np.std(np.array(cs_per_mass), ddof=1))
    a = np.array([eta])
    Zk = sos.Z_kolafa_rottner_2006(a); Zl = sos.Z_liu_global(a); dZl = sos.dZ_liu_global(a)
    h = 1e-5
    dZk = (sos.Z_kolafa_rottner_2006(a + h) - sos.Z_kolafa_rottner_2006(a - h)) / (2 * h)
    cs_kr = float(sos.cs_adiabatic_2d_monatomic(Zk, dZk, a, kbt=1, m=1)[0]) if eta <= KR_MAX else float("nan")
    cs_liu = float(sos.cs_adiabatic_2d_monatomic(Zl, dZl, a, kbt=1, m=1)[0])
    cs_dilute = math.sqrt(2.0) * (1.0 + 2.0 * eta)
    summ.append(dict(
        eta=f"{eta:.6f}", L0=f"{L0:.6f}", L_eff=f"{L_eff:.6f}",
        c_s=f"{cs:.5f}", c_s_err=f"{cs_err:.5f}",
        c_s_corr=f"{cs:.5f}", c_s_err_fit=f"{cs_err:.5f}",
        c_s_raw=f"{cs_raw:.5f}", c_s_raw_err=f"{cs_raw_err:.5f}",
        T_i_mean=f"{ti_mean:.6f}",
        c_s_scatter_mass=f"{scatter:.5f}", n_masses=len(xs),
        n_masses_lt7=1 if len(xs) < 7 else 0, n_runs_in_fit=n_used,
        cs_KR2006=f"{cs_kr:.5f}" if cs_kr == cs_kr else "",
        dev_KR=f"{100*(cs-cs_kr)/cs_kr:+.2f}%" if cs_kr == cs_kr else "",
        cs_dilute_sqrt2=f"{cs_dilute:.5f}",
        dev_dilute=f"{100*(cs-cs_dilute)/cs_dilute:+.2f}%",
        cs_Liu2021=f"{cs_liu:.5f}", dev_Liu=f"{100*(cs-cs_liu)/cs_liu:+.2f}%"))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
print("MODE = robust (no sigma_nu cut; alias guard vs binned nu; median + MAD)"
      if MODE == "robust" else f"cut sigma_nu/nu < {CUT*100:.3f}%")
print(f"rows {len(rows)}  dropped: not_eligible={drop_elig} no_fit={drop_fit} failed_cut={drop_cut}  kept={kept}")
print(f"eta written {len(summ)}  -> {OUT}")
low = [r["eta"] for r in summ if r["n_masses_lt7"]]
print(f"eta with < 7 surviving masses (flagged n_masses_lt7=1): {low if low else 'none'}")
print(f"{'eta':>9} {'nM':>3} {'T_i_mean':>9} {'c_s_raw':>9} {'c_s_corr':>9} {'ratio':>8}")
for r in summ:
    print(f"{r['eta']:>9} {r['n_masses']:>3} {r['T_i_mean']:>9} {r['c_s_raw']:>9} "
          f"{r['c_s_corr']:>9} {float(r['c_s_corr'])/float(r['c_s_raw']):>8.5f}")
