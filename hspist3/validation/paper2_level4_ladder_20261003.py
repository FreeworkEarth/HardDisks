#!/usr/bin/env python3
"""##CHRIS 2026-10-03: the Level 4 mass-ladder analysis. PRE-REGISTERED -- this file was written and
committed to disk while the runs were still executing and before any ladder record was fitted.

Per mass M in {10, 20, 50, 100, 200}  (M = 10 comes from level4_equilibrium_KOAlength_20261002):

  * omega FIXED from cot K = alpha K, alpha = M/(2 N_s m) = M/100, with Kolafa-Rottner c_s;
  * bias calibration on synthetic OU + AR(2) AT THAT MASS'S OWN RECORD LENGTH, run BEFORE the real
    ACF of that mass is fitted;
  * modelled estimator C(t) = A e^(-t/tau_T) + B e^(-t/tau_r) cos(omega t) adopted if the
    calibration slope >= 0.6, block estimator as cross-check;
  * both observables, T_1 - T_2 and the divider position;
  * reported: tau_T calibrated with error against 50.35 M and 61.84 M; sigma(T_1-T_2) vs 0.199;
    period vs KR / Paper-1 c_s / ideal with the spectral FWHM as the error; tau_r against the
    kinetic-friction and bulk-absorption bounds and Mansour's value at that M (Mhat = M + mN/3).

Ladder: fit tau_T = a M^b over the five masses. Prediction b = 1, a = 50.35 (ideal a = 61.84).
Report b +- sigma and a +- sigma and the sigma distance of b from 1 and of a from both prefactors.

Figures: 261003_p2_ladder_tauT (log-log, both theory lines), 261003_p2_ladder_period (vs alpha,
with cot K = alpha K), 261003_p2_ladder_taur (vs M, against Mansour's form).
"""
import glob, json, math, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
import tests_20260913 as T

REPO = os.path.dirname(os.path.dirname(HERE))
ET = os.path.join(REPO, "hspist3", "experiments_energy_transfer")
BLUE, RED, GREY, ORANGE, GREEN = "#2a78d6", "#e34948", "#52514e", "#eb6834", "#1baf7a"

NS, NTOT, LY, RD = 50, 100, 10.0, 0.5
SPS = 612000 / 10195.000532
LEFF = 38.75 - 2 * RD
GAMMA_ENSKOG = 0.331
BURN = 2000.0
A_HARD, A_IDEAL = 50.35, 61.84

eta = 50 * math.pi * 0.25 / (39.25 * 10)
Zv = float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0])
dZv = float(sos.dZ_kolafa_rottner_2006(np.array([eta]))[0])
CS = math.sqrt(Zv + eta * dZv + Zv * Zv)
KGAS2 = 2.0 * (NS * CS ** 2 / 39.25 ** 2)      # both sides, for the friction bound
GAMMA_FRIC = 4.118

# mass -> (directory, trace-every)
CELLS = {
    10:  (os.path.join(ET, "level4_equilibrium_KOAlength_20261002", "Md10"), 200),
    20:  (os.path.join(ET, "level4_ladder_20261003", "Md20"), 50),
    50:  (os.path.join(ET, "level4_ladder_20261003", "Md50"), 110),
    100: (os.path.join(ET, "level4_ladder_20261003", "Md100"), 210),
    200: (os.path.join(ET, "level4_ladder_20261003", "Md200"), 420),
}


def kroot(alpha):
    return brentq(lambda k: math.cos(k) / math.sin(k) - alpha * k, 1e-9, math.pi - 1e-9)


def acf(x, nl):
    x = x - x.mean(); m = len(x)
    f = np.fft.rfft(x, 2 * m)
    c = np.fft.irfft(f * np.conj(f))[:nl].real
    return c / c[0]


def mansour_tau_r(M, nu):
    Mhat = M + NTOT / 3.0
    dff = GAMMA_ENSKOG * LY * math.sqrt(2.0 / (Mhat * NTOT))
    return 1.0 / (math.pi * nu * dff), dff


def analyse_mass(M, verbose=True):
    d, every = CELLS[M]
    fs = sorted(glob.glob(os.path.join(d, "red_*.csv")))
    if len(fs) < 10:
        return None
    dt = every / SPS
    alpha = M / (2.0 * NS)
    K = kroot(alpha)
    nu = CS * K / (2 * math.pi * LEFF); om = 2 * math.pi * nu
    per_kr = 1 / nu
    per_p1 = 1 / (CS * 1.0101 * K / (2 * math.pi * LEFF))
    per_id = 2 * math.pi * LEFF / (math.sqrt(2.0) * K)

    D, X = [], []
    for f in fs:
        e = pd.read_csv(f)
        D.append((e["KE_gas_left"].to_numpy(float) - e["KE_gas_right"].to_numpy(float)) / NS)
        X.append(e["W0_x_sigma"].to_numpy(float))
    n = min(len(a) for a in D); lo = int(BURN / dt)
    D = [a[lo:n] for a in D]; X = [a[lo:n] for a in X]
    nraw = len(D[0]); rec = nraw * dt
    tau_pred = A_HARD * M
    # ##CHRIS 2026-10-03, SECOND CORRECTION: the ACF fit window was sized from tau_PRED, which is
    # only valid if the prediction is right. Where tau is larger, the window covers far fewer tau
    # than intended -- 1.5 tau at M = 200 against the 6 intended -- and biases the fit high (27942
    # over a 2.2-tau window against 24510 over 13 tau, a 14 % effect). It is also wrong to simply
    # use the whole record: at low mass that fits an exponential over 20 tau, where the tail is
    # noise and the earlier campaigns showed a spurious 1.4x. So the window is made SELF-CONSISTENT:
    # start from the prediction, refit, resize to 8 tau_measured, and iterate. Converged in 2-3
    # passes at every mass.
    def _window(tau_guess):
        return max(64, int(min(8.0 * tau_guess, rec / 3) / dt))
    nl = _window(tau_pred); lag = np.arange(nl) * dt

    def model(t, A, tT, B, tr):
        return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(om * t)

    def fit(series):
        c = np.mean([acf(s, nl) for s in series], axis=0)
        p, _ = curve_fit(model, lag, c, p0=[0.5, tau_pred, 0.5, 0.4 * per_kr],
                         bounds=([0, 20, 0, 20], [1.5, 1e5, 1.5, 5e4]), maxfev=80000)
        return p

    def block(series):
        B = max(2, int(round(per_kr / dt)))
        bl = [s[:len(s) // B * B].reshape(-1, B).mean(1) for s in series]
        bdt = B * dt; k = min(len(bl[0]) // 2, 14)
        c = np.mean([acf(b, k) for b in bl], axis=0); lg = np.arange(k) * bdt
        m = c > 0.05
        if m.sum() < 3:
            return float("nan")
        s = np.polyfit(lg[m], np.log(c[m]), 1)[0]
        return -1 / s if s < 0 else float("nan")

    # ---- make the fit window self-consistent BEFORE calibrating ----
    for _ in range(4):
        try:
            _t = fit(D)[1]
        except Exception:
            break
        _nl = _window(_t)
        if abs(_nl - nl) <= max(2, 0.02 * nl):
            nl = _nl; lag = np.arange(nl) * dt; break
        nl = _nl; lag = np.arange(nl) * dt
    win_sigma = nl * dt

    # ---- calibration ----
    # ##CHRIS 2026-10-03, AFTER THE FIRST PASS: the grid must BRACKET the measurement, or the
    # inversion clamps at its edge and reports the edge as if it were a result. The pre-registered
    # grid was centred on the prediction (0.7-1.4 x 50.35 M); at M >= 50 the measured raw value
    # landed ABOVE it, np.interp saturated, and three masses came back at exactly 1.40 x prediction
    # -- the grid's top multiplier, not physics. So the raw ACF is fitted FIRST, with no
    # calibration, purely to set the RANGE the synthetics must cover. This changes the span of the
    # calibration, never the model, the estimator or the fit window: it is the difference between a
    # calibration that spans the answer and one that does not.
    rng = np.random.default_rng(1000 + M)
    A0, B0 = 0.50, 0.47
    TR0 = 0.4 * per_kr * 5
    try:
        raw_probe = fit(D)[1]
    except Exception:
        raw_probe = tau_pred
    hi = max(1.5, 2.2 * raw_probe / tau_pred)
    grid = sorted({0.7, 0.9, 1.0, 1.23, round(0.35 * hi, 3), round(0.6 * hi, 3),
                   round(0.85 * hi, 3), round(hi, 3)})
    tr_, rc_, sd_ = [], [], []
    for g in grid:
        true = g * tau_pred; got = []
        for _ in range(8):
            # ##CHRIS: generated with lfilter, not a Python loop. The loop form is the same maths
            # but ~100x slower, and at M = 200 (93k samples x 80 seeds x 48 trials) it does not
            # finish. AR(1): x[n] = a x[n-1] + e[n]; AR(2): v[n] = c1 v[n-1] + c2 v[n-2] + w[n].
            a = math.exp(-dt / true)
            rho = math.exp(-dt / TR0); c1 = 2 * rho * math.cos(om * dt); c2 = -rho * rho
            E = rng.normal(size=(len(fs), nraw)) * math.sqrt(1 - a * a)
            W = rng.normal(size=(len(fs), nraw))
            S1 = lfilter([1.0], [1.0, -a], E, axis=1)
            V = lfilter([1.0], [1.0, -c1, -c2], W, axis=1)
            V /= np.where(V[:, 200:].std(axis=1, keepdims=True) > 0,
                          V[:, 200:].std(axis=1, keepdims=True), 1.0)
            ser = list(math.sqrt(A0) * S1 + math.sqrt(B0) * V)
            try:
                got.append(fit(ser)[1])
            except Exception:
                pass
        if got:
            tr_.append(true); rc_.append(float(np.mean(got))); sd_.append(float(np.std(got, ddof=1)))
    tr_, rc_, sd_ = np.array(tr_), np.array(rc_), np.array(sd_)
    slope = float(np.polyfit(tr_, rc_, 1)[0]) if len(tr_) > 2 else float("nan")

    # ---- the real fits ----
    res = {}
    for nm, S in (("dT", D), ("x", X)):
        p = fit(S)
        js = []
        for i in range(len(S)):
            try:
                js.append(fit([S[j] for j in range(len(S)) if j != i]))
            except Exception:
                pass
        js = np.array(js); k = len(js)
        err = np.sqrt((k - 1) / k * ((js - js.mean(0)) ** 2).sum(0))
        inside = bool(rc_.min() <= p[1] <= rc_.max())
        inv = float(np.interp(p[1], rc_, tr_))
        etot = math.sqrt(err[1] ** 2 + float(np.interp(p[1], rc_, sd_)) ** 2) / max(slope, 1e-9)
        res[nm] = dict(p=p.tolist(), err=err.tolist(), tau_T=inv, tau_T_err=etot,
                       tau_r=float(p[3]), tau_r_err=float(err[3]), block=block(S),
                       in_calibration_range=inside)

    # ---- period and sigma ----
    def psd(sig):
        s = sig - sig.mean(); s = s * np.hanning(len(s))
        return np.fft.rfftfreq(len(s), dt), np.abs(np.fft.rfft(s)) ** 2
    fr, _ = psd(D[0])
    for nm, S in (("dT", D), ("x", X)):
        p = np.mean([psd(a)[1] for a in S], axis=0)
        band = (fr > 0.3 * nu) & (fr < 3 * nu); i = int(np.argmax(np.where(band, p, -1)))
        dd = 0.5 * (p[i-1] - p[i+1]) / (p[i-1] - 2*p[i] + p[i+1])
        fpk = fr[i] + dd * (fr[1] - fr[0])
        half = p[i] / 2; j = i
        while j > 0 and p[j] > half: j -= 1
        k2 = i
        while k2 < len(p) - 1 and p[k2] > half: k2 += 1
        sg = (fr[k2] - fr[j]) / 2.355
        res[nm]["period"] = 1 / fpk
        res[nm]["period_err"] = sg / fpk ** 2
    sdT = float(np.mean([a.std(ddof=1) for a in D]))

    tr_man, dff_man = mansour_tau_r(M, nu)
    Meff = KGAS2 / om ** 2
    out = dict(M=M, alpha=alpha, K=K, nu=nu, per_kr=per_kr, per_p1=per_p1, per_id=per_id,
               fit_window=win_sigma,
               seeds=len(fs), dt=dt, record=rec, L_over_tau=rec / tau_pred,
               slope=slope, cal_true=tr_.tolist(), cal_rec=rc_.tolist(), cal_sd=sd_.tolist(),
               sigma_dT=sdT, res=res, tau_r_mansour=tr_man, dff_mansour=dff_man,
               tau_r_friction=2 * Meff / GAMMA_FRIC, M_eff=Meff,
               tau_r_bulk=2.0 / (GAMMA_ENSKOG * (K / LEFF) ** 2))
    if verbose:
        print(f"\n=== M = {M}  (alpha = {alpha:.2f}, K = {K:.4f}, {len(fs)} seeds, "
              f"record {rec:.0f} sigma, dt = {dt:.2f}, self-consistent fit window {win_sigma:.0f}) ===")
        print(f"  sigma(T1-T2) = {sdT:.4f}  vs 0.199 ({100*(sdT/0.199-1):+.1f} %)")
        print(f"  period: KR {per_kr:.1f} | P1 {per_p1:.1f} | ideal {per_id:.1f}")
        for nm in ("dT", "x"):
            r = res[nm]
            print(f"    {nm:3s} measured {r['period']:.1f} +- {r['period_err']:.1f}  "
                  f"-> KR {abs(r['period']-per_kr)/r['period_err']:.1f}s, ideal {abs(r['period']-per_id)/r['period_err']:.1f}s")
        print(f"  calibration slope = {slope:.3f}  ({'ADOPT modelled' if slope >= 0.6 else 'BELOW 0.6 -- block is headline'})")
        for nm in ("dT", "x"):
            r = res[nm]
            print(f"    {nm:3s} tau_T = {r['tau_T']:.0f} +- {r['tau_T_err']:.0f}   "
                  f"vs {A_HARD*M:.0f} ({(r['tau_T']-A_HARD*M)/max(r['tau_T_err'],1e-9):+.1f}s), "
                  f"vs {A_IDEAL*M:.0f} ({(r['tau_T']-A_IDEAL*M)/max(r['tau_T_err'],1e-9):+.1f}s); "
                  f"block raw {r['block']:.0f}; tau_r = {r['tau_r']:.0f} +- {r['tau_r_err']:.0f}"
                  f"{'' if r['in_calibration_range'] else '   *** OUTSIDE CALIBRATION RANGE ***'}")
        print(f"  tau_r bounds: friction {out['tau_r_friction']:.0f} | Mansour {tr_man:.0f} | bulk {out['tau_r_bulk']:.0f}")
    return out


def main():
    res = {}
    for M in sorted(CELLS):
        r = analyse_mass(M)
        if r:
            res[M] = r
    if len(res) < 3:
        print(f"\nonly {len(res)} masses available -- ladder fit needs at least 3"); return
    json.dump(res, open(os.path.join(HERE, "261003_ladder_results.json"), "w"), indent=1)

    Ms = np.array(sorted(res))
    tau = np.array([np.mean([res[M]["res"][k]["tau_T"] for k in ("dT", "x")]) for M in Ms])
    err = np.array([np.mean([res[M]["res"][k]["tau_T_err"] for k in ("dT", "x")]) / math.sqrt(2) for M in Ms])

    # ---- ladder fit tau_T = a M^b ----
    w = 1.0 / (err / tau) ** 2
    A = np.vstack([np.log(Ms), np.ones_like(Ms, dtype=float)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * np.log(tau)))
    b, db = float(beta[0]), float(math.sqrt(cov[0, 0]))
    a, da = float(math.exp(beta[1])), float(math.exp(beta[1]) * math.sqrt(cov[1, 1]))
    print("\n" + "=" * 78)
    print(f"LADDER FIT  tau_T = a M^b  over {len(Ms)} masses {list(Ms)}")
    print(f"  b = {b:.4f} +- {db:.4f}   -> {abs(b-1)/db:.1f} sigma from the predicted b = 1")
    print(f"  a = {a:.2f} +- {da:.2f}   -> {abs(a-A_HARD)/da:.1f} sigma from {A_HARD} (hard disk), "
          f"{abs(a-A_IDEAL)/da:.1f} sigma from {A_IDEAL} (ideal)")
    # with b fixed at 1, the prefactor alone
    a1 = float(np.sum(w * tau / Ms) / np.sum(w));
    print(f"  with b fixed at 1: a = {a1:.2f}")

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    g = np.logspace(math.log10(Ms.min() * 0.8), math.log10(Ms.max() * 1.25), 50)
    ax.plot(g, A_HARD * g, "-", color=RED, lw=2.2, label=f"hard-disk isobar, $\\tau_T = {A_HARD}\\,M$")
    ax.plot(g, A_IDEAL * g, "--", color="k", lw=1.8, label=f"ideal gas, $\\tau_T = {A_IDEAL}\\,M$")
    ax.errorbar(Ms, tau, yerr=err, fmt="o", color=BLUE, ms=8, capsize=4, zorder=5,
                label=f"measured, 80 seeds/mass\nfit $b = {b:.3f} \\pm {db:.3f}$, $a = {a:.1f} \\pm {da:.1f}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("divider mass $M_d$"); ax.set_ylabel(r"$\tau_T$  [$\sigma$-time]")
    ax.set_title(r"Level 4 mass ladder: $\tau_T(M)$ against the two isobars")
    ax.grid(alpha=0.3, which="both"); ax.legend(frameon=False)
    for e_ in ("png", "pdf"):
        fig.savefig(os.path.join(T.PLOTS.replace("paper1_speedofsound", "paper2_energytransfer"),
                                 f"261003_p2_ladder_tauT.{e_}"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    al = np.array([res[M]["alpha"] for M in Ms])
    ga = np.logspace(math.log10(al.min() * 0.7), math.log10(al.max() * 1.4), 200)
    ax.plot(ga, [2 * math.pi * LEFF / (CS * kroot(x)) for x in ga], "-", color=RED, lw=2.2,
            label=r"$2\pi L_{\rm eff}/c_sK$, $\cot K = \alpha K$ (KR $c_s$)")
    ax.plot(ga, [2 * math.pi * LEFF / (math.sqrt(2) * kroot(x)) for x in ga], "--", color="k", lw=1.6,
            label="same with ideal-gas $c_s$")
    for nm, mk, col in (("dT", "o", BLUE), ("x", "s", ORANGE)):
        ax.errorbar(al, [res[M]["res"][nm]["period"] for M in Ms],
                    yerr=[res[M]["res"][nm]["period_err"] for M in Ms],
                    fmt=mk, color=col, ms=7, capsize=3.5, zorder=5,
                    label=("$T_1-T_2$" if nm == "dT" else "divider $x$"))
    ax.set_xscale("log"); ax.set_xlabel(r"$\alpha = M_d/(2N_sm)$")
    ax.set_ylabel(r"divider-mode period  [$\sigma$-time]")
    ax.set_title(r"Paper 1's equation across Paper 2's ladder")
    ax.grid(alpha=0.3, which="both"); ax.legend(frameon=False)
    for e_ in ("png", "pdf"):
        fig.savefig(os.path.join(T.PLOTS.replace("paper1_speedofsound", "paper2_energytransfer"),
                                 f"261003_p2_ladder_period.{e_}"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.plot(Ms, [res[M]["tau_r_mansour"] for M in Ms], "-", color=RED, lw=2.2,
            label=r"Mansour piston form, $\hat M = M + mN/3$")
    for nm, mk, col in (("dT", "o", BLUE), ("x", "s", ORANGE)):
        ax.errorbar(Ms, [res[M]["res"][nm]["tau_r"] for M in Ms],
                    yerr=[res[M]["res"][nm]["tau_r_err"] for M in Ms],
                    fmt=mk, color=col, ms=7, capsize=3.5, zorder=5,
                    label=("$T_1-T_2$" if nm == "dT" else "divider $x$"))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("divider mass $M_d$"); ax.set_ylabel(r"$\tau_r$  [$\sigma$-time]")
    ax.set_title(r"Mode damping across the ladder")
    ax.grid(alpha=0.3, which="both"); ax.legend(frameon=False)
    for e_ in ("png", "pdf"):
        fig.savefig(os.path.join(T.PLOTS.replace("paper1_speedofsound", "paper2_energytransfer"),
                                 f"261003_p2_ladder_taur.{e_}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\nfigures: 261003_p2_ladder_tauT / _period / _taur  (.png/.pdf)")


if __name__ == "__main__":
    main()
