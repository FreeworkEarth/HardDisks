#!/usr/bin/env python3
"""##CHRIS 2026-10-03: robustness checks that must pass BEFORE the ladder exponent is quotable.

b = 1.43 from the five-mass fit would be a departure from the adiabatic-piston theory in an
interacting fluid, so it is tested with checks designed to kill it:

  1. per-mass table, with the BLOCK estimator calibrated too (it was raw before)
  2. a third estimator with different systematics: the zero-frequency spectral density,
     S(0) = 4 sigma^2 tau for an OU process, sigma^2 known (0.199^2 for T1-T2)
  3. the exponent from the light subset {10,20,50} and the heavy subset {50,100,200} separately
  4. degeneracy: refit with tau_r FIXED at 0.5x / 1x / 2x the free value
  5. an independent physical route at M = 50: the 260928 temperature-step relaxation, no ACF at all

eta: the PHYSICAL packing fraction of the two-compartment box is
     100 pi r^2 / ((78.5 - t) * 10) = 0.100114, NOT the 0.10134 that a compartment length of
     38.75 = L0 - 2r would give -- 38.75 is a centre-accessible length and belongs in L_eff, not
     in an area. The tau_T prefactor is 50.48 M at the physical eta (50.35 M came from 0.1013).
"""
import glob, json, math, os, sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, brentq
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
import paper2_level4_ladder_20261003 as L

ETA_PHYS = 0.100114
Zv = float(sos.Z_kolafa_rottner_2006(np.array([ETA_PHYS]))[0])
dZv = float(sos.dZ_kolafa_rottner_2006(np.array([ETA_PHYS]))[0])
CS = math.sqrt(Zv + ETA_PHYS * dZv + Zv * Zv)
A_HARD = (4 / math.sqrt(2 * math.pi)) * 38.75 / (1 + ETA_PHYS * dZv / Zv)
A_IDEAL = (4 / math.sqrt(2 * math.pi)) * 38.75
SIG_DT = 0.199
NS = 50


def load(M):
    d, every = L.CELLS[M]
    dt = every / L.SPS
    fs = sorted(glob.glob(os.path.join(d, "red_*.csv")))
    D, X = [], []
    for f in fs:
        e = pd.read_csv(f)
        D.append((e["KE_gas_left"].to_numpy(float) - e["KE_gas_right"].to_numpy(float)) / NS)
        X.append(e["W0_x_sigma"].to_numpy(float))
    n = min(len(a) for a in D); lo = int(L.BURN / dt)
    return [a[lo:n] for a in D], [a[lo:n] for a in X], dt


def mode(M):
    K = L.kroot(M / 100.0)
    nu = CS * K / (2 * math.pi * L.LEFF)
    return K, nu, 2 * math.pi * nu, 1 / nu


def synth(tau, nraw, nseed, om, dt, TR0, A0, B0, rng):
    a = math.exp(-dt / tau)
    rho = math.exp(-dt / TR0); c1 = 2 * rho * math.cos(om * dt); c2 = -rho * rho
    E = rng.normal(size=(nseed, nraw)) * math.sqrt(1 - a * a)
    W = rng.normal(size=(nseed, nraw))
    S1 = lfilter([1.0], [1.0, -a], E, axis=1)
    V = lfilter([1.0], [1.0, -c1, -c2], W, axis=1)
    sd = V[:, 200:].std(axis=1, keepdims=True)
    V /= np.where(sd > 0, sd, 1.0)
    return list(math.sqrt(A0) * S1 + math.sqrt(B0) * V)


def acf_mean(series, nl):
    return np.mean([L.acf(s, nl) for s in series], axis=0)


def fit_model(series, nl, dt, om, p0, fix_tr=None):
    lag = np.arange(nl) * dt
    c = acf_mean(series, nl)
    if fix_tr is None:
        f = lambda t, A, tT, B, tr: A*np.exp(-t/tT) + B*np.exp(-t/tr)*np.cos(om*t)
        p, _ = curve_fit(f, lag, c, p0=p0, bounds=([0,20,0,20],[1.5,1e6,1.5,5e4]), maxfev=120000)
        return p
    f = lambda t, A, tT, B: A*np.exp(-t/tT) + B*np.exp(-t/fix_tr)*np.cos(om*t)
    p, _ = curve_fit(f, lag, c, p0=[p0[0], p0[1], p0[2]],
                     bounds=([0,20,0],[1.5,1e6,1.5]), maxfev=120000)
    return np.array([p[0], p[1], p[2], fix_tr])


def fit_block(series, per, dt):
    B = max(2, int(round(per / dt)))
    bl = [s[:len(s)//B*B].reshape(-1, B).mean(1) for s in series]
    bdt = B * dt; k = min(len(bl[0])//2, 16)
    c = np.mean([L.acf(b, k) for b in bl], axis=0); lg = np.arange(k) * bdt
    m = c > 0.05
    if m.sum() < 3: return float("nan")
    s = np.polyfit(lg[m], np.log(c[m]), 1)[0]
    return -1/s if s < 0 else float("nan")


def fit_S0(series, dt, nu, var_known):
    """tau from the zero-frequency spectral density. S(f) = 4 var tau/(1+(2 pi f tau)^2)."""
    sp = []
    for s in series:
        x = s - s.mean(); w = np.hanning(len(x)); x = x * w
        norm = (w**2).sum()
        p = np.abs(np.fft.rfft(x))**2 * 2.0 * dt / norm
        sp.append(p)
    S = np.mean(sp, axis=0)
    fr = np.fft.rfftfreq(len(series[0]), dt)
    m = (fr > 0) & (fr < 0.25 * nu)
    g = lambda f, tau: 4.0 * var_known * tau / (1.0 + (2*math.pi*f*tau)**2)
    p, cov = curve_fit(g, fr[m], S[m], p0=[1000.0], bounds=([1.0], [1e7]), maxfev=40000)
    return float(p[0]), float(math.sqrt(cov[0][0]))


def calibrate(M, kind, nraw, nseed, dt, om, per, TR0, A0, B0, raw, nl):
    """Bias-calibrate ONE estimator at this mass, on a grid that brackets the measurement."""
    rng = np.random.default_rng(7000 + M + (0 if kind == "model" else 1))
    hi = max(1.6, 2.4 * raw / (A_HARD * M))
    grid = sorted({0.8, 1.0, 1.25, round(0.4*hi,3), round(0.65*hi,3), round(0.85*hi,3), round(hi,3)})
    tr, rc, sd = [], [], []
    for g in grid:
        true = g * A_HARD * M; got = []
        for _ in range(6):
            ser = synth(true, nraw, nseed, om, dt, TR0, A0, B0, rng)
            try:
                v = fit_model(ser, nl, dt, om, [A0, true, B0, TR0])[1] if kind == "model" \
                    else fit_block(ser, per, dt)
                if np.isfinite(v): got.append(v)
            except Exception:
                pass
        if len(got) > 2:
            tr.append(true); rc.append(float(np.mean(got))); sd.append(float(np.std(got, ddof=1)))
    tr, rc, sd = np.array(tr), np.array(rc), np.array(sd)
    slope = float(np.polyfit(tr, rc, 1)[0]) if len(tr) > 2 else float("nan")
    inside = bool(rc.min() <= raw <= rc.max())
    val = float(np.interp(raw, rc, tr))
    err = float(np.interp(raw, rc, sd)) / max(slope, 1e-9)
    return val, err, slope, inside


def powerlaw(Ms, tau, err):
    Ms = np.asarray(Ms, float); tau = np.asarray(tau, float); err = np.asarray(err, float)
    w = 1.0 / (err / tau) ** 2
    A = np.vstack([np.log(Ms), np.ones_like(Ms)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * np.log(tau)))
    return float(beta[0]), float(math.sqrt(cov[0,0])), float(math.exp(beta[1])), \
           float(math.exp(beta[1]) * math.sqrt(cov[1,1]))


def main():
    prev = json.load(open(os.path.join(HERE, "261003_ladder_results.json")))
    out = {}
    print(f"eta_phys = {ETA_PHYS}, c_s = {CS:.5f}, tau_T = {A_HARD:.3f} M (hard disk), "
          f"{A_IDEAL:.3f} M (ideal)\n")
    print("### 1-2. per mass: three estimators, all calibrated where calibration applies\n")
    print("| M | L/tau_T actual | cal slope | tau_T modelled | tau_T block | tau_T S(0) | "
          "tau_r | B/(A+B) | in cal range |")
    print("|---|---|---|---|---|---|---|---|---|")
    for M in (10, 20, 50, 100, 200):
        D, X, dt = load(M)
        K, nu, om, per = mode(M)
        nraw = len(D[0]); rec = nraw * dt
        p_prev = prev[str(M)]["res"]["dT"]["p"]
        TR0 = p_prev[3]; A0, B0 = p_prev[0], p_prev[2]
        nl = max(64, int(min(8.0 * p_prev[1] * 1.2, rec / 3) / dt))
        pm = fit_model(D, nl, dt, om, [A0, p_prev[1], B0, TR0])
        raw_m = pm[1]
        raw_b = fit_block(D, per, dt)
        vm, em, slm, inm = calibrate(M, "model", nraw, len(D), dt, om, per, TR0, A0, B0, raw_m, nl)
        vb, eb, slb, inb = calibrate(M, "block", nraw, len(D), dt, om, per, TR0, A0, B0, raw_b, nl)
        s0, s0e = fit_S0(D, dt, nu, SIG_DT ** 2)
        frac = pm[2] / (pm[0] + pm[2])
        out[M] = dict(tau_model=vm, tau_model_err=em, slope=slm, in_range=inm,
                      tau_block=vb, tau_block_err=eb, slope_block=slb, in_range_block=inb,
                      tau_S0=s0, tau_S0_err=s0e, tau_r=float(pm[3]), frac=float(frac),
                      raw_model=float(raw_m), raw_block=float(raw_b), rec=rec, dt=dt, nu=nu)
        print(f"| {M} | {rec/vm:.0f} | {slm:.3f} | {vm:.0f} ± {em:.0f} | {vb:.0f} ± {eb:.0f} | "
              f"{s0:.0f} ± {s0e:.0f} | {pm[3]:.0f} | {frac:.3f} | "
              f"{'yes' if inm and inb else 'NO'} |")
        print(f"|   | record {rec:.0f}, {rec/(20*vm):.1f}x the 20-tau S(0) requirement | | | | | | | |")
    json.dump({str(k): v for k, v in out.items()},
              open(os.path.join(HERE, "261003_ladder_robustness.json"), "w"), indent=1)

    print("\n### 3. exponent from subsets\n")
    print("| set | estimator | b ± σ | σ from b=1 | a ± σ |")
    print("|---|---|---|---|---|")
    for label, ms in (("all five {10,20,50,100,200}", [10,20,50,100,200]),
                      ("light {10,20,50}", [10,20,50]),
                      ("heavy {50,100,200}", [50,100,200])):
        for est, kv, ke in (("modelled", "tau_model", "tau_model_err"),
                            ("block", "tau_block", "tau_block_err"),
                            ("S(0)", "tau_S0", "tau_S0_err")):
            t = [out[m][kv] for m in ms]; e = [max(out[m][ke], 0.02*out[m][kv]) for m in ms]
            if any(not np.isfinite(v) for v in t + e): continue
            b, db, a, da = powerlaw(ms, t, e)
            print(f"| {label} | {est} | {b:.3f} ± {db:.3f} | {abs(b-1)/db:.1f} | {a:.2f} ± {da:.2f} |")

    print("\n### 4. degeneracy: tau_r fixed at 0.5x / 1x / 2x the free value\n")
    print("| M | tau_r fixed | tau_T raw | shift vs free | quoted error | degenerate? |")
    print("|---|---|---|---|---|---|")
    for M in (100, 200):
        D, X, dt = load(M)
        K, nu, om, per = mode(M)
        nraw = len(D[0]); rec = nraw * dt
        p_prev = prev[str(M)]["res"]["dT"]["p"]
        nl = max(64, int(min(8.0 * p_prev[1] * 1.2, rec / 3) / dt))
        free = fit_model(D, nl, dt, om, [p_prev[0], p_prev[1], p_prev[2], p_prev[3]])
        for mult in (0.5, 1.0, 2.0):
            q = fit_model(D, nl, dt, om, [free[0], free[1], free[2], free[3]],
                          fix_tr=free[3] * mult)
            shift = q[1] - free[1]
            qerr = out[M]["tau_model_err"] * out[M]["slope"]
            print(f"| {M} | {free[3]*mult:.0f} ({mult}x) | {q[1]:.0f} | {shift:+.0f} | "
                  f"±{qerr:.0f} (raw) | {'YES' if abs(shift) > qerr else 'no'} |")

    print("\n### 5. independent route: the 260928 temperature-step relaxation at M_d = 50\n")
    step = os.path.join(L.ET, "level4_thermal_20260928", "Md50")
    fs = sorted(glob.glob(os.path.join(step, "tr_*.csv")))
    if not fs:
        print("  data not found")
    else:
        Ds, T0 = [], None
        for f in fs:
            e = pd.read_csv(f, usecols=["Time", "KE_gas_left", "KE_gas_right"])
            Ds.append((e["KE_gas_left"].to_numpy(float) - e["KE_gas_right"].to_numpy(float)) / NS)
            if T0 is None: T0 = e["Time"].to_numpy(float)
        n = min(len(a) for a in Ds)
        m = np.mean([a[:n] for a in Ds], axis=0); t = T0[:n]
        noise = float(np.mean([a[:n].std(ddof=1) for a in Ds])) / math.sqrt(len(Ds))
        pred = A_HARD * 50
        sel = (t > 0.15 * pred) & (t < 1.5 * pred) & (np.abs(m) > 3 * noise)
        print(f"  {len(fs)} seeds, noise on the mean {noise:.4f}, "
              f"{sel.sum()} samples with |signal| > 3x noise over 0.15-1.5 tau_pred")
        if sel.sum() > 8:
            k = np.polyfit(t[sel], np.log(np.abs(m[sel])), 1)
            tau = -1/k[0] if k[0] < 0 else float("nan")
            r = np.corrcoef(t[sel], np.log(np.abs(m[sel])))[0,1]
            print(f"  tau_T(step, M=50) = {tau:.0f}   (log-linear, r = {r:.3f})")
            print(f"  equilibrium route at M=50: {out[50]['tau_model']:.0f} ± {out[50]['tau_model_err']:.0f}")
            print(f"  hard-disk prediction {pred:.0f}")
        else:
            print("  signal does not exceed the noise over enough of the window -- "
                  "no independent tau from this route, as the 260928 report already concluded")


if __name__ == "__main__":
    main()
