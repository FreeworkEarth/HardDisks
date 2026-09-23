#!/usr/bin/env python3
"""##CHRIS 2026-10-02: the PRE-REGISTERED analysis of the M_d = 10 KOA-length equilibrium run.

Registered before the data existed (prompt of 2026-10-01, and the design in 260930_..._fd2.md):

  * modelled ACF estimator, C(t) = A e^(-t/tau_T) + B e^(-t/tau_r) cos(omega t), with omega FIXED
    from cot K = alpha K at alpha = M/(2 N_s m) = 0.1 and Kolafa-Rottner c_s;
  * its bias calibrated on synthetic OU + AR(2) AT THIS RECORD LENGTH AND SEED COUNT, and the
    calibration run BEFORE any real data is fitted -- which is why step 1 below is the synthetic;
  * the block estimator as an independent cross-check;
  * both observables, T_1 - T_2 and the divider position;
  * result quoted as calibrated tau_T with error against 504 (hard-disk isobar) and 618 (ideal),
    separation in sigma WHICHEVER WAY IT GOES;
  * plus sigma(T_1 - T_2) against 0.199, the mode period against 95.16 / 94.21 / 117.4 with the
    spectral FWHM as the error, and tau_r against the 22 / 4200 bounds and Mansour's 426.

Nothing here is tuned after seeing the answer. Run:  python3 paper2_level4_Md10_long_20261002.py
"""
import glob, math, os, sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos

REPO = os.path.dirname(os.path.dirname(HERE))
P = os.path.join(REPO, "hspist3", "experiments_energy_transfer",
                 "level4_equilibrium_KOAlength_20261002", "Md10")
NS = 50
DT = 200 / (612000 / 10195.000532)      # --trace-every=200 / (steps per sigma-time) = 3.332 sigma
BURN = 2000.0                            # discard to t > 2000, as in the 2026-09-29 analysis
HARD, IDEAL = 503.5, 618.4               # 50.35 M and 61.84 M at M = 10

# ---- the mode, from cot K = alpha K -----------------------------------------------------------
eta = 50 * math.pi * 0.25 / (39.25 * 10)
Z = float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0])
dZ = float(sos.dZ_kolafa_rottner_2006(np.array([eta]))[0])
CS = math.sqrt(Z + eta * dZ + Z * Z)
K = brentq(lambda k: math.cos(k) / math.sin(k) - 0.1 * k, 1e-9, math.pi - 1e-9)
LEFF = 38.75 - 1.0
NU = CS * K / (2 * math.pi * LEFF); OM = 2 * math.pi * NU; PERIOD = 1.0 / NU
PER_IDEAL = 2 * math.pi * LEFF / (math.sqrt(2.0) * K)
PER_P1 = 1.0 / (CS * 1.0101 * K / (2 * math.pi * LEFF))


def acf(x, nl):
    x = x - x.mean(); m = len(x)
    f = np.fft.rfft(x, 2 * m)
    c = np.fft.irfft(f * np.conj(f))[:nl].real
    return c / c[0]


def model(t, A, tT, B, tr):
    return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(OM * t)


def fit_modelled(series, nl, lag):
    c = np.mean([acf(s, nl) for s in series], axis=0)
    p, _ = curve_fit(model, lag, c, p0=[0.5, 500., 0.5, 200.],
                     bounds=([0, 20, 0, 20], [1.5, 2e4, 1.5, 5e3]), maxfev=60000)
    return p


def fit_block(series, B=None):
    """Block-average away the oscillation, then one exponential. Independent cross-check."""
    B = B or max(2, int(round(PERIOD / DT)))
    blocks = []
    for s in series:
        nb = len(s) // B
        blocks.append(s[:nb * B].reshape(nb, B).mean(1))
    bdt = B * DT; nl = min(len(blocks[0]) // 2, 14)
    c = np.mean([acf(b, nl) for b in blocks], axis=0)
    lag = np.arange(nl) * bdt
    m = c > 0.05
    if m.sum() < 3:
        return float("nan")
    k = np.polyfit(lag[m], np.log(c[m]), 1)[0]
    return -1 / k if k < 0 else float("nan")


def synth(tT, nraw, A0, B0, TR0, rng):
    a = math.exp(-DT / tT)
    s = np.zeros(nraw); s[0] = rng.normal()
    e = rng.normal(size=nraw) * math.sqrt(1 - a * a)
    for j in range(1, nraw):
        s[j] = a * s[j - 1] + e[j]
    rho = math.exp(-DT / TR0); c1 = 2 * rho * math.cos(OM * DT); c2 = -rho * rho
    v = np.zeros(nraw); w = rng.normal(size=nraw)
    for j in range(2, nraw):
        v[j] = c1 * v[j - 1] + c2 * v[j - 2] + w[j]
    v /= (v[200:].std() or 1.0)
    return math.sqrt(A0) * s + math.sqrt(B0) * v


def main():
    fs = sorted(glob.glob(os.path.join(P, "red_*.csv")))
    print(f"seeds: {len(fs)}   dt = {DT:.4f} sigma-time")
    D, X = [], []
    for f in fs:
        e = pd.read_csv(f)
        D.append((e["KE_gas_left"].to_numpy(float) - e["KE_gas_right"].to_numpy(float)) / NS)
        X.append(e["W0_x_sigma"].to_numpy(float))
    n = min(len(a) for a in D); lo = int(BURN / DT)
    D = [a[lo:n] for a in D]; X = [a[lo:n] for a in X]
    nraw = len(D[0]); rec = nraw * DT
    print(f"record after burn-in: {rec:.0f} sigma-time  ({nraw} samples), L/tau ~ {rec/HARD:.0f}")

    nl = int(min(6 * HARD, rec / 3) / DT); lag = np.arange(nl) * DT
    print(f"ACF fitted over lags 0-{nl*DT:.0f}\n")

    # ---------------- 1. sigma(T1-T2) -- the premise ----------------
    sd = float(np.mean([a.std(ddof=1) for a in D]))
    print(f"[1] sigma(T_1 - T_2) = {sd:.4f}   vs Beta(N,N) prediction 0.199  "
          f"({100*(sd/0.199-1):+.1f} %)")

    # ---------------- 2. the mode period ----------------
    def psd(sig):
        s = sig - sig.mean(); s = s * np.hanning(len(s))
        return np.fft.rfftfreq(len(s), DT), np.abs(np.fft.rfft(s)) ** 2
    fr, _ = psd(D[0])
    print(f"\n[2] mode period   predicted: KR {PERIOD:.2f} | Paper-1 c_s {PER_P1:.2f} | ideal {PER_IDEAL:.2f}")
    for nm, S in (("T1-T2", D), ("divider x", X)):
        p = np.mean([psd(a)[1] for a in S], axis=0)
        band = (fr > 0.004) & (fr < 0.020); i = int(np.argmax(np.where(band, p, -1)))
        d = 0.5 * (p[i-1] - p[i+1]) / (p[i-1] - 2*p[i] + p[i+1])
        fpk = fr[i] + d * (fr[1] - fr[0])
        half = p[i] / 2; j = i
        while j > 0 and p[j] > half: j -= 1
        k2 = i
        while k2 < len(p) - 1 and p[k2] > half: k2 += 1
        sig = (fr[k2] - fr[j]) / 2.355
        per = 1 / fpk; eper = sig / fpk ** 2
        print(f"    {nm:10s}: {per:7.2f} +- {eper:5.2f}   "
              f"KR {abs(per-PERIOD)/eper:4.1f}s | P1 {abs(per-PER_P1)/eper:4.1f}s | ideal {abs(per-PER_IDEAL)/eper:4.1f}s")

    # ---------------- 3. CALIBRATION FIRST, on synthetics ----------------
    print("\n[3] bias calibration of the modelled estimator, BEFORE fitting real data")
    rng = np.random.default_rng(101)
    A0, B0, TR0 = 0.47, 0.50, 190.0
    grid = [350, 450, 503.5, 550, 618.4, 700, 850]
    rec_tab = []
    NTR = 12
    for true in grid:
        g = []
        for _ in range(NTR):
            ser = [synth(float(true), nraw, A0, B0, TR0, rng) for _ in range(len(fs))]
            try:
                g.append(fit_modelled(ser, nl, lag)[1])
            except Exception:
                pass
        g = np.array(g)
        rec_tab.append((true, g.mean(), g.std(ddof=1)))
        print(f"    true {true:6.1f} -> recovered {g.mean():7.1f} +- {g.std(ddof=1):5.1f}")
    TR = np.array([r[0] for r in rec_tab]); RC = np.array([r[1] for r in rec_tab])
    SD = np.array([r[2] for r in rec_tab])
    slope = float(np.polyfit(TR, RC, 1)[0])
    print(f"    response slope = {slope:.3f}   (0.6 was the adoption threshold; "
          f"{'ADOPT' if slope >= 0.6 else 'below threshold'})")

    # ---------------- 4. the real fits ----------------
    print("\n[4] modelled fit on the data, omega fixed")
    out = {}
    for nm, S in (("T1-T2", D), ("divider x", X)):
        p = fit_modelled(S, nl, lag)
        js = []
        for i in range(len(S)):
            try:
                js.append(fit_modelled([S[j] for j in range(len(S)) if j != i], nl, lag))
            except Exception:
                pass
        js = np.array(js); k = len(js)
        err = np.sqrt((k - 1) / k * ((js - js.mean(0)) ** 2).sum(0))
        out[nm] = (p, err)
        print(f"    {nm:10s}: A={p[0]:.3f}+-{err[0]:.3f}  tau_T={p[1]:.0f}+-{err[1]:.0f}  "
              f"B={p[2]:.3f}+-{err[2]:.3f}  tau_r={p[3]:.0f}+-{err[3]:.0f}")

    print("\n[5] CALIBRATED tau_T, and the verdict")
    for nm in ("T1-T2", "divider x"):
        p, err = out[nm]
        meas = p[1]
        inv = float(np.interp(meas, RC, TR))
        e_tot = math.sqrt(err[1] ** 2 + float(np.interp(meas, RC, SD)) ** 2) / max(slope, 1e-9)
        print(f"    {nm:10s}: recovered {meas:.0f} -> tau_T = {inv:.0f} +- {e_tot:.0f}")
        for lbl, pr in (("hard disk 503.5", HARD), ("ideal 618.4", IDEAL)):
            cr = float(np.interp(pr, TR, RC)); cs_ = float(np.interp(pr, TR, SD))
            s_ = (meas - cr) / math.sqrt(err[1] ** 2 + cs_ ** 2)
            print(f"        vs {lbl:16s}: control recovers {cr:6.1f}, measured {meas:6.1f} -> {s_:+.1f} sigma")

    print("\n[6] block estimator, independent cross-check")
    for nm, S in (("T1-T2", D), ("divider x", X)):
        print(f"    {nm:10s}: tau_block = {fit_block(S):.0f}")

    print("\n[7] tau_r against the bounds")
    for nm in ("T1-T2", "divider x"):
        p, err = out[nm]
        tr = p[3]
        dff = 1 / (math.pi * tr * NU)
        print(f"    {nm:10s}: tau_r = {tr:.0f} +- {err[3]:.0f}   Df/f = {dff:.4f}, Q = {1/dff:.1f}")
    print("        bounds: kinetic friction 22 | Mansour piston 426 | bulk absorption 4217")


if __name__ == "__main__":
    main()
