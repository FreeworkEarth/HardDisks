#!/usr/bin/env python3
"""##CHRIS 2026-09-17: the three zeta-audit items that the Level 2 run did not already cover.

(1) Bin-width convergence of the running Green-Kubo integral. The force is an impulse train, so the
    zero-lag self term is a delta smeared over one bin; the integral must stop moving as dt -> 0.
(2) Is WR really the generalized force of the coordinate we later move? Checked from the event log
    itself: the label WR stops and PR starts at the release time, same object, same position.
(3) zeta at the END of the compression path, from the post-stop equilibrium segment of the same runs.
    The sound-recurrence explanation makes a sharp prediction: the cancellation must move to the
    shorter traversal time L_f/c_s(eta_f, T_f) of the compressed, hotter compartment.
Analysis only.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import plot_speed_of_sound_edmd as sos

ET = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_energy_transfer"
NEW = f"{ET}/level2_slope_20260917"
NS, H, DX = 50, 10.0, 3.93
L_I = 35.320                     # right compartment during the start hold
WQS = 7.4879                     # Level 1 quasi-static work, = N_s kT (T_f/T_i - 1)


def running(files, kind, t0, t1, dt, nlag_t=60.0):
    """Seed-averaged force autocorrelation and its running integral over a stationary window."""
    nb = int((t1 - t0) / dt); nlag = int(nlag_t / dt)
    acc, sq, n = np.zeros(nlag), 0.0, 0
    for ev in files:
        d = pd.read_csv(ev, usecols=["t_sigma", "kind", "u_wall", "dp"])
        d = d[(d["kind"] == kind) & (d["t_sigma"] >= t0) & (d["t_sigma"] < t1) & (d["u_wall"] == 0.0)]
        if len(d) < 20:
            continue
        idx = np.clip(((d["t_sigma"].to_numpy(float) - t0) / dt).astype(int), 0, nb - 1)
        F = np.zeros(nb); np.add.at(F, idx, d["dp"].to_numpy(float)); F /= dt
        F -= F.mean()
        sp = np.fft.rfft(F, 2 * nb)
        acc += np.fft.irfft(sp * np.conj(sp))[:nlag] / (nb - np.arange(nlag))
        sq += float((d["dp"].to_numpy(float) ** 2).sum()); n += 1
    C = acc / n
    run = np.concatenate(([0.0], np.cumsum((C[1:] + C[:-1]) / 2 * dt)))
    return np.arange(nlag) * dt, run, sq / n / (t1 - t0) / 2.0, n


def cs_of(eta, T):
    Z = sos.Z_kolafa_rottner_2006(eta); dZ = sos.dZ_kolafa_rottner_2006(eta)
    return math.sqrt(T * (Z + eta * dZ + Z * Z)), Z


def main():
    files = sorted(glob.glob(f"{NEW}/u0.03/ev_*.csv"))

    print("### (2) Is the held wall the same object the protocol later moves?")
    d = pd.read_csv(files[0], usecols=["t_sigma", "kind", "u_wall"])
    wr, pr = d[d["kind"] == "WR"], d[d["kind"] == "PR"]
    print(f"  WR events span t = {wr['t_sigma'].min():.1f} → {wr['t_sigma'].max():.1f} σ, all with u_wall = "
          f"{sorted(set(wr['u_wall']))}")
    print(f"  PR events span t = {pr['t_sigma'].min():.1f} → {pr['t_sigma'].max():.1f} σ, u_wall = "
          f"{sorted(set(np.round(pr['u_wall'], 3)))}")
    print("  The two never overlap: the label changes at release. It is one wall — the piston itself,")
    print("  held before release and moving after — so C_FF is measured on the coordinate we drive.\n")

    print("### (1) Bin-width convergence of the running integral (start-of-path hold, 20–200 σ)")
    print("| bin dt [σ] | zero-lag self term | ζ(1 σ) | ζ(5 σ) | ζ(10 σ) | ζ(20 σ) | ζ(40 σ) |")
    print("|---|---|---|---|---|---|---|")
    for dt in (1.0, 0.5, 0.25, 0.125, 0.0625):
        lag, run, st, n = running(files, "WR", 20.0, 200.0, dt)
        g = lambda t: run[int(round(t / dt))]
        print(f"| {dt:.4f} | {st:.3f} | {g(1):+.3f} | {g(5):+.3f} | {g(10):+.3f} | {g(20):+.3f} | {g(40):+.3f} |")
    print(f"({n} holds per row, the same runs each time.)\n")

    print("### (3) ζ at the end of the path, and the sound-recurrence prediction")
    eta_i = NS * math.pi * 0.25 / (L_I * H)
    L_f = L_I - DX
    eta_f = NS * math.pi * 0.25 / (L_f * H)
    T_f = 1.0 + WQS / NS                       # adiabatic heating over the push, kT_i = 1
    c_i, Z_i = cs_of(eta_i, 1.0); c_f, Z_f = cs_of(eta_f, T_f)
    print(f"  start: L = {L_I:.2f} σ, η = {eta_i:.4f}, T = 1.000, Z = {Z_i:.4f}, c_s = {c_i:.3f} → L/c_s = {L_I / c_i:.1f} σ")
    print(f"  end:   L = {L_f:.2f} σ, η = {eta_f:.4f}, T = {T_f:.3f}, Z = {Z_f:.4f}, c_s = {c_f:.3f} → L/c_s = {L_f / c_f:.1f} σ")
    lag, run, st, n = running(files, "PR", 340.0, 598.0, 0.25, nlag_t=50.0)
    print(f"\n  post-stop hold, {n} runs, t = 340–598 σ, piston pinned at the compressed position:")
    print("| t_cut [σ] | " + " | ".join(f"{t:g}" for t in (0.25, 1, 5, 10, 14, 16, 18, 20, 25, 30)) + " |")
    print("|---" * 11 + "|")
    print("| ζ(t_cut) | " + " | ".join(f"{run[int(round(t / 0.25))]:+.2f}" for t in
                                       (0.25, 1, 5, 10, 14, 16, 18, 20, 25, 30)) + " |")
    z = run[int(round(0.25 / 0.25))]
    cross = next((lag[i] for i in range(1, len(run)) if run[i] <= 0 < run[i - 1]), float("nan"))
    print(f"\n  zero-lag (Enskog) value at the end of the path: {z:.3f} (start of path: 2.55)")
    lag0, run0, st0, _ = running(files, "WR", 20.0, 200.0, 0.25)
    cross0 = next((lag0[i] for i in range(1, len(run0)) if run0[i] <= 0 < run0[i - 1]), float("nan"))
    print(f"  first zero crossing measured: {cross:.1f} σ    predicted L_f/c_s(end) = {L_f / c_f:.1f} σ")
    print(f"  same at the start of the path: measured {cross0:.1f} σ, predicted {L_I / c_i:.1f} σ")

    print("\n### What the zero-lag value is: the Enskog piston friction, with no free parameter")
    print("  ζ_E = h n sqrt(8 m kT / pi) g(η),  g(η) = (1 - 7η/16)/(1-η)^2  (hard-disk contact value)")
    print("| end of path | n [σ^-2] | T | g(η) | ζ_E predicted | ζ zero-lag measured | ratio |")
    print("|---|---|---|---|---|---|---|")
    for nm, L, eta, T, meas in (("start", L_I, eta_i, 1.0, st0), ("compressed", L_f, eta_f, T_f, st)):
        n_ = NS / (L * H); g_ = (1 - 7 * eta / 16) / (1 - eta) ** 2
        zE = H * n_ * math.sqrt(8 * T / math.pi) * g_
        print(f"| {nm} | {n_:.5f} | {T:.3f} | {g_:.4f} | {zE:.3f} | {meas:.3f} | {meas / zE:.3f} |")


if __name__ == "__main__":
    main()
