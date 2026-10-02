#!/usr/bin/env python3
"""##CHRIS 2026-10-12: amendments A1-A3 to the efficiency-map pre-registration (261010 sec. 1.8).

Prints every table of sec. 1.8, before launch:
  0. a reproduction of the pre-registered sec. 1.3 reversible reference and sec. 1.5 run lengths
     (the original tables were computed inline, not by a committed script -- this closes that);
  1. tau_r per (k, M_s) from Mansour's friction with ONE gas, at the settled state;
  2. the amended run lengths max(d/u + 5P, d/u + 3 tau_r) (+ the 0.25/u gap), steps, trace rows;
  3. the cost, from a timed throwaway run of the binary (made here, into a temp dir);
  4. the floors that decide how the KE_div gate can be implemented;
  5. the SPECS block for level5_effmap_20261010.sh (--emit-specs).

Mansour, Garcia & Baras, PRE 73, 016121 (2006): Eq. 17 (p. 5) friction L_y Gamma / X_p per gas column
(two columns in their geometry; geometry C has ONE, so half the two-gas rate), M_hat = M + mN/3 (Eq. 18);
Enskog Gamma = eta_s + zeta from Eqs. 8-9 (p. 3), via paper1_linewidth_decomp_20261012.enskog.
Convention (methods sec. 13): tau_r is the AMPLITUDE decay time, energy rate Gamma_E = 2/tau_r.
"""
import math, os, subprocess, sys, tempfile, time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
from paper1_linewidth_decomp_20261012 import enskog

N, H, R, D, L0 = 100, 10.0, 0.5, 7.96, 78.4998
ETA0 = N * math.pi * R * R / (H * L0)
KS = (0.25, 0.5, 1.0); MS = (50, 200); US = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
XEQ_RUN = {0.25: 36.7996, 0.5: 33.65, 1.0: 32.0749}       # as in the runner (k = 0.5 keeps Level 3's 33.65)
STEPS_PER_SIGMA, HOLD = 60, 12000                        # 0.4 internal dt, 24 px/sigma; --steps counts post-hold
BIN = os.path.join(os.path.dirname(HERE), "00ALLINONE")

def Z(e): return float(sos.Z_kolafa_rottner_2006(np.array([e]))[0])
def dZ(e): return float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
def eta(L): return ETA0 * L0 / L
def T_ad(L):  # KR isentrope, d ln T = -Z d ln L
    return math.exp(-quad(lambda y: Z(eta(math.exp(y))), math.log(L0), math.log(L))[0])
def F(L): return N * T_ad(L) * Z(eta(L)) / L          # P h on the adiabat
def x_of(L, d): return 109.0 - d - L                     # sec. 1.3 geometry
def period(M, k): return 2 * math.pi * math.sqrt(M / k)
def steps(T): return int(math.ceil(STEPS_PER_SIGMA * T / 1000.0) * 1000)

def reversible():
    Ph = F(L0); out = {}
    for k in KS:
        xeq = 30.5 + Ph / k
        Lf = brentq(lambda L: k * (xeq - x_of(L, D)) - F(L), 60, L0)
        Eqs = N * (T_ad(Lf) - 1); s = 30.5 - x_of(Lf, D)
        Esp = 0.5 * k * (x_of(Lf, D) - xeq) ** 2 - 0.5 * k * (30.5 - xeq) ** 2
        out[k] = dict(xeq=xeq, Lf=Lf, s=s, Eqs=Eqs, Esp=Esp, W=Eqs + Esp, eps=Esp / (Eqs + Esp), Tf=T_ad(Lf))
    return Ph, out

def tau_r(k, M, rev):
    r = rev[k]; ef = eta(r["Lf"]); Tf = r["Tf"]
    _, _, es, ze, _ = enskog(ef, Tf)
    gam = H * (es + ze) / r["Lf"]                       # ONE gas column (Eq. 17, one term)
    Mhat = M + N / 3.0
    return dict(eta_f=ef, Tf=Tf, G=es + ze, gam=gam, Mhat=Mhat, GE=gam / Mhat, tr=2 * Mhat / gam)

def time_binary():
    with tempfile.TemporaryDirectory() as td:
        cmd = [BIN, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0",
               "--seed-drift-order=drift-first", f"--energy-transfer-summary={td}/s.csv",
               f"--energy-transfer-trace={td}/t.csv", "--trace-every=188", "--particles=100",
               "--particles-boxes=0,100", "--particle-radius=0.5", "--l0=54.75", "--height=10",
               "--num-walls=1", "--wall-positions=30.5", "--wall-mass-factors=200", "--spring-k-sigma=0.5",
               "--spring-wall=0", "--spring-eq=33.65", "--eff-output=spring",
               "--piston-right-protocol-mode=step", "--velocity-right-piston-step=0.05",
               "--max-right-piston-travel=7.96", "--auto-piston-step", f"--wall-hold-steps={HOLD}",
               "--steps=240000", "--fixed-dt=0.4", "--kbt1", "--seed=1"]
        t0 = time.perf_counter(); subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return (time.perf_counter() - t0) / (240000 + HOLD)

def main():
    emit = "--emit-specs" in sys.argv
    Ph, rev = reversible()
    TR = {(k, M): tau_r(k, M, rev) for k in KS for M in MS}
    cells = []
    for k in KS:
        for M in MS:
            P = period(M, k); every = int(STEPS_PER_SIGMA * P / 40); tr = TR[(k, M)]["tr"]
            for u in US:
                T_old = D / u + 5 * P; T_new = 0.25 / u + D / u + max(5 * P, 3 * tr)
                cells.append(dict(tag=f"k{k}_M{M}_u{u}", k=k, M=M, u=u, P=P, every=every, T_old=T_old,
                                  T_new=T_new, s_old=steps(T_old), s_new=steps(T_new), push=1))
            c = [x for x in cells if x["k"] == k and x["M"] == M and x["u"] == 0.01][0]
            cells.append(dict(c, tag=f"ctrl_k{k}_M{M}", push=0))
    if emit:
        for c in cells:
            print(f'  "{c["tag"]} {c["k"]} {XEQ_RUN[c["k"]]} {c["M"]} {c["u"]} {c["s_new"]} {c["every"]} {c["push"]}"')
        return
    print("### 0. Reproduction of the pre-registered sec. 1.3 (reversible reference)\n")
    print(f"eta_0 = {ETA0:.8f}, P h = F(L_0) = {Ph:.5f}\n")
    print("| k | x_eq | L_f | divider moves | E_qs(gas) | E_spring | W_rev | epsilon_rev | T_ad(L_f) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for k in KS:
        r = rev[k]
        print(f"| {k} | {r['xeq']:.4f} | {r['Lf']:.4f} | {r['s']:+.4f} | {r['Eqs']:.4f} | {r['Esp']:.4f} | "
              f"{r['W']:.4f} | {r['eps']:.4f} | {r['Tf']:.5f} |")
    print("\n### 1. tau_r per (k, M_s): Mansour friction, ONE gas column, at the settled state\n")
    print("gamma_1 = L_y (eta_s + zeta)/L_f (Eq. 17, one term = half the two-gas rate); M_hat = M_s + N m/3 (Eq. 18);")
    print("Enskog at eta_f = eta_0 L_0/L_f and T_f = T_ad(L_f) (KR isentrope); tau_r = 2 M_hat/gamma_1 (amplitude time).\n")
    print("| k | M_s | eta_f | T_f | eta_s + zeta | gamma_1 | M_hat | Gamma_E = gamma_1/M_hat | **tau_r** | 3 tau_r | 5 P (spring) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            t = TR[(k, M)]
            print(f"| {k} | {M} | {t['eta_f']:.5f} | {t['Tf']:.4f} | {t['G']:.4f} | {t['gam']:.5f} | {t['Mhat']:.2f} | "
                  f"{t['GE']:.3e} | **{t['tr']:.0f}** | {3*t['tr']:.0f} | {5*period(M, k):.0f} |")
    print("\nMansour is an UPPER bound on tau_r here: in all 7 cells of the 261006 ladder the measured damping")
    print("exceeded Mansour's friction (260913 REPORT, 1b), so 3 tau_r(Mansour) is the conservative length.\n")
    print("### 2. Run lengths and steps: old (d/u + 5P) vs amended (0.25/u + d/u + max(5P, 3 tau_r))\n")
    print("| cell | P | T old | steps old | T new | steps new | trace-every | rows/trace |")
    print("|---|---|---|---|---|---|---|---|")
    for c in cells:
        print(f"| {c['tag']} | {c['P']:.1f} | {c['T_old']:.0f} | {c['s_old']} | {c['T_new']:.0f} | {c['s_new']} | "
              f"{c['every']} | {c['s_new']//c['every']} |")
    old = 8 * sum(c["s_old"] for c in cells); new = 8 * sum(c["s_new"] for c in cells)
    sps = time_binary()
    print(f"\n### 3. Cost\n")
    print(f"timed here: one throwaway 240 000-step run of {os.path.basename(BIN)} in this geometry -> "
          f"**{1e6*sps:.2f} us per step** (incl. the {HOLD}-step hold)")
    print(f"total steps (8 seeds, 36 cells + 6 controls): old {old/1e6:.1f} M -> **new {new/1e6:.1f} M**")
    print(f"CPU: **{new*sps/3600:.2f} core-h**; wall at 9 jobs ~ **{new*sps/3600/9*60:.0f} min** "
          f"(longest single run {max(c['s_new'] for c in cells)*sps/60:.1f} min)")
    print("\n### 4. Floors that decide how the KE_div gate can be implemented\n")
    print("Equipartition gives the divider a thermal <KE_div> = kT_f/2 in its one degree of freedom, and the")
    print("instantaneous E_spring a thermal kT_f/2 on top of its static value. The coherent-energy floor of an")
    print("8-seed average is kT_f/(2 x 8). Sudden-limit bound on the coherent energy averaged over the last-tau_r")
    print("window of a run >= 3 tau_r long: (1 + k_S/k) e^-4 (1 - e^-2)/2 x Delta E_spring, k_S = N m c_s^2/L_f^2.\n")
    print("| k | Delta E_spring (rev) | 0.01 Delta E_spring | kT_f/2 (thermal KE_div) | ratio | kT_f/16 (8-seed floor) | k_S/k | sudden-limit E_coh/Delta E_spring |")
    print("|---|---|---|---|---|---|---|---|")
    for k in KS:
        r = rev[k]; ef = eta(r["Lf"]); kS = N * r["Tf"] * (Z(ef) + ef * dZ(ef) + Z(ef) ** 2) / r["Lf"] ** 2
        b = (1 + kS / k) * math.exp(-4) * (1 - math.exp(-2)) / 2
        print(f"| {k} | {r['Esp']:.4f} | {0.01*r['Esp']:.4f} | {r['Tf']/2:.4f} | {r['Tf']/2/(0.01*r['Esp']):.0f}x | "
              f"{r['Tf']/16:.4f} | {kS/k:.4f} | {b:.4f} |")
    print(f"\nA2 boundary: 2 L_0/c_s(eta_0) = {2*L0/math.sqrt(Z(ETA0)+ETA0*dZ(ETA0)+Z(ETA0)**2):.1f} sigma-time "
          f"-> tau_push = d/u below it for u > {D/(2*L0/math.sqrt(Z(ETA0)+ETA0*dZ(ETA0)+Z(ETA0)**2)):.4f}")

if __name__ == "__main__":
    main()
