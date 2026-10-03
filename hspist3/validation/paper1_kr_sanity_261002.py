#!/usr/bin/env python3
"""##CHRIS 2026-10-02 (Task L1): an INDEPENDENT evaluation of the Kolafa-Rottner hard-disk equation of state and of the
adiabatic sound speed built from it, against the project's module plot_speed_of_sound_edmd.py (tracked since 4f5fd31).
The coefficients below are typed from the paper, not copied from the module.

[SOURCE: J. Kolafa and M. Rottner, Mol. Phys. 104, 3435-3441 (2006), DOI 10.1080/00268970600967963; PDF on disk at
 ZZZ_PAPER/Kolafa and Rottner - 2006 - Simulation-based equation of state of the hard disk fluid ....pdf]
  - Eq. (7), p. 3437:  Z(y) = sum_{i=0}^{k} A_i (y/(1-y))^i,  y = A_HD rho = the packing fraction eta.
  - Section 3.2, p. 3438: "Note that x = y/(1 - y), where y is the packing fraction", and three fitted equations,
    for rho_max = 0.88 (s = 0.724, p. 3438), rho_max = 0.89 (s = 0.966, p. 3439) and rho_max = 0.90 (s = 0.927, p. 3439).
    rho = N sigma^2 / A is the reduced number density (p. 3435), so eta = (pi/4) rho and rho_max = 0.88 <=> eta = 0.6912.
  - The paper gives no table of A_i; the "table" is the three equations of section 3.2. Table 2 lists virial coefficients.
REFERENCE VERSION, fixed before the module was read: rho_max = 0.88 -- the project's documented validity "eta ~ 0.69"
(writeup/papers/README.md) is exactly its range, (pi/4)*0.88 = 0.6912. The other two versions are printed for information.

[DERIVATION] c_s for a 2D monatomic fluid: c_s^2 = (1/m)(dp/drho)_S with p = rho kT Z(eta), c_v = k per particle:
  (dp/drho)_T = kT (Z + eta Z'),  T (dp/dT)_rho^2 / (rho^2 c_v) = kT Z^2   =>   c_s^2 = (kT/m)(Z + eta Z' + Z^2).
Z' = dZ/deta = (dZ/dx)/(1 - eta)^2, evaluated analytically and, as a second check, by a complex-step derivative.
PASS: module and independent values agree to 1e-10 (relative) for Z, Z' and c_s at every eta. Nothing is edited.
usage: python3 hspist3/validation/paper1_kr_sanity_261002.py
"""
import inspect, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))

# powers of x -> coefficient, typed from section 3.2 of the paper
KR = {
    0.88: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504951668, 6: -0.0325433846,
           7: 0.0133946531, 8: 0.00174265604, 9: -0.00944632202, 10: 0.00851111768, 11: -0.0035963525,
           12: 0.000577345106, 19: -1.06399127e-7},
    0.89: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504963915, 6: -0.0325578581,
           7: 0.0134816028, 8: 0.00129187484, 9: -0.00808881628, 10: 0.00669011963, 11: -0.00250795961,
           12: 0.000336036442, 22: -5.15282664e-9, 57: 5.77730095e-23},
    0.90: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504960168, 6: -0.0325537792,
           7: 0.0134578632, 8: 0.00140888182, 9: -0.00834273601, 10: 0.00694127367, 11: -0.00262254723,
           12: 0.000355746352, 22: -5.24672938e-9, 57: 5.88054639e-23},
}
REF = 0.88
ETAS = (0.05, 0.30, 0.50, 0.65)

def Z(eta, v=REF):
    x = eta / (1 - eta); return sum(a * x ** i for i, a in KR[v].items())

def dZ(eta, v=REF):
    x = eta / (1 - eta); return sum(i * a * x ** (i - 1) for i, a in KR[v].items() if i) / (1 - eta) ** 2

def dZ_cstep(eta, v=REF, h=1e-30):
    return (Z(complex(eta, h), v)).imag / h

def cs(eta, z, dz): return math.sqrt(z + eta * dz + z * z)          # kT = m = 1

def rel(a, b): return abs(a - b) / max(abs(b), 1e-300)

def main():
    import plot_speed_of_sound_edmd as sos
    import tests_20260913 as T
    print("module functions:", ", ".join(f"{n}{inspect.signature(getattr(sos, n))}" for n in
          ("Z_kolafa_rottner_2006", "dZ_kolafa_rottner_2006", "cs_adiabatic_2d_monatomic")))
    print(f"independent reference: rho_max = {REF} version (eta_max = {math.pi / 4 * REF:.4f})\n")
    print("| eta | Z paper | Z module | rel diff | Z' paper (analytic) | Z' complex-step | Z' module | rel diff | "
          "c_s paper | c_s module | rel diff | tests_20260913.kr_cs | rel diff vs paper |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    worst = 0.0; worst_krcs = 0.0
    for e in ETAS:
        z, d, dc = Z(e), dZ(e), dZ_cstep(e); c = cs(e, z, d)
        zm = float(sos.Z_kolafa_rottner_2006(np.array([e]))[0]); dm = float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
        cm = float(sos.cs_adiabatic_2d_monatomic(np.array([zm]), np.array([dm]), np.array([e]), kbt=1, m=1)[0])
        kc = float(T.kr_cs(e))
        r = (rel(zm, z), rel(dm, d), rel(cm, c)); worst = max(worst, *r, rel(dc, d)); worst_krcs = max(worst_krcs, rel(kc, c))
        print(f"| {e:.2f} | {z:.15g} | {zm:.15g} | {r[0]:.1e} | {d:.15g} | {dc:.15g} | {dm:.15g} | {r[1]:.1e} | "
              f"{c:.15g} | {cm:.15g} | {r[2]:.1e} | {kc:.15g} | {rel(kc, c):.1e} |")
    print(f"\nlargest relative difference module vs paper (Z, Z', c_s; and analytic vs complex-step Z'): {worst:.1e}")
    note = " (above 1e-10; reported, not part of the verdict -- see how kr_cs forms Z')" if worst_krcs > 1e-10 else ""
    print(f"tests_20260913.kr_cs vs paper: largest relative difference {worst_krcs:.1e}{note}")
    print("\nfor information (NOT the verdict) -- the other two published versions against the module, largest relative "
          "difference over the four eta:")
    print("| version | Z | Z' (analytic) | c_s |\n|---|---|---|---|")
    for v in (0.89, 0.90):
        dz = dzp = dc = 0.0
        for e in ETAS:
            zm = float(sos.Z_kolafa_rottner_2006(np.array([e]))[0]); dm = float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
            cm = float(sos.cs_adiabatic_2d_monatomic(np.array([zm]), np.array([dm]), np.array([e]), kbt=1, m=1)[0])
            z, d = Z(e, v), dZ(e, v)
            dz, dzp, dc = max(dz, rel(zm, z)), max(dzp, rel(dm, d)), max(dc, rel(cm, cs(e, z, d)))
        print(f"| rho_max = {v} | {dz:.1e} | {dzp:.1e} | {dc:.1e} |")
    print("module comment, plot_speed_of_sound_edmd.py:620: '# Kolafa & Rottner (2006), rho_max=0.90 fit.  Their x is eta/(1-eta).'")
    ok = worst <= 1e-10
    print(f"\n**VERDICT: {'PASS' if ok else 'FAIL -- STOP, the module is NOT edited'}** (criterion: relative difference <= 1e-10)")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
