#!/usr/bin/env python3
"""##CHRIS 2026-09-20: why the Level 3 pilot could not have passed, and where Level 3 has to live.

Pen and paper first, then the campaign. Two results, both derived here and both printed with the
numbers of the master box so the design is checked rather than asserted.

1. THE SERIES-SPRING CEILING. In the quasi-static limit the gas and the spring are two springs in
   series. The gas stiffness against a wall bounding a compartment of length L with N disks is

       F(L) = N kT Z(eta) / L ,   eta ~ 1/L
       k_gas = -dF/dx_w = (N kT / L^2) (Z + eta Z')

   A piston compression Delta x divides between them, the spring taking

       x_eq = Delta x * k_gas / (k + k_gas) ,   E_qs = 1/2 k x_eq^2 .

   d/dk of k/(k+k_gas)^2 vanishes at k = k_gas, so the capture is MAXIMAL at k = k_gas with

       E_max = k_gas Delta x^2 / 8 = N kT f^2 (Z + eta Z') / 8 ,   f = Delta x / L .

   The capture is SECOND order in the compression fraction f while W_qs is FIRST order. At f = 0.1
   that is a tenth of a kT however the spring is chosen -- and the wall carries kT/2 of thermal
   energy in that same coordinate. No choice of k rescues it. The 2026-09-19 pilot was not noisy,
   it was below the pedestal by construction.

2. WHERE IT DOES WORK. Level 2 measured what a stepped piston launches beyond W_qs: A u^2 with
   A_step = 25.9 +- 4.3. That is coherent acoustic energy, which is exactly the "ordered energy
   before it thermalizes" the ladder is about, and at u >= 0.2 it clears kT/2. The wall then acts
   as a mass-loaded acoustic termination; the share a pulse can deposit is the impedance transfer

       T_imp = 4 Z_w Z_g / (Z_w + Z_g)^2 ,   Z_g = rho c_s H ,   Z_w = M_s omega_w ,

   maximal at Z_w = Z_g. So Level 3 lives at u = 0.1 .. 1.0 with k scanned across BOTH special
   points -- k = k_gas (quasi-static optimum) and the k that matches impedance -- and the primary
   observable is the ENSEMBLE-MEAN E_spring(t), whose coherent peak rides above the kT/2 pedestal
   because the pedestal is incoherent and averages to a constant.
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

R, H, MS = 0.5, 10.0, 200.0
KT = 1.0


def gas(N, L):
    """Everything the gas contributes: eta, Z, k_gas, c_s, rho, acoustic impedance."""
    eta = N * math.pi * R * R / (L * H)
    Z = sos.Z_kolafa_rottner_2006(eta)
    dZ = sos.dZ_kolafa_rottner_2006(eta)
    k_gas = N * KT * (Z + eta * dZ) / L ** 2
    cs = math.sqrt(Z + eta * dZ + Z * Z)          # adiabatic, kT = m = 1
    rho = N / (L * H)
    return dict(eta=eta, Z=Z, dZ=dZ, stiff=Z + eta * dZ, k_gas=k_gas, cs=cs, rho=rho,
                Zg=rho * cs * H, L=L, N=N)


def e_qs(k, k_gas, dx):
    return 0.5 * k * (dx * k_gas / (k + k_gas)) ** 2


def t_imp(k, k_gas, Zg, Ms=MS):
    Zw = Ms * math.sqrt((k + k_gas) / Ms)          # M omega, omega = sqrt((k + k_gas)/M)
    return 4 * Zw * Zg / (Zw + Zg) ** 2, Zw


def main():
    print("## Level 3, designed on paper\n")

    # ---- 1. the pilot could not have passed -----------------------------------------------
    print("### 1. The series-spring ceiling\n")
    print("| geometry | N | L [σ] | η | Z + ηZ′ | k_gas | Δx | f | E_max = k_gas Δx²/8 |")
    print("|---|---|---|---|---|---|---|---|---|")
    pilot = gas(50, 38.75)
    newC = gas(100, 78.5)
    for lab, g, dx in (("pilot (2026-09-19)", pilot, 3.93), ("master box, geometry C", newC, 7.96)):
        f = dx / g["L"]
        em = g["k_gas"] * dx ** 2 / 8
        print(f"| {lab} | {g['N']} | {g['L']:.2f} | {g['eta']:.4f} | {g['stiff']:.4f} | "
              f"{g['k_gas']:.5f} | {dx:.2f} | {f:.4f} | **{em:.4f} kT** |")
    print()
    print(f"    cross-check of the closed form at the pilot: N f^2 (Z + eta Z')/8 = "
          f"{pilot['N'] * (3.93/pilot['L'])**2 * pilot['stiff'] / 8:.4f} kT")
    print(f"    the pilot ran k = 5, a factor {5/pilot['k_gas']:.0f} stiffer than k_gas, giving "
          f"E_qs = {e_qs(5, pilot['k_gas'], 3.93):.4f} kT")
    print(f"    the wall's thermal energy in that coordinate is kT/2 = 0.5, i.e. "
          f"{0.5/e_qs(5, pilot['k_gas'], 3.93):.0f}x the signal. NO k rescues it: the ceiling over "
          f"all k is {pilot['k_gas'] * 3.93**2 / 8:.4f} kT, still {0.5/(pilot['k_gas']*3.93**2/8):.0f}x below.")

    # what compression WOULD be needed
    print("\n    What quasi-static capture would need (E_max = N f^2 (Z + eta Z')/8 > 0.5 kT):\n")
    print("| N | f needed for E_max = kT/2 | η after that compression |")
    print("|---|---|---|")
    for N, L in ((100, 78.5), (500, 392.5)):
        g = gas(N, L)
        # solve N f^2 (Z + eta Z')/8 = 0.5 with the stiffness evaluated at the START density
        f = math.sqrt(0.5 * 8 / (N * g["stiff"]))
        eta_f = N * math.pi * R * R / ((L * (1 - f)) * H)
        print(f"| {N} | {f:.3f} | {eta_f:.4f} |")
    print("\n    (the stiffness rises during the push, so these f are upper bounds -- the point is")
    print("     the order of magnitude: tens of per cent, not ten.)")

    # ---- 2. the acoustic route -------------------------------------------------------------
    print("\n### 2. The acoustic route: what a stepped piston launches\n")
    A_STEP, dA = 25.91, 4.28
    g = newC
    print(f"    geometry C: N = {g['N']}, L = {g['L']} σ, η = {g['eta']:.4f}, "
          f"c_s = {g['cs']:.3f}, ρ = {g['rho']:.5f}, Z_g = ρ c_s H = {g['Zg']:.3f}")
    print(f"    k_gas = {g['k_gas']:.5f};  the impedance-matched stiffness solves M ω = Z_g:")
    k_match = MS * (g["Zg"] / MS) ** 2 - g["k_gas"]
    print(f"      M sqrt((k + k_gas)/M) = Z_g  ->  k = Z_g²/M - k_gas = {k_match:.5f}")
    print()
    print("| k | ω_w | T_w [σ-t] | Z_w = M ω_w | T_imp | E_qs (Δx = 7.96) | A u² at u = 0.2 / 0.5 / 1.0 |")
    print("|---|---|---|---|---|---|---|")
    KS = sorted({round(g["k_gas"], 5), 0.05, round(k_match, 4), 0.5, 5.0})
    for k in KS:
        T, Zw = t_imp(k, g["k_gas"], g["Zg"])
        w = math.sqrt((k + g["k_gas"]) / MS)
        tag = ""
        if abs(k - g["k_gas"]) < 1e-4: tag = "  ← k_gas"
        if abs(k - k_match) < 1e-3: tag = "  ← impedance match"
        print(f"| {k:.4f}{tag} | {w:.4f} | {2*math.pi/w:.1f} | {Zw:.2f} | {T:.3f} | "
              f"{e_qs(k, g['k_gas'], 7.96):.4f} | "
              + " / ".join(f"{T * A_STEP * u**2:.2f}" for u in (0.2, 0.5, 1.0)) + " |")
    print()
    print(f"    A_step = {A_STEP} ± {dA} (Level 2, stepped piston, per u²).")
    print(f"    Raw launched energy A u²: u = 0.2 -> {A_STEP*0.04:.2f} kT, 0.5 -> {A_STEP*0.25:.2f}, "
          f"1.0 -> {A_STEP:.2f} kT, against the kT/2 pedestal.")
    print("    T_imp × A u² is the OPTIMISTIC deposit: it assumes the whole excess arrives as one")
    print("    pulse at the wall. It is an upper estimate and is reported as such.")

    # ---- 3. the campaign that follows -------------------------------------------------------
    print("\n### 3. What that makes the campaign\n")
    print(f"    k ladder     {', '.join(f'{k:g}' for k in KS)}  (spans k_gas = {g['k_gas']:.4f} and")
    print(f"                 the impedance match {k_match:.4f}; 0.5 and 5 are the stiff wing)")
    print("    u ladder     0.1, 0.2, 0.5, 1.0   (A u² = "
          + ", ".join(f"{A_STEP*u**2:.2f}" for u in (0.1, 0.2, 0.5, 1.0)) + " kT)")
    print(f"    trace        200 σ-time after the piston stops; slowest wall period is "
          f"{2*math.pi/math.sqrt((min(KS)+g['k_gas'])/MS):.0f} σ-t")
    print(f"    sound        L/c_s = {g['L']/g['cs']:.1f} σ-t one way, so the first return is "
          f"{2*g['L']/g['cs']:.1f} σ-t after the step")
    print("    pass         E_coh > 3σ above the pedestal in every cell with u ≥ 0.2, and")
    print("                 eps_coh maximal at the k nearest the impedance match")


if __name__ == "__main__":
    main()
