#!/usr/bin/env python3
"""##CHRIS 2026-10-12: Paper 1 confinement campaign -- PRE-REGISTRATION tables (no runs).

Prints every table of 261012_paper1_confinement.md sec. 1: the cell geometry, the predictions of
hypotheses A / B / C, the held-divider (method A) derivative budget, and the cost per scan. Nothing
here launches anything; the campaign waits for the go.

Inputs read from disk (not typed):
  * per-face force noise: Level 3 F(L) c0 event logs (held wall, N = 100, H = 10, eta = 0.10005);
  * run cost: A1v2 run.log wall times at eta = 0.392699 and 0.112200 (N = 100);
  * anchors: the 260919 A1v2 table (eta = 0.392699) and CS_P1 = 1.0101 x KR at eta ~ 0.10
    (paper2_level4_mode_ladder_20261006.py, the value the 261010 prompt rounds to +1.0 %).
"""
import csv, glob, math, os, re, statistics as st, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)
import plot_speed_of_sound_edmd as sos
ET = os.path.join(os.path.dirname(HERE), "experiments_energy_transfer")

R, TW, GRID = 0.5, 0.05, 24                 # disk radius, Paper 1 wall thickness, 1/24-sigma grid
MASSES = T.A1_MASSES                        # alpha set = A1_MASSES / (2 x 50) = 0.5 ... 20
SEEDS_B = 25; NPER = 200                    # A1v2: 25 seeds, 200 oscillations per trace
DX_ANCHOR = {"0.10": 1.0101 - 1.0}          # fractional excess at H = 10 (CS_P1 in the 261006 script)
T_SEED_A = 5000.0                           # method A: max record per seed [sigma-time]
M_HOLD = 1e9                                # held divider mass factor (Level 3 convention)
BIAS_MAX, NOISE_MAX = 0.001, 0.009          # derivative budget: bias <= 0.1 %, noise <= 0.9 % -> < 1 %

def Z(e): return float(sos.Z_kolafa_rottner_2006(np.array([e]))[0])
def dZ(e): return float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
def cs2(e): return Z(e) + e * dZ(e) + Z(e) ** 2
def eta_of(Ns, H, L0): return Ns * math.pi * R * R / (H * L0)
def leff(L0): return L0 - 2 * R - TW / 2
def kroot(a): return brentq(lambda k: math.cos(k) / math.sin(k) - a * k, 1e-9, math.pi - 1e-9)
def ongrid(v): return abs(v * GRID - round(v * GRID)) < 1e-9
def g2(e): return (1 - 7 * e / 16) / (1 - e) ** 2      # Mansour Eq. 11 with n pi/4 = eta

def noise_eps0():
    """Relative per-face force noise x sqrt(T), Level 3 c0 (N = 100 one gas, H = 10)."""
    import pandas as pd
    F, Tr = [], []
    for f in sorted(glob.glob(os.path.join(ET, "level3_FofL_20260925", "c0", "ev_*.csv"))):
        e = pd.read_csv(edmd_acc_guard.guard(f), usecols=["t_sigma", "kind", "dp"]); d = e[e["kind"] == "D0"]
        t0, t1 = e["t_sigma"].min(), e["t_sigma"].max()
        F.append(d["dp"].abs().sum() / (t1 - t0)); Tr.append(t1 - t0)
    F = np.array(F); return F.std(ddof=1) / F.mean() * math.sqrt(np.mean(Tr)), len(F), float(np.mean(Tr)), float(F.mean())

def cost_rate():
    """ms of CPU per sigma-time per 100 particles, from A1v2 run.log wall times."""
    out = {}
    for tag, e in (("0p392699", 0.392699), ("0p112200", 0.112200)):
        secs = tot = 0.0
        for d in glob.glob(os.path.join(T.DROOT, f"eta_{tag}", "m_*")):
            s = open(edmd_acc_guard.guard(os.path.join(d, "run.log")), errors="ignore").read()
            secs += sum(int(x) for x in re.findall(r"##RUN .*?\((\d+) s\)", s))
            tot += sum(float(x) for x in re.findall(r"T=([\d.]+) sigma-time", s))
        out[e] = 1000 * secs / tot
    return out

def rate_at(e, rates):
    """Interpolate the cost per sigma-time in the Enskog collision frequency ~ eta g2(eta)."""
    (e1, r1), (e2, r2) = sorted(rates.items())
    f = lambda x: x * g2(x); return r1 + (r2 - r1) * (f(e) - f(e1)) / (f(e2) - f(e1))

def anchors():
    ref = [r for r in csv.DictReader(open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv")))
           if abs(float(r["eta"]) - 0.392699) < 1e-6][0]
    DX_ANCHOR["0.39"] = float(ref["dev_KR_pct"]) / 100
    return {"0.10": dict(L0=39.25, H=10.0, Ns=50, sig=0.005294 / 1.81155),   # A1v2 eta=0.1122 relative scaled error
            "0.39": dict(L0=10.0, H=10.0, Ns=50, sig=float(ref["c_s_err_scaled"]) / float(ref["c_s"]))}

def cells(A):
    out = []
    for lab, a in A.items():
        for H in (5.0, 10.0, 20.0, 40.0):
            out.append(dict(scan="H", eta_lab=lab, H=H, L0=a["L0"], Ns=int(round(a["Ns"] * H / a["H"]))))
        for f in (0.5, 2.0):
            out.append(dict(scan="L", eta_lab=lab, H=a["H"], L0=a["L0"] * f, Ns=int(round(a["Ns"] * f))))
        area = a["H"] * a["L0"]
        for r in (1, 2, 4, 8):
            L0 = round(math.sqrt(area * r) * GRID) / GRID; H = round(math.sqrt(area / r) * GRID) / GRID
            out.append(dict(scan="aspect", eta_lab=lab, H=H, L0=L0, Ns=a["Ns"]))
    return out

def main():
    eps0, nse, Tr, F0 = noise_eps0(); rates = cost_rate(); A = anchors(); C = cells(A)
    print("### Inputs read from disk\n")
    print(f"- per-face force noise, Level 3 c0 ({nse} seeds, {Tr:.0f} sigma-time each, F = {F0:.4f}): "
          f"sigma_F/F = **{eps0:.4f}/sqrt(T)** for N = 100 behind the face, H = 10, eta = 0.10005")
    for e, r in sorted(rates.items()):
        print(f"- CPU cost, A1v2 run.log (N = 100): eta = {e}: **{r:.3f} ms per sigma-time**")
    print(f"- anchors (fractional c_s excess over KR at H = 10, N_s = 50): eta ~ 0.10: {100*DX_ANCHOR['0.10']:+.2f} %; "
          f"eta = 0.392699: {100*DX_ANCHOR['0.39']:+.3f} % (260919 table)")
    print(f"- anchor relative error on c_s (scaled): eta ~ 0.10: {100*A['0.10']['sig']:.3f} % (A1v2 eta = 0.1122 used as proxy); "
          f"eta = 0.392699: {100*A['0.39']['sig']:.3f} %\n")

    print("### Table G -- cell geometry (gas | divider | gas; r = 0.5, t = 0.05; eta = N_s pi r^2 / (H L_0))\n")
    print("| scan | eta anchor | H | L_0 | N_s per side | eta exact | L_eff | L_0/H | grid-exact (1/24) | masses M (alpha = 0.5 ... 20) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for c in C:
        c["eta"] = eta_of(c["Ns"], c["H"], c["L0"]); c["Le"] = leff(c["L0"])
        c["Ms"] = [int(round(M * c["Ns"] / 50)) for M in MASSES]
        g = "yes" if ongrid(c["L0"]) and ongrid(c["H"]) else "NO"
        print(f"| {c['scan']} | {c['eta_lab']} | {c['H']:.4f} | {c['L0']:.4f} | {c['Ns']} | {c['eta']:.6f} | {c['Le']:.4f} | "
              f"{c['L0']/c['H']:.3f} | {g} | {c['Ms'][0]} ... {c['Ms'][-1]} |")

    print("\n### Table P -- predicted fractional c_s excess over KR, by hypothesis\n")
    print("A: Delta = a (2/H + 2/L_0), a fixed by the H = 10 anchor.  B: Delta = Delta_10 (10/H), no L dependence.")
    print("C (INFERENCE, heavy-divider estimate, no free parameter): thermal-amplitude anharmonicity,")
    print("   Delta_C = (q+1)(q+2) / (16 N_s q Z),  q = c_s^2/(Z kT/m) = (Z + eta Z' + Z^2)/Z.\n")
    print("| scan | eta | H | L_0 | N_s | A [%] | B [%] | C [%] | (A-B)/sigma | (B-C_scaled)/sigma | sigma assumed [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for c in C:
        a0 = A[c["eta_lab"]]; d10 = DX_ANCHOR[c["eta_lab"]]
        aA = d10 / (2 / a0["H"] + 2 / a0["L0"])
        c["dA"] = aA * (2 / c["H"] + 2 / c["L0"]); c["dB"] = d10 * a0["H"] / c["H"]
        e = c["eta"]; z = Z(e); q = cs2(e) / z
        c["dC"] = (q + 1) * (q + 2) / (16 * c["Ns"] * q * z)
        c["q"] = q; c["Z"] = z
        dC_scaled = d10 * a0["Ns"] / c["Ns"]       # C's SHAPE (~1/N_s) anchored like B, for the separation column
        s = a0["sig"]
        print(f"| {c['scan']} | {e:.6f} | {c['H']:.3f} | {c['L0']:.3f} | {c['Ns']} | {100*c['dA']:+.3f} | {100*c['dB']:+.3f} | "
              f"{100*c['dC']:+.3f} | {(c['dA']-c['dB'])/s:+.1f} | {(c['dB']-dC_scaled)/s:+.1f} | {100*s:.3f} |")

    print("\n### Table K -- bulk thermodynamic ratios at the two anchors (KR)\n")
    print("| eta | Z | eta Z' | c_s^2 | q = c_s^2/Z | gamma = 1 + Z^2/(Z+eta Z') | k_T/k_S (bulk) |")
    print("|---|---|---|---|---|---|---|")
    for lab, a in A.items():
        e = eta_of(a["Ns"], a["H"], a["L0"]); z = Z(e); dz = e * dZ(e)
        print(f"| {e:.6f} | {z:.5f} | {dz:.5f} | {cs2(e):.5f} | {cs2(e)/z:.4f} | {1+z*z/(z+dz):.5f} | {(z+dz)/cs2(e):.5f} |")

    print("\n### Table I -- the identity test per cell: expected sigma, and what hypothesis C would do to it\n")
    print("rho_I = [N_s m c_s^2/L_eff^2 - k_T - F^2/(N_s kT)] / (N_s m c_s^2/L_eff^2), predicted 0.")
    print("sigma(rho_I)^2 = (2 sigma_cs)^2 + ((k_T/k_S) sigma_kT)^2, with sigma_cs = the anchor's relative error and")
    print(f"sigma_kT = {100*NOISE_MAX:.1f} % (Table A budget, noise part); the F^2 term's own error (< 0.05 %) is neglected.")
    print("Under C the mode frequency carries the thermal-amplitude shift and the local k_T does not: rho_I ~ 2 Delta_C.\n")
    print("| scan | eta | H | L_0 | N_s | k_T/k_S (bulk) | sigma(rho_I) [%] | rho_I under C [%] | rho_I(C)/sigma |")
    print("|---|---|---|---|---|---|---|---|---|")
    for c in C:
        e = c["eta"]; z = c["Z"]; dz = e * dZ(e); kr = (z + dz) / cs2(e); s = A[c["eta_lab"]]["sig"]
        sI = math.hypot(2 * s, kr * NOISE_MAX); rC = 2 * c["dC"]
        print(f"| {c['scan']} | {e:.4f} | {c['H']:.3f} | {c['L0']:.3f} | {c['Ns']} | {kr:.4f} | {100*sI:.3f} | "
              f"{100*rC:+.3f} | {rC/sI:+.1f} |")

    print("\n### Table A -- held divider (method A): stencil, derivative budget, record, cost\n")
    print(f"Rule: delta_L = max(1/24, grid-rounded sigma_x/2), sigma_x = (kT L_eff^2/(2 N_s m c_s^2))^(1/2) the free divider's")
    print(f"thermal rms excursion, so the stencil L_0 + {{0, +-dL, +-2dL}} spans what the free divider samples.")
    print(f"k_T = -[F(-2) - 8F(-1) + 8F(+1) - F(+2)]/(12 dL); noise factor sqrt(130)/12; per L-point the two faces of the")
    print(f"mirror runs (+x, -x) are averaged. Noise model eps(N_s) = eps0 (100/N_s)^(1/2) (ASSUMED; calibrated at eta = 0.39")
    print(f"by a 4-seed pilot before launch). Seeds of <= {T_SEED_A:.0f} sigma-time, held mass {M_HOLD:.0e}.\n")
    print("| scan | eta | H | L_0 | N_s | sigma_x | delta_L | 2dL/L_0 | bias (KR model) | T per position needed | seeds/position | divider drift/seed | CPU [core-h] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in C:
        e = c["eta"]; cs = math.sqrt(cs2(e)); z = c["Z"]; dz = e * dZ(e)
        sx = math.sqrt(c["Le"] ** 2 / (2 * c["Ns"] * cs * cs)); dL = max(1 / GRID, round(sx / 2 * GRID) / GRID)
        Ns, H = c["Ns"], c["H"]
        F = lambda L: Ns * Z(eta_of(Ns, H, L)) / L
        exact = Ns * (Z(eta_of(Ns, H, c["L0"])) + eta_of(Ns, H, c["L0"]) * dZ(eta_of(Ns, H, c["L0"]))) / c["L0"] ** 2
        x0 = c["L0"]
        st5 = -(F(x0 - 2 * dL) - 8 * F(x0 - dL) + 8 * F(x0 + dL) - F(x0 + 2 * dL)) / (12 * dL)
        bias = st5 / exact - 1
        eps = eps0 * math.sqrt(100 / Ns); FoverK = F(x0) / exact
        Tpos = ((math.sqrt(130) / 12) * eps * FoverK / (dL * NOISE_MAX)) ** 2 / 2
        nseed = math.ceil(Tpos / T_SEED_A)
        nu_c = 2 * H * (4 * e / math.pi) * z * math.sqrt(1 / (2 * math.pi))          # impacts per sigma-time, both faces
        drift = math.sqrt(nu_c * T_SEED_A * 8) * T_SEED_A / (math.sqrt(3) * M_HOLD)
        cpu = 5 * nseed * T_SEED_A * rate_at(e, rates) * (2 * Ns / 100) / 1000 / 3600
        c.update(sx=sx, dL=dL, cpuA=cpu, nseed=nseed, Tpos=Tpos)
        print(f"| {c['scan']} | {e:.4f} | {c['H']:.3f} | {c['L0']:.3f} | {Ns} | {sx:.3f} | {dL:.4f} ({round(dL*GRID)}/24) | "
              f"{2*dL/c['L0']:.4f} | {bias:+.1e} | {Tpos:.3g} | {nseed} | {drift:.1e} | {cpu:.2f} |")

    print("\n### Table B -- free divider (method B): 9 masses x 25 seeds x 200 periods, cost\n")
    print("| scan | eta | H | L_0 | N_s | period range (alpha 0.5 ... 20) | CPU [core-h] | wall [h] at 9 jobs |")
    print("|---|---|---|---|---|---|---|---|")
    for c in C:
        e = c["eta"]; cs = math.sqrt(cs2(e)); r = rate_at(e, rates)
        pers = [2 * math.pi * c["Le"] / (cs * kroot(M / (2 * c["Ns"]))) for M in c["Ms"]]
        cpu = sum(SEEDS_B * NPER * p * r * (2 * c["Ns"] / 100) / 1000 / 3600 for p in pers)
        c["cpuB"] = cpu
        print(f"| {c['scan']} | {e:.4f} | {c['H']:.3f} | {c['L0']:.3f} | {c['Ns']} | {min(pers):.1f} ... {max(pers):.1f} | {cpu:.2f} | {cpu/9:.2f} |")

    print("\n### Cost per scan (A + B; the H = 10 anchor cell is counted once, in the H-scan)\n")
    print("| scan | eta anchor | cells | CPU A [core-h] | CPU B [core-h] | total [core-h] | wall [h] at 9 jobs |")
    print("|---|---|---|---|---|---|---|")
    for scan in ("H", "L", "aspect"):
        for lab in ("0.10", "0.39"):
            cc = [c for c in C if c["scan"] == scan and c["eta_lab"] == lab]
            if scan == "aspect":   # L_0/H = 1 at eta = 0.39 IS the anchor cell -- not rerun
                cc = [c for c in cc if not (abs(c["L0"] - A[lab]["L0"]) < 1e-9 and abs(c["H"] - A[lab]["H"]) < 1e-9)]
            a_ = sum(c["cpuA"] for c in cc); b_ = sum(c["cpuB"] for c in cc)
            print(f"| {scan} | {lab} | {len(cc)} | {a_:.1f} | {b_:.1f} | {a_+b_:.1f} | {(a_+b_)/9:.1f} |")
    return C, A

def amend_c1():
    """##CHRIS 2026-10-12, amendment C1 (261012 sec. 1.9): the identity in its length-free form."""
    from paper1_populate_cs_err_20261002 import cell
    A = anchors(); eps0 = None
    print("### Table L -- the length convention the c_s/L form would depend on\n")
    print("| anchor | eta | L_0 (geometric) | L_eff = L_0 - 2r - t/2 | (L_0/L_eff)^2 |")
    print("|---|---|---|---|---|")
    for lab, a in A.items():
        print(f"| {lab} | {eta_of(a['Ns'], a['H'], a['L0']):.6f} | {a['L0']:.4f} | {leff(a['L0']):.4f} | {(a['L0']/leff(a['L0']))**2:.4f} |")
    print("\n### Table H -- heavy-divider form vs the exact standing wave, per alpha (box-independent)\n")
    print("Exact: omega^2 = c_s^2 K^2/L^2 with cot K = alpha K. Heavy form: omega^2 = 2 k_S / M_hat, k_S = N_s m c_s^2/L^2,")
    print("M_hat = M + 2 N_s m/3, i.e. omega^2 = (c_s^2/L^2)/(alpha + 1/3). The plan's M + N_s m/3 is shown for comparison.\n")
    print("| alpha | K | heavy/exact omega^2, M_hat = M + 2N_s m/3 | same with M + N_s m/3 | used in primary |")
    print("|---|---|---|---|---|")
    for M in MASSES:
        al = M / 100.0; K = kroot(al)
        r2 = 1 / (K * K * (al + 1 / 3)); r1 = 1 / (K * K * (al + 1 / 6))
        print(f"| {al:g} | {K:.5f} | {r2:.5f} | {r1:.5f} | {'yes' if al >= 5 else 'check only'} |")
    print("\n### Table I-omega -- expected sigma of the identity residual in the length-free form\n")
    print("Primary (alpha >= 5): k_S^dyn = M_hat omega_1^2 / 2 per mass, inverse-variance mean over the five heavy masses;")
    print("sigma(k_S^dyn)/k_S = 2 sigma_nu/nu (M_hat exact). Per-mass sigma_nu/nu = seed SE of the A1v2 cell (canonical")
    print("estimator) at the anchor (eta = 0.1122 stands in for 0.10). Static side as in Table I (k_T noise 0.9 %).")
    print("Standing-wave check (all alpha): k_S^SW = N_s m omega_1^2 / K(alpha)^2, per mass.\n")
    print("| anchor | per-mass 2 sigma_nu/nu, alpha = 0.5 ... 20 [%] | heavy combined 2 sigma_nu/nu [%] | k_T/k_S | sigma(rho_I) omega-form [%] | sigma(rho_I) c_s/L form (Table I) [%] | rho_I under C [%] | rho_I(C)/sigma |")
    print("|---|---|---|---|---|---|---|---|")
    for lab, tag, eA in (("0.10", "eta_0p112200", 0.112200), ("0.39", "eta_0p392699", 0.392699)):
        a = A[lab]; L0d = {"0.112200": 34.9999, "0.392699": 10.0}[f"{eA:.6f}"]
        rel = []
        for M in MASSES:
            c = cell((eA, L0d, M, T.cell_runs(os.path.join(T.DROOT, tag, f"m_{M}"), M)))
            rel.append(2 * c["sd"] / math.sqrt(c["n"]) / c["nu"])
        heavy = [r for M, r in zip(MASSES, rel) if M / 100.0 >= 5]
        comb = 1 / math.sqrt(sum(1 / r ** 2 for r in heavy))
        e = eta_of(a["Ns"], a["H"], a["L0"]); z = Z(e); dz = e * dZ(e); kr = (z + dz) / cs2(e)
        sI = math.hypot(comb, kr * NOISE_MAX); sI_old = math.hypot(2 * a["sig"], kr * NOISE_MAX)
        q = cs2(e) / z; dC = (q + 1) * (q + 2) / (16 * a["Ns"] * q * z)
        print(f"| {lab} | {' / '.join(f'{100*r:.2f}' for r in rel)} | {100*comb:.3f} | {kr:.4f} | {100*sI:.3f} | {100*sI_old:.3f} | "
              f"{200*dC:+.3f} | {2*dC/sI:+.1f} |")

if __name__ == "__main__":
    if "--c1" in sys.argv: amend_c1()
    else: main()
