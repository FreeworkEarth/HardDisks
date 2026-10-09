#!/usr/bin/env python3
"""##CHRIS 2026-10-04 (Task X): the PRE-REGISTERED analysis of the Paper 1 confinement campaign, applied ONCE.

Registration: 0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md sec. 1 with amendments C1-C3. Nothing here is
a new estimator; every quantity is the one the registration names, computed by the code it names:
  (B) free divider -> c_s: the canonical estimator (paper1_populate_cs_err_20261002: per-(cell, M) mean of the per-run
      peak frequency, TD = 200, X_EDGE = 2.5, which reduce_B.py wrote per run into red_nu.csv; through-origin
      slope_with_errors over the nine masses; x = K(alpha)/(2 pi L_eff,true), cot K = alpha K). alpha = M/(2 N_s m) with
      the CELL's N_s (sec. 1.4: "the masses are scaled with N_s ... every cell runs the same set alpha = 0.5 ... 20");
      the A1v2 helper T.x_of hard-codes N_s = 50 and is therefore not called. L_eff,true = L_0 - 2r - t/2 - delta/2 and
      eta_true = eta_rec L_0/(L_0 - delta/2) with delta = box_delta(L_0) (methods sec. 14; 0 for these grid-exact boxes).
      Delta = c_s/c_s^KR(eta_true) - 1, sigma = c_s_err_scaled/c_s^KR (sec. 1.5 "weighted by the scaled sigma_i").
  Confinement verdict (sec. 1.5), per eta, on the campaign's own cells: A: 2/H + 2/L_0, B: 1/H, C: 1/N_s, each with ONE
      free amplitude; excluded if p(chi2) < 0.01; one survivor = result; several = "not separated" (Delta chi2 given);
      none = the two-term forms b/H + c'/N_s and a(2/H + 2/L_0) + c'/N_s, exploratory. C also at its fixed amplitude
      Delta_C = (q+1)(q+2)/(16 N_s q Z), q = (Z + eta Z' + Z^2)/Z.
  (A) held divider -> F, k_T, T: per position the seed means of F_L, F_R (reduce_A.py); F(L_0 + x) = mean of F_L of run +x
      and F_R of run -x; k_T = -[F(-2) - 8F(-1) + 8F(+1) - F(+2)]/(12 dL); kT = KE/N on the x = 0 seeds. Symmetry checks
      (each within 2 sigma): F_L(L_0) = F_R(L_0), and F_L(L_0 + x) of run +x = F_R(L_0 + x) of run -x.
  Identity (amendment C1, primary): k_S^dyn = M_hat omega_1^2/2, M_hat = M + (2/3) N_s m, per heavy mass alpha = 5, 7.5,
      10, 15, 20, inverse-variance mean; compared with k_T + F^2/(N_s kT); rho_I = (k_S^dyn - static)/k_S^dyn;
      sigma(rho_I)^2 = (sigma_kS/k_S)^2 + (sigma_kT/k_S)^2 (the F^2 term's own error neglected, as registered; its size
      is printed). Verdict: agreement within 2 sigma at every cell. Check at all alpha (no verdict):
      k_S^SW = N_s m omega_1^2/K(alpha)^2.  gamma_box = k_S^dyn/k_T against the bulk 1 + Z^2/(Z + eta Z') and against
      1 + F^2/(N_s kT k_T) (A alone); no pass/fail.
  Damping (sec. 1.4 (B), methods sec. 13): per (cell, M) the seed-averaged position ACF (acf_runs.npz, reduce_B.py) fitted
      as C(t) = A e^{-t/tau_T} + B e^{-t/tau_r} cos(omega t) with omega free -- the model, normalisation and bounds of
      paper2_level4_mode_ladder_20261006.fit_free, delete-group jackknife (10 groups) as its jack() -- except that the
      lower bounds of tau_T, tau_r are 2 dt instead of fit_free's absolute 20 sigma-time (set there for periods of
      100-600; here the shortest period is 6). Gamma = 2/tau_r, Gamma^-1 = tau_r/2, P_1 = B [sigma^2].
Exclusion: a trajectory with a health line is not used (sec. 1.7 item 6); conf_worker.sh never moves such a trace into
the cell, so it is absent from red_nu.csv (one: epi8_H_H40_L10, M = 3000, run 5).
Inputs: the summaries fetched from KOA (fetch_confinement.sh): B <cell>/m_<M>/{red_nu.csv, acf_runs.npz, run.log};
A <cell>/x_<pos>/{red_<seed>.csv, run_<seed>.log, summary_<seed>.csv}; the pilot cell for the determinism gate.
Outputs: tables on stdout; paper1_speedofsound/experiments/final/261004_p1_confinement_{shift,identity}.{png,pdf},
261004_p1_confinement_cells.csv, 261004_p1_confinement_damping.csv.
usage (from hspist3/): python3 validation/paper1_confinement_results_261004.py
"""
import contextlib, filecmp, glob, io, math, os, re, sys
import numpy as np, pandas as pd
from scipy.stats import chi2 as CHI2
from scipy.optimize import curve_fit
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import tests_20260913 as T
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)
import paper1_confinement_prereg_20261012 as PR
from paper1_populate_cs_err_20261002 import slope_with_errors, box_delta, TD, X_EDGE

REL_B = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013"
REL_A = "experiments_energy_transfer/paper1_confinement_A_20261013"
CONF = os.path.join(HS, "cluster", "confinement_20261013")
OUT = T.PLOTS
RD, TW, BUILD = T.RDISK, T.WALL_T, "279282b"
POS = (("m2", -2), ("m1", -1), ("0", 0), ("p1", 1), ("p2", 2))
ALPHAS = (0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0)
HEAVY = (5.0, 7.5, 10.0, 15.0, 20.0)
BLUE, BLUE2, RED = "#2a78d6", "#0b3d91", "#e34948"
HEALTH = re.compile(r"EDMD-HEALTH")
P_EXCL = 0.01
# ##CHRIS 2026-10-05 (engine gate G-E6, 261012 sec. 4.4): build-generation guard. This is the registered analysis of the 279282b
# campaign; it must never read data of a later build generation without an explicit flag. Every build of the
# engine-divider-resched generation prints "[EDMD-RESCHED]" into each run log (00ALLINONE.c, edmd_backend_create); 279282b never
# did. method_B and method_A stop on a cell whose logs carry it, unless ALLOW_NEW_BUILD is set: env HD_ALLOW_BUILD_MIX=1, or the
# caller sets the attribute (validation/resched_gate_261005.py does, to compare the two generations). No effect on 279282b data.
ALLOW_NEW_BUILD = os.environ.get("HD_ALLOW_BUILD_MIX") == "1"


def refuse_new_build(logs):
    if ALLOW_NEW_BUILD: return
    for p in logs:
        if os.path.exists(p) and "[EDMD-RESCHED]" in open(edmd_acc_guard.guard(p), errors="ignore").read():
            sys.exit(f"STOP: {p} was written by a post-279282b build ([EDMD-RESCHED] line); the registered 279282b analysis "
                     "does not mix build generations (HD_ALLOW_BUILD_MIX=1 overrides explicitly)")


def cells():
    """The 19 cells from the array lists, geometry from the task files (what ran), cross-checked with the registration."""
    with contextlib.redirect_stdout(io.StringIO()):
        C, _ = PR.main()
    reg = {}
    for c in C:
        reg.setdefault((c["eta_lab"], round(c["H"], 6), round(c["L0"], 6)), c)
    out = []
    for lab, g in (("0.10", "0.10"), ("0.39", "0.39")):
        for cid in open(os.path.join(CONF, f"cells_B_{g}.tsv")).read().split():
            fb = [l.split() for l in open(os.path.join(CONF, f"tasks_B_{cid}.txt"))]
            fa = [l.split() for l in open(os.path.join(CONF, f"tasks_A_{cid}.txt"))]
            L0, H, Ns = float(fb[0][5]), float(fb[0][6]), int(fb[0][7])
            xw = {}
            for f in fa:
                xw.setdefault(os.path.basename(f[1])[2:], float(f[2]))
            seeds = {}
            for f in fa:
                seeds.setdefault(os.path.basename(f[1])[2:], []).append(int(f[3]))
            dL = xw["p1"] - xw["0"]
            r = reg[(lab, round(H, 6), round(L0, 6))]
            scan = cid.split("_")[1]
            out.append(dict(cid=cid, lab=lab, scan=scan, L0=L0, H=H, Ns=Ns, dL=dL, xw=xw, seeds=seeds,
                            Ms=sorted({int(f[2]) for f in fb}), reg=r, n_B=len(fb), n_A=len(fa)))
    return out


# ------------------------------------------------------------------------------------------------ X1 inventory
def inventory(CS):
    print("### X1 -- inventory of the fetched summaries (expected = task file)\n")
    print("| eta | cell | B masses | B trajectories (exp.) | B files, MB | B nu rows with n missing | B wall_x = 200 + 24 L0 px "
          "| A positions x seeds (exp.) | A files, MB | A window min..max | health lines (B run.log / A run logs) "
          "| A build | A geometry (L0, H, 2N_s, t, box, eta, x_wall) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    ok_all = True
    for c in CS:
        bd = os.path.join(HS, REL_B, c["cid"]); ad = os.path.join(HS, REL_A, c["cid"])
        nB = nmiss = hB = 0; bf = bsz = 0; wx_bad = 0; nwx = 0
        for M in c["Ms"]:
            d = os.path.join(bd, f"m_{M}")
            for f in ("red_nu.csv", "acf_runs.npz", "run.log"):
                p = os.path.join(d, f)
                if os.path.exists(p):
                    bf += 1; bsz += os.path.getsize(p)
            r = pd.read_csv(edmd_acc_guard.guard(os.path.join(d, "red_nu.csv"))); nB += len(r); nmiss += int(r["n"].isna().sum())
            lg = open(edmd_acc_guard.guard(os.path.join(d, "run.log")), errors="ignore").read(); hB += len(HEALTH.findall(lg))
            for v in re.findall(r"Initial wall_x = ([\d.]+)", lg):
                nwx += 1; wx_bad += abs(float(v) - (200 + 24 * c["L0"])) > 1e-3
        nA = hA = af = asz = 0; wins = []; builds = set(); geo_bad = 0
        for lab, j in POS:
            d = os.path.join(ad, f"x_{lab}")
            for s in c["seeds"][lab]:
                for f in (f"red_{s}.csv", f"run_{s}.log", f"summary_{s}.csv"):
                    p = os.path.join(d, f)
                    if os.path.exists(p):
                        af += 1; asz += os.path.getsize(p)
                rp = os.path.join(d, f"red_{s}.csv")
                if os.path.exists(rp) and os.path.getsize(rp) > 0:
                    nA += 1; wins.append(float(pd.read_csv(edmd_acc_guard.guard(rp))["window"].iloc[0]))
                hA += len(HEALTH.findall(open(edmd_acc_guard.guard(os.path.join(d, f"run_{s}.log")), errors="ignore").read()))
                sm = pd.read_csv(edmd_acc_guard.guard(os.path.join(d, f"summary_{s}.csv")))
                if len(sm) != 1: geo_bad += 1
                sm = sm.iloc[-1]; builds.add(str(sm["build_git"]))
                eta_rec = c["Ns"] * math.pi * RD ** 2 / (c["H"] * c["L0"])
                # ##CHRIS 2026-10-04: tolerances at the summary's print precision -- L0, height and the wall position are
                # printed to 4 decimals (e.g. '19.7917'), so 1/24-grid values differ by up to 3.3e-5 (first run used 1e-5).
                g = [abs(sm["L0"] - c["L0"]) < 5e-5, abs(float(sm["height"]) - c["H"]) < 5e-5,
                     int(sm["particles_total"]) == 2 * c["Ns"], abs(sm["wall_thickness_sigma"] - TW) < 1e-6,
                     abs(sm["box_width_sigma"] - 2 * c["L0"]) < 1e-5, abs(sm["eta_nominal"] - eta_rec) < 1e-6,
                     abs(float(sm["wall_positions_cli"]) - c["xw"][lab]) < 5e-5]
                geo_bad += not all(g)
        ok = (nB == c["n_B"] - (1 if c["cid"] == "epi8_H_H40_L10" else 0) and nA == c["n_A"] and nmiss == 0
              and hB == 0 and hA == 0 and builds == {BUILD} and geo_bad == 0 and wx_bad == 0
              and max(abs(w - 5000) for w in wins) <= 50)
        ok_all &= ok
        print(f"| {c['lab']} | {c['cid']} | {len(c['Ms'])} | {nB} ({c['n_B']}) | {bf}, {bsz / 1e6:.2f} | {nmiss} | "
              f"{nwx - wx_bad}/{nwx} | 5 x {nA // 5} ({c['n_A']}) | {af}, {asz / 1e6:.2f} | {min(wins):.1f}..{max(wins):.1f} | "
              f"{hB} / {hA} | {','.join(sorted(builds))} | {'all as registered' if geo_bad == 0 else f'**{geo_bad} differ**'} |")
    print(f"\ninventory: {'every cell complete and as registered, except the one excluded B trajectory' if ok_all else '**NOT COMPLETE -- see table**'}")
    return ok_all


def reduction_gate():
    """sec. 1.10 and the reduce_B.py docstring: on the full KOA pilot cell (fetched), reduce_B.py's red_nu.csv must equal the
    canonical cell()'s per-mass nu. reduce_B.py is run on a temporary copy, so nothing is written into the data tree."""
    import shutil, subprocess, tempfile
    from paper1_populate_cs_err_20261002 import cell
    src = os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10")
    print("\n### Reduction gate (sec. 1.10): reduce_B.py vs the canonical cell() on the full KOA pilot cell (pi/8 anchor, 1 run per mass)\n")
    print("| M | nu, reduce_B.py (red_nu.csv) | nu, cell() | abs. difference |\n|---|---|---|---|")
    worst = 0.0
    with tempfile.TemporaryDirectory() as tmp:
        d0 = os.path.join(tmp, "cell")
        shutil.copytree(edmd_acc_guard.guard(src), d0, ignore=shutil.ignore_patterns("_determinism"))
        subprocess.run([sys.executable, os.path.join(CONF, "reduce_B.py"), d0], check=True, capture_output=True, cwd=HS)
        for d in sorted(glob.glob(os.path.join(d0, "m_*")), key=lambda q: int(q.rsplit("_", 1)[1])):
            M = int(d.rsplit("_", 1)[1]); r = pd.read_csv(os.path.join(d, "red_nu.csv"))
            runs = [(int(q.rsplit("_run", 1)[1][:-4]), q, False) for q in sorted(glob.glob(os.path.join(d, "wall_x_positions_L0_*_run*.csv")))]
            cc = cell((0.392699, 10.0, M, runs)); diff = abs(float(r["nu"].mean()) - cc["nu"]); worst = max(worst, diff)
            print(f"| {M} | {r['nu'].mean():.17g} | {cc['nu']:.17g} | {diff:.1e} |")
    print(f"\nreduction gate: max abs. difference {worst:.1e} -> {'PASS' if worst < 1e-12 else 'FAIL'}")
    return worst < 1e-12


def determinism_gate():
    a = os.path.join(HS, REL_A, "epi8_H_H10_L10"); p = os.path.join(HS, REL_A, "pilot_epi8_H_H10_L10")
    print("\n### Determinism gate (sec. 1.12): anchor cell of conf_A_0.39 vs the pilot, red_970[0-3].csv, cmp\n")
    print("| position | seed | anchor bytes | pilot bytes | cmp |\n|---|---|---|---|---|")
    n_ok = 0
    for lab, _ in POS:
        for s in range(9700, 9704):
            fa, fp = os.path.join(a, f"x_{lab}", f"red_{s}.csv"), os.path.join(p, f"x_{lab}", f"red_{s}.csv")
            same = filecmp.cmp(fa, fp, shallow=False); n_ok += same
            print(f"| x_{lab} | {s} | {os.path.getsize(fa)} | {os.path.getsize(fp)} | {'IDENTICAL' if same else '**DIFFERENT**'} |")
    print(f"\ndeterminism gate: {n_ok}/20 IDENTICAL -> {'PASS' if n_ok == 20 else 'FAIL'}")
    return n_ok == 20


# ------------------------------------------------------------------------------------------------ method B
def method_B(c):
    rows = []
    refuse_new_build([os.path.join(HS, REL_B, c["cid"], f"m_{M}", "run.log") for M in c["Ms"]])
    for M in c["Ms"]:
        d = pd.read_csv(edmd_acc_guard.guard(os.path.join(HS, REL_B, c["cid"], f"m_{M}", "red_nu.csv")))
        nu = d["nu"].to_numpy(float); n = len(nu)
        al = M / (2.0 * c["Ns"])
        rows.append(dict(M=M, alpha=al, K=T.k_root(al), nu=nu.mean(), sd=nu.std(ddof=1), n=n, se=nu.std(ddof=1) / math.sqrt(n),
                         nus=nu, dt=float(d["dt"].iloc[0]), nu_pred=float(d["nu_pred"].iloc[0])))
    delta = box_delta(c["L0"]); LeT = T.l_eff(c["L0"]) - delta / 2
    x = np.array([r["K"] / (2 * math.pi * LeT) for r in rows]); y = np.array([r["nu"] for r in rows])
    sy = np.array([r["se"] for r in rows])
    s, err, errs, chi2r = slope_with_errors(x, y, sy)
    eta_rec = c["Ns"] * math.pi * RD ** 2 / (c["H"] * c["L0"]); etaT = eta_rec * c["L0"] / (c["L0"] - delta / 2)
    KR = T.kr_cs(etaT)
    Z, dZ = PR.Z(etaT), PR.dZ(etaT); q = (Z + etaT * dZ + Z * Z) / Z
    c.update(B=rows, delta=delta, LeT=LeT, x=x, eta_rec=eta_rec, eta=etaT, cs=s, cs_err=err, cs_errs=errs, chi2r=chi2r, KR=KR,
             D=s / KR - 1, sD=errs / KR, DC=(q + 1) * (q + 2) / (16 * c["Ns"] * q * Z), Z=Z, gbulk=1 + Z * Z / (Z + etaT * dZ),
             nB=sum(r["n"] for r in rows))


SHAPES = {"A": lambda c: 2 / c["H"] + 2 / c["L0"], "B": lambda c: 1 / c["H"], "C": lambda c: 1 / c["Ns"]}


def fit_one(D, S, f):
    w = 1 / S ** 2; a = float((w * f * D).sum() / (w * f * f).sum()); sa = float(1 / math.sqrt((w * f * f).sum()))
    ch = float((w * (D - a * f) ** 2).sum()); dof = len(D) - 1
    return dict(a=a, sa=sa, chi2=ch, dof=dof, p=float(CHI2.sf(ch, dof)))


def fit_two(D, S, f1, f2):
    w = 1 / S ** 2; X = np.vstack([f1, f2]).T; W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X); b = cov @ X.T @ W @ D
    ch = float((w * (D - X @ b) ** 2).sum()); dof = len(D) - 2
    return dict(b=b, sb=np.sqrt(np.diag(cov)), chi2=ch, dof=dof, p=float(CHI2.sf(ch, dof)))


def confinement_verdict(CS):
    res = {}
    print("\n### Confinement: per-cell shift from method B (canonical estimator, campaign cells only)\n")
    print("| eta | cell | scan | H | L_0 | N_s | eta_true | c_s | c_s_err | chi2_red (9 masses) | c_s_err_scaled | c_s^KR(eta_true) "
          "| Delta = c_s/c_s^KR - 1 [%] | sigma [%] | shape A: 2/H + 2/L_0 | shape B: 1/H | shape C: 1/N_s | Delta_C (fixed) [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in CS:
        print(f"| {c['lab']} | {c['cid']} | {c['scan']} | {c['H']:g} | {c['L0']:g} | {c['Ns']} | {c['eta']:.6f} | {c['cs']:.5f} | "
              f"{c['cs_err']:.5f} | {c['chi2r']:.2f} | {c['cs_errs']:.5f} | {c['KR']:.5f} | {100 * c['D']:+.3f} | {100 * c['sD']:.3f} | "
              f"{SHAPES['A'](c):.4f} | {SHAPES['B'](c):.4f} | {SHAPES['C'](c):.4f} | {100 * c['DC']:+.3f} |")
    for lab in ("0.10", "0.39"):
        cc = [c for c in CS if c["lab"] == lab]
        D = np.array([c["D"] for c in cc]); S = np.array([c["sD"] for c in cc])
        fits = {h: fit_one(D, S, np.array([SHAPES[h](c) for c in cc])) for h in SHAPES}
        DCv = np.array([c["DC"] for c in cc]); chC = float((((D - DCv) / S) ** 2).sum())
        fits["C_fixed"] = dict(a=1.0, sa=0.0, chi2=chC, dof=len(D), p=float(CHI2.sf(chC, len(D))))
        print(f"\n#### eta {lab}: one-amplitude fits over its {len(cc)} cells (sec. 1.5 decision rule; excluded if p < {P_EXCL})\n")
        print("| hypothesis | shape | amplitude | chi2 | dof | p(chi2) | excluded? |\n|---|---|---|---|---|---|---|")
        names = {"A": "a (2/H + 2/L_0)", "B": "b/H", "C": "c'/N_s", "C_fixed": "Delta_C, no free parameter"}
        for h, f in fits.items():
            amp = "fixed" if h == "C_fixed" else f"{f['a']:.5f} +- {f['sa']:.5f}"
            print(f"| {h} | {names[h]} | {amp} | {f['chi2']:.2f} | {f['dof']} | {f['p']:.3g} | {'YES' if f['p'] < P_EXCL else 'no'} |")
        surv = [h for h in ("A", "B", "C") if fits[h]["p"] >= P_EXCL]
        if len(surv) == 1:
            verdict = f"{surv[0]} (the only survivor)"
        elif len(surv) > 1:
            best = min(surv, key=lambda h: fits[h]["chi2"])
            verdict = "not separated: " + ", ".join(surv) + " survive; Delta chi2 to the best (" + best + "): " + \
                      ", ".join(f"{h} {fits[h]['chi2'] - fits[best]['chi2']:+.2f}" for h in surv if h != best)
        else:
            verdict = "none survives"
        print(f"\n**Verdict, eta {lab}: {verdict}.**")
        two = {}
        if not surv:
            fA = np.array([SHAPES["A"](c) for c in cc]); fB = np.array([SHAPES["B"](c) for c in cc]); fC = np.array([SHAPES["C"](c) for c in cc])
            two = {"b/H + c'/N_s": fit_two(D, S, fB, fC), "a(2/H + 2/L_0) + c'/N_s": fit_two(D, S, fA, fC)}
            print("\nExploratory two-term forms (registered as exploratory only, not a verdict):\n")
            print("| form | first amplitude | c' | chi2 | dof | p(chi2) |\n|---|---|---|---|---|---|")
            for k, f in two.items():
                print(f"| {k} | {f['b'][0]:.5f} +- {f['sb'][0]:.5f} | {f['b'][1]:.5f} +- {f['sb'][1]:.5f} | {f['chi2']:.2f} | {f['dof']} | {f['p']:.3g} |")
        res[lab] = dict(fits=fits, surv=surv, verdict=verdict, two=two)
    return res


# ------------------------------------------------------------------------------------------------ method A
def method_A(c):
    ad = os.path.join(HS, REL_A, c["cid"]); st = {}
    refuse_new_build([os.path.join(ad, f"x_{lab}", f"run_{s}.log") for lab, _ in POS for s in c["seeds"][lab]])
    for lab, j in POS:
        R = pd.concat([pd.read_csv(edmd_acc_guard.guard(os.path.join(ad, f"x_{lab}", f"red_{s}.csv"))) for s in c["seeds"][lab]], ignore_index=True)
        n = len(R)
        st[j] = dict(n=n, FL=R["F_L"].mean(), FR=R["F_R"].mean(), sFL=R["F_L"].std(ddof=1) / math.sqrt(n),
                     sFR=R["F_R"].std(ddof=1) / math.sqrt(n), T=float(((R["T_L"] + R["T_R"]) / 2).mean()))
    F = {j: 0.5 * (st[j]["FL"] + st[-j]["FR"]) for j in (-2, -1, 0, 1, 2)}
    sF = {j: 0.5 * math.sqrt(st[j]["sFL"] ** 2 + st[-j]["sFR"] ** 2) for j in (-2, -1, 0, 1, 2)}
    dL = c["dL"]
    kT_ = -(F[-2] - 8 * F[-1] + 8 * F[1] - F[2]) / (12 * dL)
    s_kT = math.sqrt(sF[-2] ** 2 + 64 * sF[-1] ** 2 + 64 * sF[1] ** 2 + sF[2] ** 2) / (12 * dL)
    temp = st[0]["T"]
    z = {"x=0": (st[0]["FL"] - st[0]["FR"]) / math.hypot(st[0]["sFL"], st[0]["sFR"])}
    for j in (-2, -1, 1, 2):
        z[f"x={j:+d}dL"] = (st[j]["FL"] - st[-j]["FR"]) / math.hypot(st[j]["sFL"], st[-j]["sFR"])
    F2 = F[0] ** 2 / (c["Ns"] * temp)
    c.update(A=st, F=F, sF=sF, F0=F[0], sF0=sF[0], kT=kT_, s_kT=s_kT, temp=temp, F2term=F2, s_F2=2 * F[0] * sF[0] / (c["Ns"] * temp),
             static=kT_ + F2, sym=z, nA=sum(st[j]["n"] for j in st))


def identity(c):
    heavy = [r for r in c["B"] if any(abs(r["alpha"] - a) < 1e-9 for a in HEAVY)]
    k, s = [], []
    for r in heavy:
        om = 2 * math.pi * r["nu"]; Mh = r["M"] + 2.0 * c["Ns"] / 3.0
        k.append(Mh * om * om / 2.0); s.append(Mh * om * om / 2.0 * 2 * r["se"] / r["nu"])
    k, s = np.array(k), np.array(s); w = 1 / s ** 2
    kS = float((w * k).sum() / w.sum()); s_kS = float(1 / math.sqrt(w.sum()))
    chi2h = float((w * (k - kS) ** 2).sum()) / (len(k) - 1)
    rho = (kS - c["static"]) / kS; s_rho = math.sqrt((s_kS / kS) ** 2 + (c["s_kT"] / kS) ** 2)
    sw = {r["alpha"]: c["Ns"] * (2 * math.pi * r["nu"]) ** 2 / r["K"] ** 2 for r in c["B"]}
    g = kS / c["kT"]; s_g = g * math.sqrt((s_kS / kS) ** 2 + (c["s_kT"] / c["kT"]) ** 2)
    c.update(kS=kS, s_kS=s_kS, chi2h=chi2h, kS_each=k, kS_sig=s, rho=rho, s_rho=s_rho, zrho=rho / s_rho, sw=sw, gbox=g, s_gbox=s_g,
             gA=1 + c["F2term"] / c["kT"])


# ------------------------------------------------------------------------------------------------ damping
def _g(t, A, tT, B, tr, om):
    return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(om * t)


def _fit_acf(C, dt, nu):
    lag = np.arange(len(C)) * dt; c = C / C[0]; om0 = 2 * math.pi * nu; P = 1 / nu
    p0 = [0.3, 2 * P, 0.7, 5 * P, om0]; lo = [0, 2 * dt, 0, 2 * dt, 0.2 * om0]; hi = [1.5, 1e7, 1.5, 1e7, 5.0 * om0]
    p, _ = curve_fit(_g, lag, c, p0=p0, bounds=(lo, hi), maxfev=200000)
    return p


def damping(c):
    out = []
    for r in c["B"]:
        z = np.load(edmd_acc_guard.guard(os.path.join(HS, REL_B, c["cid"], f"m_{r['M']}", "acf_runs.npz")))
        runs = [z[k].astype(float) for k in sorted(z.files, key=lambda s: int(s[3:]))]
        L = min(len(a) for a in runs); A = np.array([a[:L] for a in runs]); Cm = A.mean(0)
        try:
            p = _fit_acf(Cm, r["dt"], r["nu"])
            jk = []
            for gi in range(10):
                keep = [i for i in range(len(A)) if i % 10 != gi]
                try:
                    jk.append(_fit_acf(A[keep].mean(0), r["dt"], r["nu"]))
                except Exception:
                    pass
            if len(jk) >= 3:
                J = np.array(jk); G = len(jk); e = np.sqrt((G - 1) / G * ((J - J.mean(0)) ** 2).sum(0))
            else:
                e = np.full(5, np.nan)
            out.append(dict(M=r["M"], alpha=r["alpha"], n=len(A), tau_r=p[3], s_tau_r=e[3], Ginv=p[3] / 2, s_Ginv=e[3] / 2,
                            Gamma=2 / p[3], P1=p[2] * Cm[0], s_P1=e[2] * Cm[0], tau_T=p[1], om_fit=p[4], nu=r["nu"],
                            Q=math.pi * r["nu"] * p[3], window=L * r["dt"], ok=True))
        except Exception as ex:
            out.append(dict(M=r["M"], alpha=r["alpha"], n=len(A), ok=False, err=str(ex)[:60]))
    c["damp"] = out


# ------------------------------------------------------------------------------------------------ X2 bound
def exclusion_bound(c):
    rows = c["B"]; i = [k for k, r in enumerate(rows) if r["M"] == 3000][0]; r = rows[i]
    sxx = float((c["x"] ** 2).sum())
    dev_obs = float(np.max(np.abs(r["nus"] - r["nu"]))); dev3 = 3 * r["sd"]
    print("\n### X2 -- the excluded trajectory: how far could it have moved the cell?\n")
    print(f"cell {c['cid']}, M = 3000 (alpha = {r['alpha']:g}): {r['n']} of 25 runs used; nu = {r['nu']:.7f}, sd over seeds = {r['sd']:.7f}")
    for lab, dev in (("largest deviation among the 24 used runs", dev_obs), ("3 sd", dev3)):
        dnu = dev / 25.0; dcs = c["x"][i] * dnu / sxx
        print(f"- a 25th run {lab} ({dev:.7f}) away from the mean moves that mass's nu by {dnu:.2e} and c_s by {dcs:.2e} = "
              f"{dcs / c['cs_errs']:.3f} of the cell's c_s_err_scaled ({c['cs_errs']:.5f})")
    k = [j for j, rr in enumerate([q for q in rows if any(abs(q['alpha'] - a) < 1e-9 for a in HEAVY)]) if rr["M"] == 3000][0]
    dk = 2 * 3 * r["sd"] / 25 / r["nu"] * c["kS_each"][k]; w = 1 / c["kS_sig"] ** 2; dcell = w[k] / w.sum() * dk
    print(f"- the same mass is one of the five heavy masses of the identity: its k_S^dyn = {c['kS_each'][k]:.5f}; "
          f"a 3-sd 25th run changes it by {dk:.2e}; with its inverse-variance weight {w[k] / w.sum():.3f} the cell's combined "
          f"k_S^dyn moves by {dcell:.2e} = {dcell / c['s_kS']:.2f} of its sigma ({c['s_kS']:.2e}) and rho_I by {100 * dcell / c['kS']:.3f} %")


# ------------------------------------------------------------------------------------------------ figures
def fig_shift(CS, res):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(15, 9.2), sharey="row")
    sty = {"A": ("#eb6834", "--"), "B": ("#1baf7a", "-."), "C": ("#8f4fd1", ":")}
    for row, lab in enumerate(("0.10", "0.39")):
        cc = [c for c in CS if c["lab"] == lab]; F = res[lab]["fits"]
        anc = [c for c in cc if c["scan"] == "H" and abs(c["H"] - 10) < 1e-9][0]
        eta, L0a, Ha = anc["eta"], anc["L0"], anc["H"]
        Z, dZ = PR.Z(eta), PR.dZ(eta); q = (Z + eta * dZ + Z * Z) / Z
        Ns = lambda H, L: 4 * eta * H * L / math.pi
        pan = [("H-scan at L_0 = %g:  Delta vs 1/H" % L0a, [c for c in cc if c["scan"] == "H"], lambda c: 1 / c["H"],
                np.linspace(1 / 45, 1 / 4.5, 120), lambda u: (1 / u, L0a), "1/H  [1/sigma]"),
               ("L-scan at H = 10:  Delta vs 1/L_0", [c for c in cc if c["scan"] == "L"] + [anc], lambda c: 1 / c["L0"],
                np.linspace(1 / (8.6 * L0a / 4), 1 / (0.45 * L0a), 120), lambda u: (Ha, 1 / u), "1/L_0  [1/sigma]"),
               ("aspect at N_s = 50, fixed area:  Delta vs L_0/H", [c for c in cc if c["scan"] == "aspect"] + ([anc] if lab == "0.39" else []),
                lambda c: c["L0"] / c["H"], np.linspace(0.8, 8.6, 120),
                lambda u: (math.sqrt(anc["H"] * anc["L0"] / u), math.sqrt(anc["H"] * anc["L0"] * u)), "L_0/H")]
        for col, (title, pts, xf, grid, geo, xl) in enumerate(pan):
            ax = axs[row, col]
            ax.axhline(0, color=RED, lw=2.0, label="KR 2006 (rho_max 0.90 fit): Delta = 0" if col == 0 else None)
            for h in ("A", "B", "C"):
                col_, ls = sty[h]; f = F[h]
                ys = []
                for u in grid:
                    H, L = geo(u); cdict = dict(H=H, L0=L, Ns=(50 if col == 2 else Ns(H, L)))
                    ys.append(100 * f["a"] * SHAPES[h](cdict))
                ax.plot(grid, ys, color=col_, ls=ls, lw=1.6,
                        label=f"{h}: {'2/H+2/L_0' if h == 'A' else ('1/H' if h == 'B' else '1/N_s')}, fit "
                              f"(chi2/dof {f['chi2']:.1f}/{f['dof']}, p = {f['p']:.2g})" if col == 0 else None)
            ysC = [100 * (q + 1) * (q + 2) / (16 * (50 if col == 2 else Ns(*geo(u))) * q * Z) for u in grid]
            ax.plot(grid, ysC, color=sty["C"][0], lw=0.9, alpha=0.7,
                    label=f"C, fixed Delta_C (chi2/dof {F['C_fixed']['chi2']:.1f}/{F['C_fixed']['dof']})" if col == 0 else None)
            X = np.array([xf(c) for c in pts]); Y = np.array([100 * c["D"] for c in pts]); E = np.array([100 * c["sD"] for c in pts])
            ax.errorbar(X, Y, yerr=E, fmt="o", color=BLUE, ms=6, capsize=3, lw=1.3, zorder=5,
                        label="campaign, method B (25 seeds x 9 masses per cell)" if col == 0 else None)
            for c, xx, yy in zip(pts, X, Y):
                ax.annotate(f"N_s {c['Ns']}", (xx, yy), textcoords="offset points", xytext=(5, 6), fontsize=7.5, color="0.35")
            ax.set_title(("eta = 0.100   " if lab == "0.10" else "eta = pi/8   ") + title, fontsize=10)
            ax.set_xlabel(xl); ax.grid(True, ls=":", alpha=0.6)
            if col == 0: ax.set_ylabel("Delta = c_s / c_s^KR(eta_true) - 1  [%]")
        axs[row, 0].legend(fontsize=7.6, loc="upper left", framealpha=0.95)
    fig.suptitle("Confinement shift of the divider-mode sound speed (pre-registered campaign, 261012 sec. 1; "
                 "error bars = c_s_err_scaled; amplitudes fitted over all cells of each eta)", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"261004_p1_confinement_shift.{ext}"), dpi=200)
    print("\nfigure -> 261004_p1_confinement_shift.png/.pdf")


def fig_identity(CS):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6.2), gridspec_kw=dict(width_ratios=[1, 1.25]))
    for lab, col, mk, nm in (("0.10", BLUE, "o", "eta = 0.100"), ("0.39", BLUE2, "s", "eta = pi/8")):
        cc = [c for c in CS if c["lab"] == lab]
        a1.errorbar([c["static"] for c in cc], [c["kS"] for c in cc], xerr=[c["s_kT"] for c in cc], yerr=[c["s_kS"] for c in cc],
                    fmt=mk, color=col, ms=6, capsize=2.5, label=nm + " (method B vs method A)")
    lo = min(c["static"] for c in CS) / 1.6; hi = max(c["static"] for c in CS) * 1.6
    a1.plot([lo, hi], [lo, hi], color="k", ls="--", lw=1.2, label="identity  k_S = k_T + F^2/(N_s kT)")
    a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlim(lo, hi); a1.set_ylim(lo, hi)
    a1.set_xlabel("static, method A:  k_T + F^2/(N_s kT)  [kT/sigma^2]"); a1.set_ylabel("dynamic, method B:  k_S^dyn = M_hat omega_1^2 / 2")
    a1.grid(True, which="both", ls=":", alpha=0.5); a1.legend(fontsize=8.5, loc="upper left")
    a1.set_title("Identity, length-free form (amendment C1)", fontsize=10.5)
    order = sorted(CS, key=lambda c: (c["lab"], c["scan"], c["Ns"], c["L0"]))
    xs = np.arange(len(order))
    for i, c in enumerate(order):
        col = BLUE if c["lab"] == "0.10" else BLUE2; mk = "o" if c["lab"] == "0.10" else "s"
        a2.errorbar(i, 100 * c["rho"], yerr=200 * c["s_rho"], fmt=mk, color=col, ms=6, capsize=0, lw=0.8, alpha=0.45)
        a2.errorbar(i, 100 * c["rho"], yerr=100 * c["s_rho"], fmt=mk, color=col, ms=6, capsize=3, lw=1.4)
        a2.plot(i, 200 * c["DC"], marker="_", color="#8f4fd1", ms=14, mew=2)
    a2.axhline(0, color="k", ls="--", lw=1.0)
    a2.plot([], [], marker="_", color="#8f4fd1", ls="none", ms=14, mew=2, label="hypothesis C: rho_I = 2 Delta_C")
    a2.errorbar([], [], yerr=[], fmt="o", color=BLUE, capsize=3, label="rho_I, 1 sigma (thin: 2 sigma)")
    a2.set_xticks(xs); a2.set_xticklabels([c["cid"].replace("e0p10_", "0.10 ").replace("epi8_", "pi/8 ") for c in order],
                                          rotation=70, ha="right", fontsize=7)
    a2.set_ylabel("rho_I = (k_S^dyn - k_T - F^2/(N_s kT)) / k_S^dyn  [%]"); a2.grid(True, ls=":", alpha=0.5)
    a2.legend(fontsize=8.5, loc="upper left"); a2.set_title("Residual per cell; verdict rule: |rho_I| <= 2 sigma at every cell", fontsize=10.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"261004_p1_identity.{ext}"), dpi=200)
    print("figure -> 261004_p1_identity.png/.pdf")


# ------------------------------------------------------------------------------------------------ main
def main():
    CS = cells()
    print(f"cells: {len(CS)} ({sum(c['lab'] == '0.10' for c in CS)} at eta 0.10, {sum(c['lab'] == '0.39' for c in CS)} at pi/8); "
          f"build {BUILD}; registration 261012 sec. 1 + C1-C3\n")
    inv_ok = inventory(CS)
    red_ok = reduction_gate()
    det_ok = determinism_gate()
    for c in CS:
        method_B(c); method_A(c); identity(c); damping(c)
    # ##CHRIS 2026-10-04: the task files carry L_0, H and the wall positions to 6 decimals, so compare at 1e-6 (first run: 1e-9)
    reg_bad = [c["cid"] for c in CS if abs(c["eta_rec"] - c["reg"]["eta"]) > 1e-6 or abs(c["dL"] - c["reg"]["dL"]) > 1e-6]
    print(f"max |eta_rec - eta_reg| = {max(abs(c['eta_rec'] - c['reg']['eta']) for c in CS):.1e}, "
          f"max |dL - dL_reg| = {max(abs(c['dL'] - c['reg']['dL']) for c in CS):.1e}")
    print(f"\ngeometry vs registration (eta, dL from paper1_confinement_prereg_20261012): {'all equal' if not reg_bad else reg_bad}; "
          f"box shortfall delta: max {max(c['delta'] for c in CS):.2e} sigma (grid-exact boxes)")
    res = confinement_verdict(CS)

    print("\n### Method A per cell: F at L_0, k_T (5-point stencil, mirror faces averaged), kT, symmetry checks\n")
    print("| eta | cell | dL | seeds/position | F(L_0) | sigma_F | kT (x = 0 seeds) | k_T | sigma(k_T) | sigma(k_T)/k_T [%] "
          "| symmetry z: x=0; x=+1,+2,-1,-2 dL | |z| > 2 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    nz = 0
    for c in CS:
        zz = c["sym"]; big = sum(abs(v) > 2 for v in zz.values()); nz += big
        print(f"| {c['lab']} | {c['cid']} | {c['dL']:.4f} | {c['A'][0]['n']} | {c['F0']:.5f} | {c['sF0']:.5f} | {c['temp']:.6f} | "
              f"{c['kT']:.5f} | {c['s_kT']:.5f} | {100 * c['s_kT'] / c['kT']:.2f} | "
              f"{zz['x=0']:+.2f}; {zz['x=+1dL']:+.2f}, {zz['x=+2dL']:+.2f}, {zz['x=-1dL']:+.2f}, {zz['x=-2dL']:+.2f} | {big} |")
    ntest = 5 * len(CS)
    print(f"\nsymmetry checks beyond 2 sigma: {nz} of {ntest} (registered: each within 2 sigma; "
          f"expected by chance if all hold: {ntest * 0.0455:.1f})")

    print("\n### Identity, length-free form (amendment C1): k_S^dyn (heavy masses) vs k_T + F^2/(N_s kT)\n")
    print("| eta | cell | N_s | k_S^dyn | sigma | chi2_red of the 5 heavy masses | k_T | F^2/(N_s kT) | its own sigma (neglected) "
          "| static = k_T + F^2/(N_s kT) | rho_I [%] | sigma(rho_I) [%] | rho_I/sigma | within 2 sigma? | 2 Delta_C (hyp. C) [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fails = []
    for c in CS:
        okc = abs(c["zrho"]) <= 2; fails += [] if okc else [c["cid"]]
        print(f"| {c['lab']} | {c['cid']} | {c['Ns']} | {c['kS']:.6f} | {c['s_kS']:.6f} | {c['chi2h']:.2f} | {c['kT']:.6f} | "
              f"{c['F2term']:.6f} | {c['s_F2']:.2e} | {c['static']:.6f} | {100 * c['rho']:+.3f} | {100 * c['s_rho']:.3f} | "
              f"{c['zrho']:+.2f} | {'yes' if okc else '**NO**'} | {200 * c['DC']:+.3f} |")
    chi_all = sum(c["zrho"] ** 2 for c in CS)
    print(f"\n**Identity verdict (C1 rule: agreement within 2 sigma at every cell): "
          f"{'PASS' if not fails else 'FAIL'}** -- {len(CS) - len(fails)} of {len(CS)} cells within 2 sigma"
          f"{'' if not fails else '; outside: ' + ', '.join(fails)}.")
    print(f"(information, not the registered rule: sum of (rho_I/sigma)^2 over the {len(CS)} cells = {chi_all:.1f}, "
          f"p = {CHI2.sf(chi_all, len(CS)):.3g}; P(all {len(CS)} within 2 sigma | identity exact) = {0.9545 ** len(CS):.2f})")

    print("\n### Standing-wave check at all alpha (C1, no verdict): k_S^SW(alpha)/k_S^dyn, k_S^SW = N_s m omega_1^2/K(alpha)^2\n")
    print("| eta | cell | " + " | ".join(f"alpha {a:g}" for a in ALPHAS) + " |\n|---|---|" + "---|" * len(ALPHAS))
    for c in CS:
        print(f"| {c['lab']} | {c['cid']} | " + " | ".join(f"{c['sw'][a] / c['kS']:.4f}" for a in ALPHAS) + " |")

    print("\n### gamma_box = k_S^dyn / k_T (no pass/fail; sec. 1.5)\n")
    print("| eta | cell | scan | H | L_0 | N_s | gamma_box | sigma | bulk 1 + Z^2/(Z + eta Z') | 1 + F^2/(N_s kT k_T) (A alone) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for c in CS:
        print(f"| {c['lab']} | {c['cid']} | {c['scan']} | {c['H']:g} | {c['L0']:g} | {c['Ns']} | {c['gbox']:.4f} | {c['s_gbox']:.4f} | "
              f"{c['gbulk']:.5f} | {c['gA']:.4f} |")
    for lab in ("0.10", "0.39"):
        cc = [c for c in CS if c["lab"] == lab]; g = np.array([c["gbox"] for c in cc]); s = np.array([c["s_gbox"] for c in cc])
        w = 1 / s ** 2; m = float((w * g).sum() / w.sum()); sm = float(1 / math.sqrt(w.sum()))
        ch = float((w * (g - m) ** 2).sum()); gA = np.array([c["gA"] for c in cc])
        print(f"\neta {lab}: gamma_box inverse-variance mean {m:.4f} +- {sm:.4f} (chi2 {ch:.1f} / {len(g) - 1} dof), "
              f"range {g.min():.4f} .. {g.max():.4f}; bulk {cc[0]['gbulk']:.5f}; A-alone mean {gA.mean():.4f}")

    print("\n### Damping per (cell, M): Gamma^-1 = tau_r/2 [sigma-time] +- jackknife (methods sec. 13 model, ACF to 20 periods)\n")
    print("| eta | cell | " + " | ".join(f"alpha {a:g}" for a in ALPHAS) + " |\n|---|---|" + "---|" * len(ALPHAS))
    nfail = 0
    for c in CS:
        cells_ = []
        for d in c["damp"]:
            if d["ok"]:
                cells_.append(f"{d['Ginv']:.3g} +- {d['s_Ginv']:.2g}")
            else:
                cells_.append("fit failed"); nfail += 1
        print(f"| {c['lab']} | {c['cid']} | " + " | ".join(cells_) + " |")
    print(f"\nACF fits: {sum(len(c['damp']) for c in CS) - nfail} of {sum(len(c['damp']) for c in CS)} converged")

    exclusion_bound([c for c in CS if c["cid"] == "epi8_H_H40_L10"][0])

    print("\n### X3 (a) -- per-cell table\n")
    print("| eta | cell | H | L_0 | N_s | eta_true | c_s^B +- err (scaled) | k_S^dyn +- | k_T +- | F(L_0) | k_T + F^2/(N_s kT) "
          "| gamma_box = k_S/k_T | Gamma^-1 at alpha = 5 (M = 10 N_s) +- |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for c in CS:
        d5 = [d for d in c["damp"] if abs(d["alpha"] - 5) < 1e-9][0]
        print(f"| {c['lab']} | {c['cid']} | {c['H']:g} | {c['L0']:g} | {c['Ns']} | {c['eta']:.6f} | {c['cs']:.4f} +- {c['cs_errs']:.4f} | "
              f"{c['kS']:.5g} +- {c['s_kS']:.2g} | {c['kT']:.5g} +- {c['s_kT']:.2g} | {c['F0']:.5g} | {c['static']:.5g} | "
              f"{c['gbox']:.3f} +- {c['s_gbox']:.3f} | {d5.get('Ginv', float('nan')):.4g} +- {d5.get('s_Ginv', float('nan')):.2g} |")
        rows.append(dict(eta_lab=c["lab"], cell=c["cid"], scan=c["scan"], H=c["H"], L0=c["L0"], N_s=c["Ns"], eta_rec=c["eta_rec"],
                         delta_sigma=c["delta"], eta_true=c["eta"], L_eff_true=c["LeT"], c_s=c["cs"], c_s_err=c["cs_err"],
                         c_s_err_scaled=c["cs_errs"], chi2_red=c["chi2r"], KR=c["KR"], Delta=c["D"], sigma_Delta=c["sD"],
                         Delta_C=c["DC"], trajectories_B=c["nB"], trajectories_A=c["nA"], F_L0=c["F0"], sigma_F_L0=c["sF0"],
                         kT=c["temp"], k_T=c["kT"], sigma_k_T=c["s_kT"], F2_over_NkT=c["F2term"], static=c["static"],
                         k_S_dyn=c["kS"], sigma_k_S_dyn=c["s_kS"], chi2_red_heavy=c["chi2h"], rho_I=c["rho"],
                         sigma_rho_I=c["s_rho"], gamma_box=c["gbox"], sigma_gamma_box=c["s_gbox"], gamma_bulk=c["gbulk"],
                         gamma_A_alone=c["gA"]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "261004_p1_confinement_cells.csv"), index=False)
    pd.DataFrame([dict(cell=c["cid"], eta_lab=c["lab"], **{k: v for k, v in d.items()}) for c in CS for d in c["damp"]]).to_csv(
        os.path.join(OUT, "261004_p1_confinement_damping.csv"), index=False)
    print("\ntables -> 261004_p1_confinement_cells.csv, 261004_p1_confinement_damping.csv")
    fig_shift(CS, res); fig_identity(CS)
    print(f"\ngates: inventory {'PASS' if inv_ok else 'FAIL'}, reduction {'PASS' if red_ok else 'FAIL'}, "
          f"determinism {'PASS' if det_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
