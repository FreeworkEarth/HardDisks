#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.4, decision 3 item 1; EXPLORATORY, POST HOC, no verdict rule): the within-mass aging check
of sec. 4.6 item 4b. Existing data only (the recorded traces and per-run psi6 summaries); nothing is run.

CELLS: those of sec. 4.6 item 4b -- every cell paper1_window_explore_261008.py analyses with eta_true in [0.69, 0.725] whose runs
wrote per-run psi6 summaries (A1 v2, the transition run, the three ladders, route B). The window is eta_true in [0.6995, 0.7175]
(sec. 4.6). Trajectories: as sec. 4.6 (bad runs and failed trace checks excluded; the canonical argmax frequency nu over TD
periods, TD = 200 where the record has 200 planned periods, else its planned periods rounded down).

PER MASS: Spearman(nu, psi6 run mean) over the seeds of one mass (n = 25; 10 in the ladders). Within one mass the
per-trajectory sound-speed estimate c = nu / x_M is nu times a constant, so it has the same ranks. COMBINED per cell: Fisher's z
of the per-mass Spearman values with the Spearman variance 1.06 / (n - 3) (Fieller, Hartley and Pearson, Biometrika 44, 470
(1957)), weights (n - 3) / 1.06, a 95 % interval, and Cochran's Q over the masses with its chi-square p (do the masses agree?).
POOLED (sec. 4.6 item 4b): Spearman of nu/nu_M - 1 against psi6 - <psi6>_M over all trajectories of the cell, recomputed here
with sec. 4.6's code path and compared with its recorded output (exploratory_261008_window/261008_window_explore_output.txt):
the gate requires equality to the printed digits.

SPLIT (first half against second half of each trajectory): possible from the existing data, because the frequency is the
argmax of the spectrum of the recorded divider position; psi6 has no time series (only hold, end and run mean per run). Per
trajectory nu1, nu2 from the first and second TD/2 periods of the same record and the same estimator; dnu = nu2/nu1 - 1. The
argmax resolution of one half is 2/TD (1 % at TD = 200, 8 % at TD = 25), printed per cell. Per cell: the mean of dnu with its
SE (trajectories), the fraction with dnu = 0, and the within-mass Spearman of dnu against dpsi6 = psi6(end) - psi6(hold),
combined over the masses as above.
usage (from hspist3/):  python3 validation/paper1_window_aging_261009.py [--workers 10]
"""
import math, os, re, sys
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.stats import spearmanr, chi2, norm
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import tests_20260913 as T
import paper1_window_explore_261008 as W
import edmd_acc_guard
from paper1_populate_cs_err_20261002 import box_delta

REC = os.path.join(W.OUT, "261008_window_explore_output.txt")
VAR_SPEARMAN = 1.06            # Fieller, Hartley and Pearson (1957): var(atanh r_s) ~ 1.06 / (n - 3)


def traj2(p):
    """sec. 4.6's per-trajectory estimator (W.traj: argmax over TD periods) and the same estimator on each half of that record."""
    try:
        h = W.header(p); per_t = float(h["Planned_Duration"]) * float(h["Predicted_Frequency"])
        TD = 200 if per_t >= 199.999 else int(per_t + 1e-9)
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD)
        if n is None: return None
        def argmax_nu(xx, td):
            P, df = T._spectrum(xx, dt); k = int(round(td / W.X_EDGE)); return (k + int(np.argmax(P[k:]))) * df
        nu = argmax_nu(x[:n], TD)
        h2 = n // 2
        return dict(nu=nu, nu1=argmax_nu(x[:h2], TD / 2.0), nu2=argmax_nu(x[h2:2 * h2], TD / 2.0), TD=TD)
    except edmd_acc_guard.AcceleratedRunError:
        return "refused"
    except Exception:
        return None


def fisher(rs, ns):
    """combined Spearman over masses: Fisher z with var 1.06/(n-3); 95 % interval; Cochran's Q and its p"""
    rs, ns = np.array(rs, float), np.array(ns, float)
    ok = np.isfinite(rs) & (ns > 3) & (np.abs(rs) < 1)
    if ok.sum() < 2: return (float("nan"),) * 5
    z = np.arctanh(rs[ok]); w = (ns[ok] - 3) / VAR_SPEARMAN
    zc = float((w * z).sum() / w.sum()); se = 1 / math.sqrt(w.sum())
    Q = float((w * (z - zc) ** 2).sum()); dof = int(ok.sum()) - 1
    return math.tanh(zc), math.tanh(zc - 1.96 * se), math.tanh(zc + 1.96 * se), Q, float(chi2.sf(Q, dof))


def recorded_pooled():
    """sec. 4.6 item 4b's printed pooled value per (campaign, eta): '+0.84 (p 5.2e-62), 225'"""
    out = {}; on = False
    for l in open(REC, errors="ignore"):
        if l.startswith("## Item 4b"): on = True; continue
        if on and l.startswith("## "): break
        if on and l.startswith("| ") and not l.startswith("| campaign |") and not l.startswith("|---"):
            p = [x.strip() for x in l.strip().strip("|").split("|")]
            m = re.match(r"([+-]\d\.\d\d) \(p (\S+)\), (\d+)", p[-1])
            if m: out[(p[0], p[1])] = (m.group(1), m.group(2), int(m.group(3)))
    return out


def main():
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 10
    print("# Decision 3 item 1: the within-mass aging check (261012 sec. 4.7.4; EXPLORATORY, POST HOC, no verdict rule), printed by "
          "validation/paper1_window_aging_261009.py\n")
    # ---------------------------------------------------------------- the cells and trajectories of sec. 4.6 (its own code path)
    cells = W.discover(); inv = []
    for cell, ms in cells.items():
        M0 = sorted(ms)[0]; mt = W.meta(cell, sorted(ms[M0])[0][1])
        if mt is None: continue
        L0, en, per, Ns = mt["L0"], mt["eta"], mt["per"], mt["Ns"]
        et = en * L0 / (L0 - box_delta(L0) / 2)
        if not 0.69 <= et <= 0.725: continue
        d = W.date_of(cell)
        if not (d[:10] >= W.NEW_CORE and len(ms) >= 3 and mt["complete"] and Ns > 0): continue
        pr = W.psi6_runs(cell)
        if not pr: continue
        bad, _ = W.bad_runs(cell)
        inv.append(dict(cell=cell, camp=W.campaign_of(cell), eta=et, Ns=Ns, ms=ms, bad=bad, pr=pr))
    inv.sort(key=lambda q: (q["camp"], q["eta"]))
    jobs = []
    for q in inv:
        for M, runs in q["ms"].items():
            for r, p in sorted(runs):
                if (M, r) in q["bad"] or not T.trace_check(p)[0]: continue
                jobs.append((q["cell"], M, r, p))
    with Pool(workers) as pool:
        out = pool.map(traj2, [p for _, _, _, p in jobs], chunksize=8)
    res = defaultdict(lambda: defaultdict(list)); refused = defaultdict(int)
    for (cell, M, r, p), o in zip(jobs, out):
        if o == "refused": refused[cell] += 1; continue
        if o is not None: o["r"] = r; res[cell][M].append(o)
    if refused:
        print("refused by the loader provenance guard (accelerated backend; sec. 4.7.1, 4.7.5): " +
              ", ".join(f"{W.campaign_of(c)} {os.path.basename(c)}: {n}" for c, n in sorted(refused.items())) + "\n")
    rec = recorded_pooled()
    # ---------------------------------------------------------------- table 1: per mass, combined, pooled
    print("## 1. Per mass: Spearman(per-trajectory nu, per-trajectory psi6 run mean); combined over the masses; the pooled value of sec. 4.6 item 4b\n")
    print("| campaign | eta_true | in the window | seeds per mass | Spearman per mass, M ascending | combined over masses [95 % interval] | "
          "Cochran Q (dof, p) | pooled, recomputed (n) | pooled, recorded in sec. 4.6 |\n|---|---|---|---|---|---|---|---|---|")
    gate_bad = 0; rows_split = []; compared = set()
    for q in inv:
        R = res[q["cell"]]; rs, ns, pairs_nu, pairs_ps = [], [], [], []
        sp_d, n_d, dnu_all, dnu_zero = [], [], [], 0
        for M in sorted(R):
            tr = [o for o in R[M] if (M, o["r"]) in q["pr"]]
            if len(R[M]) < 3: continue
            if len(tr) < 3: rs.append(float("nan")); ns.append(len(tr)); continue
            nu = np.array([o["nu"] for o in tr]); ps = np.array([q["pr"][(M, o["r"])][2] for o in tr], float)
            dps = np.array([q["pr"][(M, o["r"])][1] - q["pr"][(M, o["r"])][0] for o in tr], float)
            s = spearmanr(nu, ps); rs.append(float(s.correlation)); ns.append(len(tr))
            pairs_nu += list(nu / nu.mean() - 1); pairs_ps += list(ps - ps.mean())
            dn = np.array([o["nu2"] / o["nu1"] - 1 for o in tr]); dnu_all += list(dn); dnu_zero += int((dn == 0).sum())
            sd = spearmanr(dn, dps); sp_d.append(float(sd.correlation)); n_d.append(len(tr))
        if not rs: continue
        rc, lo, hi, Q, pQ = fisher(rs, ns)
        sp = spearmanr(pairs_nu, pairs_ps)
        key = (q["camp"], f"{q['eta']:.4f}"); rr = rec.get(key)
        mine = f"{sp.correlation:+.2f}"
        if rr is not None:
            compared.add(key)
            if mine != rr[0] or len(pairs_nu) != rr[2]: gate_bad += 1
        inside = W.WIN_LO <= q["eta"] <= W.WIN_HI
        nsd = sorted(set(ns))
        print(f"| {q['camp']} | {q['eta']:.4f} | {'yes' if inside else 'no'} | {nsd[0]}{'' if len(nsd) == 1 else '-' + str(nsd[-1])} | "
              f"{', '.join('%+.2f' % v if v == v else 'nan' for v in rs)} | {rc:+.2f} [{lo:+.2f}, {hi:+.2f}] | {Q:.1f} ({len([v for v in rs if v == v]) - 1}, {pQ:.2g}) | "
              f"{mine} (p {sp.pvalue:.1e}), {len(pairs_nu)} | {('%s (p %s), %d' % rr) if rr else 'not in the recorded table'} |")
        TDs = sorted({o["TD"] for M in R for o in R[M]})
        dn = np.array(dnu_all); dc, dlo, dhi, dQ, dpQ = fisher(sp_d, n_d)
        rows_split.append((q, TDs, dn, dnu_zero, (dc, dlo, dhi, dQ, dpQ), len(sp_d)))
    gone = sorted(k for k in rec if k not in compared)
    print(f"\ngate: the pooled values recomputed here against sec. 4.6's printed table: {'IDENTICAL to the printed digits in every cell compared' if gate_bad == 0 else f'**{gate_bad} cells DIFFER**'} "
          f"({len(compared)} of its {len(rec)} rows compared)")
    if gone:
        why = (" (refused by the loader provenance guard: the accelerated backend, withdrawn in sec. 4.7.1)"
               if all(c.startswith("validate_acc") for c, e in gone) else " (**not explained**)")
        print("rows of sec. 4.6's table not recomputed here: " + ", ".join(f"{c} {e}" for c, e in gone) + why)
    # ---------------------------------------------------------------- table 2: the split
    print("\n## 2. First half against second half of each trajectory: dnu = nu2/nu1 - 1 (same estimator on TD/2 periods each)\n")
    print("| campaign | eta_true | in the window | TD (periods per half) | argmax resolution of a half | trajectories | mean dnu [%] (SE) | "
          "dnu exactly 0 (one bin) | within-mass Spearman(dnu, psi6 end - hold), combined [95 % interval] | Cochran Q (dof, p) |\n"
          "|---|---|---|---|---|---|---|---|---|---|")
    for q, TDs, dn, nz, f, nm in rows_split:
        inside = W.WIN_LO <= q["eta"] <= W.WIN_HI
        se = dn.std(ddof=1) / math.sqrt(len(dn)) if len(dn) > 1 else float("nan")
        print(f"| {q['camp']} | {q['eta']:.4f} | {'yes' if inside else 'no'} | {', '.join(str(t) for t in TDs)} ({', '.join('%g' % (t / 2) for t in TDs)}) | "
              f"{', '.join('%.1f %%' % (200.0 / t) for t in TDs)} | {len(dn)} | {100 * dn.mean():+.2f} ({100 * se:.2f}) | {nz} ({100.0 * nz / max(1, len(dn)):.0f} %) | "
              f"{f[0]:+.2f} [{f[1]:+.2f}, {f[2]:+.2f}] | {f[3]:.1f} ({nm - 1}, {f[4]:.2g}) |")
    print("\n(EXPLORATORY: no verdict rule. A positive mean dnu means the frequency rises from the first to the second half of the record; "
          "psi6 end - hold is the structural change over the whole record, the only per-run time information psi6 has.)")


if __name__ == "__main__":
    main()
