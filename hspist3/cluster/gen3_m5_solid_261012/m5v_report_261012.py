#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 6; 261012 sec. 4.7.81): REPORT M5-V -- the registered M5 rule (sec. 4.7, gate plan item 5)
on K2-2's two cells and the solid cell together, and the verdict for stage 1 gate G6 (b). Read-only, nothing runs.
  the cells    K2-2 part (d) (job 15043030, cn-09-30-02, Gold 6342; fetched, experiments_gen3_gate_koa_261009/rate_ld_k2 of the main worktree):
               N 100 at pi/8, M 300 (fluid) and N 400 at eta 0.70, M 500 (the fluid limit by the phase constants); the M5 job (job
               15051694, cn-05-03-03, E5-2670 v3; fetched, experiments_gen3_m5_solid_koa_261012 of the engine-gen3 worktree): N 400 at eta
               0.74, M 500 (solid). Four seeds each, the koa (double) and the koa-ld (80-bit) build of a643bc7, the same commands and seeds.
  the rule     per seed pair, by K2-2's reader (koa_k2_m5_rates_261010: contact(), ledger(), psi6(), decorrelation(); its reading of the
               energy part is CC's [INTERPRETATION], printed as such): 1. the contact maxima (disk-disk, outer walls, divider) of both runs
               at or below item 3's floor 1e-6 px, and their difference too; 2. the energy drift |koa| <= |koa-ld| + the koa ledger's
               rounding scale; 3. the phase; 4. the decorrelation time (reported).
  the verdict  G6 (b) asks that the M5 verdict "does not say DIFFERENCE". Here the verdict says NO DIFFERENCE when components 1 and 2 hold
               for every seed pair of the fluid cell and of the solid cell (all runs clean), DIFFERENCE otherwise, naming the pair.
  the proviso  decision 39 item 5 b: the M5 job on the a643bc7 pair is the gate PROVIDED every source compiled into the binary except
               00ALLINONE.c is byte-identical between a643bc7 and f72d764: the conclusion line of
               experiments_gen3_s1_261012/g6b_identity_output.txt is printed.
  argmax nu    per cell and build (reduce_B's estimator), koa-ld minus koa with z: information, not a component of the rule.
usage (from hspist3/): python3 cluster/gen3_m5_solid_261012/m5v_report_261012.py
"""
import math, os, re, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); MAIN_HS = os.path.join(os.path.dirname(os.path.dirname(HS)), "HardDisks", "hspist3")
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import koa_k2_m5_rates_261010 as K
import g3_rate_ld_koa as RL
K2 = os.path.join(MAIN_HS, "experiments_gen3_gate_koa_261009", "rate_ld_k2", "ld")
M5 = os.path.join(HS, "experiments_gen3_m5_solid_koa_261012", "job_15051694", "ld")
CELLS = [("fluid", 100, math.pi / 8, 300, os.path.join(K2, "N100_M300"), "K2-2, job 15043030"),
         ("the fluid limit", 400, 0.70, 500, os.path.join(K2, "N400_M500"), "K2-2, job 15043030"),
         ("solid", 400, 0.74, 500, os.path.join(M5, "N400_M500_eta0.74"), "M5, job 15051694")]


def load(d):
    out = {}
    for s in range(4):
        for tag in ("koa", "koa-ld"):
            p = os.path.join(d, f"seed{s}_{tag}"); f = os.path.join(p, "run.log")
            t = open(f, errors="ignore").read() if os.path.exists(f) else ""
            out[(s, tag)] = dict(rec=RL.record(p) if t else None, contact=K.contact(t), ledger=K.ledger(t), psi6=K.psi6(p), nu=RL.argmax_nu(p) if t else math.nan)
    return out


def main():
    print("# Report M5-V: the registered M5 rule on the fluid cells (K2-2) and the solid cell together (decision 39 item 6), printed by "
          "cluster/gen3_m5_solid_261012/m5v_report_261012.py\n")
    print("The registered rule (261012 sec. 4.7, '5. Gate plan', item 5, verbatim):\n> " + K.RULE + "\nwith item 3, its floor:\n> " + K.ITEM3 + "\n")
    print("| cell | phase (constants: fluid < 0.70, solid >= 0.72) | source | seed | build | clean | events | contact max disk / wall / divider [px] | "
          "energy residual (ratio to its scale) | gap / u_t | argmax nu | psi6 over the record |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    verdict, comps, info = [], [], []
    for ph, N, eta, M, d, src in CELLS:
        R = load(d)
        for (s, tag), q in R.items():
            r, c, l, p = q["rec"], q["contact"], q["ledger"], q["psi6"]
            rp = [x for (_, phs, x) in (p or []) if phs != "hold"]
            print(f"| N {N}, eta {eta:.4g}, M {M} | {ph} | {src} | {s} | {tag} | {r['clean'] if r else '-'} | {r['events'] if r else '-'} | "
                  f"{'/'.join(f'{x:.2e}' for x in c[1:]) if c else '-'} | {f'{l[0]:+.2e} ({l[2]:.2g})' if l else '-'} | "
                  f"{(r['gap'] / r['u_t']) if r else float('nan'):.3g} | {q['nu']:.6f} | {statistics.mean(rp) if rp else float('nan'):.3f} |")
        clean = all(q["rec"] and q["rec"]["clean"] == "1" for q in R.values())
        oc, oe, dec = [], [], []
        for s in range(4):
            a, b = R[(s, "koa")], R[(s, "koa-ld")]
            if a["contact"] and b["contact"]:
                dc = max(abs(x - y) for x, y in zip(a["contact"][1:], b["contact"][1:])); mx = max(a["contact"][1:] + b["contact"][1:])
                oc.append((s, dc, mx, dc <= K.FLOOR and mx <= K.FLOOR))
            if a["ledger"] and b["ledger"]:
                oe.append((s, a["ledger"][0], b["ledger"][0], a["ledger"][1], abs(a["ledger"][0]) <= abs(b["ledger"][0]) + a["ledger"][1]))
            if a["psi6"] and b["psi6"]: dec.append((s,) + K.decorrelation(a["psi6"], b["psi6"]))
        c_ok = len(oc) == 4 and all(x[3] for x in oc); e_ok = len(oe) == 4 and all(x[4] for x in oe)
        comps.append((ph, N, eta, clean, c_ok, e_ok, oc, oe, dec))
        v = {t: [R[(s, t)]["nu"] for s in range(4)] for t in ("koa", "koa-ld")}
        if all(math.isfinite(x) for x in v["koa"] + v["koa-ld"]):
            m = {t: statistics.mean(x) for t, x in v.items()}; se = {t: statistics.stdev(x) / 2 for t, x in v.items()}
            dd = (m["koa-ld"] - m["koa"]) / m["koa"] * 100; sd = math.hypot(se["koa"], se["koa-ld"]) / m["koa"] * 100
            info.append((ph, N, eta, m["koa"], se["koa"], m["koa-ld"], se["koa-ld"], dd, dd / sd))
    print("\n## The rule, component by component, per cell\n")
    print("| cell | phase | all 8 runs clean | 1. contact maxima: largest difference / largest maximum [px] (floor 1e-6) | 2. energy drift: koa | koa-ld | "
          "the koa ledger's scale | within | 4. decorrelation: first printed difference / first difference >= 0.05 [sigma-time] |\n|---|---|---|---|---|---|---|---|---|")
    for ph, N, eta, clean, c_ok, e_ok, oc, oe, dec in comps:
        print(f"| N {N}, eta {eta:.4g} | {ph} | {'yes' if clean else '**NO**'} | {max(x[1] for x in oc):.2e} / {max(x[2] for x in oc):.2e}: "
              f"{'within' if c_ok else '**NOT within**'} | {', '.join(f'{x[1]:+.1e}' for x in oe)} | {', '.join(f'{x[2]:+.1e}' for x in oe)} | "
              f"{max(x[3] for x in oe):.1e} | {'within' if e_ok else '**NOT within**'} [INTERPRETATION] | " + "; ".join(f"seed {x[0]}: {x[1]} / {x[2]}" for x in dec) + " |")
    print("\n## argmax nu (information, not a component of the rule)\n\n| cell | phase | koa mean +- SE | koa-ld mean +- SE | koa-ld minus koa [%] | z |\n|---|---|---|---|---|---|")
    for ph, N, eta, mk, sk, ml, sl, dd, z in info:
        print(f"| N {N}, eta {eta:.4g} | {ph} | {mk:.6f} +- {sk:.6f} | {ml:.6f} +- {sl:.6f} | {dd:+.2f} | {z:+.2f} |")
    g6 = open(os.path.join(HS, "experiments_gen3_s1_261012", "g6b_identity_output.txt")).read()
    line = next(l for l in g6.splitlines() if l.startswith("EVERY COMPILED SOURCE"))
    print(f"\nthe proviso of decision 39 item 5 b (experiments_gen3_s1_261012/g6b_identity_output.txt): {line}")
    need = [c for c in comps if c[0] in ("fluid", "solid")]
    ok = all(c[3] and c[4] and c[5] for c in need)
    bad = [f"N {c[1]}, eta {c[2]:.4g}" for c in need if not (c[3] and c[4] and c[5])]
    print(f"\nTHE M5 VERDICT (gate G6 b), the registered rule at one fluid and one solid cell: "
          + ("NO DIFFERENCE -- the double run's contact maxima and energy drift are within the 80-bit run's for every seed pair of the fluid cell "
             "(N 100, pi/8) and of the solid cell (N 400, eta 0.74), all runs clean" if ok else "DIFFERENCE at " + ", ".join(bad))
          + "; the cell at the fluid limit (N 400, eta 0.70): " + ("within" if all(c[3] and c[4] and c[5] for c in comps if c[0] == "the fluid limit") else "NOT within"))


if __name__ == "__main__":
    main()
