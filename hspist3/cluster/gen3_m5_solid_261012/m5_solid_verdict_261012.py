#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 38 item 6; 261012 sec. 4.7.80): the registered M5 rule at the SOLID cell, read on the Mac
from the fetched job folder (fetch_m5_solid_koa.sh), without running anything. The rule and its reading are K2-2's
(cluster/gen3_koa_261009/koa_k2_m5_rates_261010.py, item 4, whose functions are imported, not copied):
  1. the contact maxima ([EDMD-CONTACT]: disk-disk, outer walls, divider) of the koa (double) and the koa-ld (80-bit) run of each seed:
     both at or below item 3's floor 1e-6 px, and their difference too;
  2. the energy drift ([EDMD3-LEDGER] energy residual E - E0 - W): |koa| <= |koa-ld| + the koa ledger's own rounding scale (item 3's
     floor is a length, so this reading is CC's [INTERPRETATION], as in K2-2; the residuals are printed for any other reading);
  3. the phase: psi6_global over the records, against the phase constants (fluid < 0.70, solid >= 0.72);
  4. the decorrelation time: the first psi6(t) sample where the two runs of a seed differ, and where they differ by 0.05;
  then argmax nu per build (reduce_B's estimator; information, not a component of the rule, as in K2-2).
usage (from hspist3/): python3 cluster/gen3_m5_solid_261012/m5_solid_verdict_261012.py [experiments_gen3_m5_solid_koa_261012]
"""
import glob, math, os, re, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import koa_k2_m5_rates_261010 as K                                   # RULE, ITEM3, FLOOR, contact(), ledger(), psi6(), decorrelation()
import g3_rate_ld_koa as RL                                           # record(), argmax_nu()
N, M, ETA, SEEDS = 400, 500, 0.74, 4


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HS, "experiments_gen3_m5_solid_koa_261012"))
    jobs = sorted(glob.glob(os.path.join(root, "job_*")))
    if len(jobs) != 1: sys.exit(f"STOP: {root} holds {len(jobs)} job folders ({', '.join(os.path.basename(j) for j in jobs) or 'none'}), not one")
    J = jobs[0]; cell = os.path.join(J, "ld", f"N{N}_M{M}_eta{ETA:.2f}")
    log = open(os.path.join(J, "job_log.out"), errors="ignore").read() if os.path.exists(os.path.join(J, "job_log.out")) else ""
    node = re.search(r"^M5 solid cell on (.*?): out ", log, re.M); bins = re.search(r"^binaries: (.*)$", log, re.M)
    print(f"# M5, the solid cell (N {N}, eta {ETA}, M {M}), read from {os.path.relpath(J, HS)}\n")
    print(f"node: {node.group(1) if node else '(no job log)'}; {bins.group(1) if bins else ''}\n")
    print("The registered rule (261012 sec. 4.7, '5. Gate plan', item 5, verbatim):\n> " + K.RULE + "\nwith item 3, its floor:\n> " + K.ITEM3 + "\n")
    runs = {}
    for s in range(SEEDS):
        for tag in ("koa", "koa-ld"):
            d = os.path.join(cell, f"seed{s}_{tag}"); p = os.path.join(d, "run.log")
            t = open(p, errors="ignore").read() if os.path.exists(p) else ""
            runs[(s, tag)] = {"rec": RL.record(d) if t else None, "contact": K.contact(t), "ledger": K.ledger(t), "psi6": K.psi6(d),
                              "nu": RL.argmax_nu(d) if t else math.nan}
    print("| seed | build | clean | events | run [s] | contact max disk / wall / divider [px] | energy residual (ratio to its rounding scale) | "
          "gap / u_t | argmax nu | psi6 mean over the record |\n|---|---|---|---|---|---|---|---|---|---|")
    for (s, tag), q in runs.items():
        r, c, l, p = q["rec"], q["contact"], q["ledger"], q["psi6"]
        if not r: print(f"| {s} | {tag} | no run record | | | | | | | |"); continue
        rp = [x for (_, ph, x) in (p or []) if ph != "hold"]
        print(f"| {s} | {tag} | {r['clean']} | {r['events']} | {float(r['run_s']):.1f} | {'/'.join(f'{x:.2e}' for x in c[1:]) if c else '-'} | "
              f"{f'{l[0]:+.2e} ({l[2]:.2g})' if l else '-'} | {r['gap'] / r['u_t']:.3g} | {q['nu']:.6f} | {f'{statistics.mean(rp):.3f}' if rp else '-'} |")
    complete = all(q["rec"] and q["contact"] and q["ledger"] and q["psi6"] for q in runs.values())
    clean = complete and all(q["rec"]["clean"] == "1" for q in runs.values())
    print("\n## The rule, component by component\n")
    oc, oe, dec = [], [], []
    for s in range(SEEDS):
        a, b = runs[(s, "koa")], runs[(s, "koa-ld")]
        if a["contact"] and b["contact"]:
            dc = max(abs(x - y) for x, y in zip(a["contact"][1:], b["contact"][1:])); mx = max(a["contact"][1:] + b["contact"][1:])
            oc.append((s, dc, mx, dc <= K.FLOOR and mx <= K.FLOOR))
        if a["ledger"] and b["ledger"]:
            oe.append((s, a["ledger"][0], b["ledger"][0], a["ledger"][1], abs(a["ledger"][0]) <= abs(b["ledger"][0]) + a["ledger"][1]))
        if a["psi6"] and b["psi6"]: dec.append((s,) + K.decorrelation(a["psi6"], b["psi6"]))
    c_ok = len(oc) == SEEDS and all(x[3] for x in oc); e_ok = len(oe) == SEEDS and all(x[4] for x in oe)
    print(f"1. Contact maxima: the largest koa minus koa-ld difference over the {SEEDS} seed pairs is "
          f"{max((x[1] for x in oc), default=math.nan):.2e} px and the largest maximum {max((x[2] for x in oc), default=math.nan):.2e} px, "
          f"against the floor {K.FLOOR:g} px: {'within, every pair' if c_ok else 'NOT within: ' + str([x for x in oc if not x[3]])}.")
    print(f"2. Energy drift: koa {', '.join(f'{x[1]:+.1e}' for x in oe)}; koa-ld {', '.join(f'{x[2]:+.1e}' for x in oe)}; read as |koa| <= "
          f"|koa-ld| + the koa ledger's rounding scale ({', '.join(f'{x[3]:.1e}' for x in oe)}) [INTERPRETATION]: "
          f"{'within, every pair' if e_ok else 'NOT within: ' + str([x for x in oe if not x[4]])}.")
    ps = [x for q in runs.values() for (_, ph, x) in (q["psi6"] or []) if ph != "hold"]
    kind = "fluid" if ETA < K.ETA_FLUID_MAX else ("solid" if ETA >= K.ETA_SOLID else "between the phase constants")
    print(f"3. The cell is {kind} by the phase constants (fluid < {K.ETA_FLUID_MAX}, solid >= {K.ETA_SOLID}); psi6_global over the records: "
          + (f"mean {statistics.mean(ps):.3f}, range {min(ps):.3f}-{max(ps):.3f} ({len(ps)} samples)." if ps else "no samples."))
    print("4. Decorrelation time (psi6_global(t), both builds from the same start): first printed difference / first difference >= 0.05, "
          "sigma-time: " + ("; ".join(f"seed {x[0]}: {x[1]} / {x[2]}" for x in dec) if dec else "none"))
    print("\n## argmax nu (information, not a component of the rule)\n")
    v = {tag: [runs[(s, tag)]["nu"] for s in range(SEEDS)] for tag in ("koa", "koa-ld")}
    if all(math.isfinite(x) for x in v["koa"] + v["koa-ld"]):
        m = {k: statistics.mean(x) for k, x in v.items()}; se = {k: statistics.stdev(x) / math.sqrt(SEEDS) for k, x in v.items()}
        dif = (m["koa-ld"] - m["koa"]) / m["koa"] * 100; sed = math.hypot(se["koa"], se["koa-ld"]) / m["koa"] * 100
        print(f"koa {m['koa']:.6f} +- {se['koa']:.6f}; koa-ld {m['koa-ld']:.6f} +- {se['koa-ld']:.6f}; koa-ld minus koa {dif:+.2f} % "
              f"(z {dif / sed if sed > 0 else math.nan:+.2f}, n = {SEEDS} each)")
    else:
        print("argmax nu missing in some run")
    print(f"\nM5 AT THE SOLID CELL, BY THE REGISTERED RULE: " + (
        "the double run's contact maxima and energy drift are within the 80-bit run's, every seed pair (all 8 runs clean)" if clean and c_ok and e_ok
        else "NOT SHOWN: " + ", ".join(w for w, ok in (("a run is missing or not clean", clean), ("the contact maxima", c_ok), ("the energy drift", e_ok)) if not ok)))


if __name__ == "__main__":
    main()
