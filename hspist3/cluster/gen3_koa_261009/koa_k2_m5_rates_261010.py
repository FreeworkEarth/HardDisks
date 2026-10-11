#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 22, items 4 and 5; 261012 sec. 4.7.49): K2-2's output (job 15043030, gate parts (c) and (d), on
cn-09-30-02, Intel Xeon Gold 6342 2.80 GHz per decision 22), fetched to the Mac, read without running anything:
  ITEM 4  the registered M5 rule (sec. 4.7, gate plan item 5, quoted below) on the 16 rows of part (d): per seed the koa (double) run
          against the koa-ld (80-bit) run -- the contact maxima ([EDMD-CONTACT], item 3's floor 1e-6 px), the energy drift ([EDMD3-LEDGER]),
          the cell's phase (psi6 over the record), the decorrelation time (the first psi6(t) sample where the two runs differ, and where
          they differ by 0.05), clean, gap / u_t per build ([EDMD3-GAP]); then argmax nu (reduce_B's estimator, information): per cell and
          build mean and SE over the 4 seeds, the difference and its z, and the seeds per build a comparison would need to resolve a
          stated difference, with its cost on one KOA core.
  ITEM 5  the rate table of part (c) as measured on the Gold node; the E5/Gold ratio at N 400, eta 0.70 against the step-4 cross-node
          cell on an E5-2670 v3 node (Chris's run record, sec. 4.7.47: 18.9 M events in 51.5 s of engine time, 0.37 M/s; whole run
          0.32 M/s); the slowest real case against the 2e5 events/s planning value.
usage (from hspist3/): python3 cluster/gen3_koa_261009/koa_k2_m5_rates_261010.py [experiments_gen3_gate_koa_261009/rate_ld_k2]
"""
import glob, math, os, re, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import g3_rate_ld_koa as RL                                            # record() and argmax_nu() as K2-2 used them

RULE = ("5. **Long-double spot check:** the double run's contact maxima and energy drift within the 80-bit run's by the floor of item 3, at "
        "one fluid and one solid cell; decorrelation time reported.")
ITEM3 = "3. **Contact audit:** ≤ 1e-6 px in every class and run (expected ~1e-11, as today)."
FLOOR = 1e-6                                                           # px, item 3
E5 = {"engine": 18.9e6 / 51.5, "whole": 0.32e6}                        # step 4, N 400, eta 0.70, M 500 (sec. 4.7.47)
PLAN = 2e5                                                             # the planning value (events/s)
ETA_FLUID_MAX, ETA_SOLID = 0.700, 0.720                                # plot_speed_of_sound_edmd.py's phase constants


def contact(t):
    m = re.search(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), outer walls (\S+), divider (\S+),", t)
    return (int(m.group(1)), *[float(x) for x in m.groups()[1:]]) if m else None


def ledger(t):
    m = re.search(r"\[EDMD3-LEDGER\] .*?energy: residual (\S+), scale (\S+), ratio (\S+) \(E0 (\S+), E (\S+), W (\S+)\)", t)
    return [float(x) for x in m.groups()] if m else None


def psi6(d):
    fs = glob.glob(os.path.join(d, "psi6_t_*_run0.csv")) + glob.glob(os.path.join(d, "psi6_t_*_run0.csv.gz"))
    if not fs: return None
    import gzip
    op = gzip.open if fs[0].endswith(".gz") else open
    rows = [l.strip().split(",") for l in op(fs[0], "rt")][1:]
    return [(float(r[0]), r[1], float(r[2])) for r in rows if len(r) >= 3]


def decorrelation(a, b):
    """the first sample time where the psi6_global series of two runs differ (at the printed 6 decimals), and where they differ by 0.05"""
    t1 = t5 = None
    for (ta, pa, xa), (tb, pb, xb) in zip(a, b):
        if ta != tb: return "sample times differ", None
        if t1 is None and xa != xb: t1 = ta
        if t5 is None and abs(xa - xb) >= 0.05: t5 = ta; break
    return t1, t5


def main():
    D = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HS, "experiments_gen3_gate_koa_261009", "rate_ld_k2"))
    if not os.path.isdir(D): sys.exit(f"STOP: {D} is not there (fetch K2-2's folder first)")
    print(f"# K2-2 (job 15043030; cn-09-30-02, Intel Xeon Gold 6342 2.80 GHz per decision 22) read from {D}\n")
    # ------------------------------------------------------------------ item 5: the rate table on the Gold node
    print("## Item 5: part (c), the rate table as measured on the Gold node (one run at a time; 500 sigma-time held, 100 released, M = 300)\n")
    print("| N | eta | events | run [s] | engine [s] | driver share | events/s (run) | events/s (engine) | clean |\n|---|---|---|---|---|---|---|---|---|")
    rate = {}
    for N in (100, 400, 900, 1600):
        for lab in ("pi8", "0.70", "0.78"):
            d = os.path.join(D, "rate", f"N{N}_{lab}"); r = RL.record(d) if os.path.exists(os.path.join(d, "run.log")) else None
            if not r: print(f"| {N} | {lab} | no run record | | | | | | |"); continue
            ev, rs, es = r["events"], float(r["run_s"]), float(r["engine_s"]); rate[(N, lab)] = (ev / rs, ev / es)
            print(f"| {N} | {lab.replace('pi8', 'pi/8')} | {ev} | {rs:.1f} | {es:.1f} | {float(r['driver_share']):.3f} | {ev / rs:.3g} | {ev / es:.3g} | {r['clean']} |")
    # ------------------------------------------------------------------ item 4: the 16 rows of part (d)
    print("\n## Item 4: M5 on KOA, the registered rule on the 16 rows of part (d)\n")
    print("The registered rule (261012 sec. 4.7, '5. Gate plan', item 5, verbatim):\n> " + RULE + "\nwith item 3, its floor:\n> " + ITEM3 + "\n")
    rows, runs = [], {}
    for N, M in ((100, 300), (400, 500)):
        for s in range(4):
            for tag in ("koa", "koa-ld"):
                d = os.path.join(D, "ld", f"N{N}_M{M}", f"seed{s}_{tag}"); t = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
                runs[(N, M, s, tag)] = {"rec": RL.record(d) if t else None, "contact": contact(t), "ledger": ledger(t), "psi6": psi6(d), "nu": RL.argmax_nu(d) if t else math.nan}
    print("| N | M | seed | build | clean | events | run [s] | contact max disk / wall / divider [px] | energy residual (ratio to its rounding scale) | gap / u_t | argmax nu | psi6 mean over the record |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (N, M, s, tag), q in runs.items():
        r, c, l, p = q["rec"], q["contact"], q["ledger"], q["psi6"]
        rec_psi = [x for (_, ph, x) in (p or []) if ph != "hold"]
        print(f"| {N} | {M} | {s} | {tag} | {r['clean'] if r else '-'} | {r['events'] if r else '-'} | {float(r['run_s']):.1f} | "
              f"{'/'.join(f'{x:.2e}' for x in c[1:]) if c else '-'} | {f'{l[0]:+.2e} ({l[2]:.2g})' if l else '-'} | "
              f"{r['gap'] / r['u_t']:.3g} | {q['nu']:.6f} | {f'{statistics.mean(rec_psi):.3f}' if rec_psi else '-'} |" if r else f"| {N} | {M} | {s} | {tag} | no run record | | | | | | | |")
    print("\n### The rule, component by component\n")
    ok_c, ok_e, dec = [], [], []
    for N, M in ((100, 300), (400, 500)):
        for s in range(4):
            a, b = runs[(N, M, s, "koa")], runs[(N, M, s, "koa-ld")]
            if a["contact"] and b["contact"]:
                dc = max(abs(x - y) for x, y in zip(a["contact"][1:], b["contact"][1:])); mx = max(a["contact"][1:] + b["contact"][1:])
                ok_c.append((N, s, dc, mx, dc <= FLOOR and mx <= FLOOR))
            if a["ledger"] and b["ledger"]:
                ok_e.append((N, s, a["ledger"][0], b["ledger"][0], a["ledger"][1], abs(a["ledger"][0]) <= abs(b["ledger"][0]) + a["ledger"][1]))
            if a["psi6"] and b["psi6"]: dec.append((N, s) + decorrelation(a["psi6"], b["psi6"]))
    print(f"1. Contact maxima (disk-disk, outer walls, divider): the largest koa minus koa-ld difference over the 8 seed pairs is "
          f"{max(x[2] for x in ok_c):.2e} px and the largest maximum {max(x[3] for x in ok_c):.2e} px, against the floor {FLOOR:g} px: "
          f"{'within, every pair' if all(x[4] for x in ok_c) and len(ok_c) == 8 else 'NOT within: ' + str([x for x in ok_c if not x[4]])}.")
    print(f"2. Energy drift ([EDMD3-LEDGER] energy residual, E - E0 - W): koa {', '.join(f'{x[2]:+.1e}' for x in ok_e)}; koa-ld "
          f"{', '.join(f'{x[3]:+.1e}' for x in ok_e)}. Read as |koa| <= |koa-ld| + the koa ledger's own rounding scale "
          f"(item 3's floor is a length; this reading is CC's [INTERPRETATION]): {'within, every pair' if all(x[5] for x in ok_e) and len(ok_e) == 8 else 'NOT within: ' + str([x for x in ok_e if not x[5]])}.")
    for N, M, eta in ((100, 300, math.pi / 8), (400, 500, 0.70)):
        ps = [x for s in range(4) for tag in ("koa", "koa-ld") for (_, ph, x) in (runs[(N, M, s, tag)]["psi6"] or []) if ph != "hold"]
        kind = "fluid" if eta < ETA_FLUID_MAX else ("solid" if eta >= ETA_SOLID else "the fluid limit (eta = ETA_FLUID_MAX)")
        print(f"3. Cell N {N}, eta {eta:.4g}: {kind} by the phase constants (fluid < {ETA_FLUID_MAX}, solid >= {ETA_SOLID}); psi6_global over the "
              f"records: mean {statistics.mean(ps):.3f}, range {min(ps):.3f}-{max(ps):.3f} ({len(ps)} samples).")
    print("4. Decorrelation time (psi6_global(t), both runs from the same start, 0.25 sigma-time samples): first printed difference / "
          "first difference >= 0.05, sigma-time: " + "; ".join(f"N {x[0]} seed {x[1]}: {x[2]} / {x[3]}" for x in dec))
    # ------------------------------------------------------------------ argmax nu (information), and the size of a larger comparison
    print("\n### argmax nu (information in part (d); not a component of the registered rule)\n")
    print("| N | build | mean argmax nu | SE (n = 4) | koa-ld minus koa | z | per-run SD (pooled) | seeds per build for SE(diff) = 0.25 % | mean run [s] on Gold, koa / koa-ld |")
    print("|---|---|---|---|---|---|---|---|---|")
    cost = 0.0
    for N, M in ((100, 300), (400, 500)):
        v = {tag: [runs[(N, M, s, tag)]["nu"] for s in range(4)] for tag in ("koa", "koa-ld")}
        m = {k: statistics.mean(x) for k, x in v.items()}; se = {k: statistics.stdev(x) / 2 for k, x in v.items()}
        dif = (m["koa-ld"] - m["koa"]) / m["koa"] * 100; sed = math.hypot(se["koa"], se["koa-ld"]) / m["koa"] * 100
        sd = math.sqrt((statistics.variance(v["koa"]) + statistics.variance(v["koa-ld"])) / 2) / m["koa"] * 100
        n = math.ceil(2 * sd ** 2 / 0.25 ** 2)
        t = {k: statistics.mean(float(runs[(N, M, s, k)]["rec"]["run_s"]) for s in range(4)) for k in ("koa", "koa-ld")}
        cost += n * (t["koa"] + t["koa-ld"])
        for k in ("koa", "koa-ld"):
            print(f"| {N} | {k} | {m[k]:.6f} | {se[k]:.6f} | " + (f"{dif:+.2f} % | {dif / sed:+.2f} | {sd:.2f} % | {n} | {t['koa']:.1f} / {t['koa-ld']:.1f} |" if k == "koa-ld" else "| | | | |"))
    f = rate.get((400, "0.70"))
    g = statistics.mean(runs[(400, 500, s, "koa")]["rec"]["events"] / float(runs[(400, 500, s, "koa")]["rec"]["engine_s"]) for s in range(4))
    print(f"\nCost of that comparison on one core: {cost / 3600:.2f} core-hours on the Gold node"
          + (f", about {cost / 3600 * g / E5['engine']:.2f} on an E5-2670 v3 node (the ratio below)" if g else ""))
    # ------------------------------------------------------------------ item 5: the E5/Gold ratio, the slowest case, the planning value
    print("\n## Item 5: E5-2670 v3 against Gold 6342 at N 400, eta 0.70\n")
    print(f"- E5 (step 4, the cross-node cell: M 500, 2000 sigma-time held, 200 periods): {E5['engine']:.3g} events/s of engine time, {E5['whole']:.3g} whole run.")
    print(f"- Gold, part (d), koa build, the same N, eta, M and seeding (2000 steps held, 200 periods), mean of 4 seeds: {g:.3g} events/s of engine time. "
          f"E5/Gold = {E5['engine'] / g:.3f} (Gold {g / E5['engine']:.2f}x faster).")
    if f: print(f"- Gold, part (c), N 400 eta 0.70 (M 300, 500 sigma-time held, 100 released): {f[1]:.3g} events/s of engine time, {f[0]:.3g} whole run; "
                f"E5/Gold = {E5['engine'] / f[1]:.3f} (engine), {E5['whole'] / f[0]:.3f} (whole run).")
    if rate:
        (Ns, ls), (r_run, r_eng) = min(rate.items(), key=lambda kv: kv[1][0])
        r5 = r_run * E5["engine"] / g
        print(f"- The slowest real case of part (c) on Gold: N {Ns}, eta {ls}, {r_run:.3g} events/s (whole run); on an E5 node about {r5:.3g} "
              f"(scaled by the ratio above) against the planning value {PLAN:.0e}: {r5 / PLAN:.2f}x of it.")


if __name__ == "__main__":
    main()
