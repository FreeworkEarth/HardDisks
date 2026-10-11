#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 20, task 1; 261012 sec. 4.7.47): the unit test of g3_rate_ld_koa.record() on
  (1) a REAL run record and [EDMD3-GAP] line printed by the gen3 driver (engine-gen3's experiments_gen3_m3_261009/evidence_f42befb/
      acc4/afix/run.log, stage A, HD_CONTACT_AUDIT=1), copied below unchanged;
  (2) THE LINE THAT STOPPED KOA's gate job 15040861: the same line with the divider's value 1.45e-12 -- KOA printed the failing string
      '1.45e-12,'; its run log stays on KOA until the fetch, so the line is rebuilt around that value in the driver's printf format;
  (3) every position (pair, wall, divider, piston, u_t) with the other values printf's %.3g / %.6g can print: 0, 1e-05, 1.5e+03, inf,
      nan, 3.68e-12.
For each line it checks that the OLD pattern (a643bc7, copied below) fails as on KOA -- its last group takes the comma, so float()
raises "could not convert string to float: '<divider>,'" -- and that record() of this commit returns the largest of pair, wall and
divider, the engine's u_t and the run record's fields. The test writes its run logs into a new folder under the temporary directory
(TMPDIR) and leaves it there.
usage: python3 hspist3/cluster/gen3_koa_261009/test_g3_rate_ld_record_261010.py        (exit 0 and "ALL PASS" when every check holds)
"""
import math, os, re, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g3_rate_ld_koa as R

OLD = r"\[EDMD3-GAP\] .*?max gap \[px\] pair (\S+), wall (\S+), divider (\S+)"      # g3_rate_ld_koa.py line 44 at a643bc7
REAL_HEALTH = ("[EDMD3-HEALTH] energy-transfer seed=9700: clean=1 overlap_repair=0 wall_overdue=0 past_event=0 clamp_repair=0 cell_repair=0 "
               "grid_escape=0 stagnation=0 obj_overlap_repair=0 body_findings=0 local_findings=0 full_findings=0 contact_now=0 obj_contact_now=0 "
               "wall_contact_now=0 local_checks=3141958 full_checks=10456 local_worst=-3.5e-12 full_worst=-2.37e-12 contact_c_min=0 "
               "obj_contact_gap_min=0 cross_residual_max=3.35e-12 ev_pair=1097449 ev_wall=199212 ev_cross=681479 ev_div=66369 ev_piston=0 "
               "ev_band=0 ev_stale=5594624 origin_shifts=15 syncs=15 heap_compactions=0 heap_max=1809 validator_every=60 hash=4428249ad3af787f "
               "t_end=125279.60186681151 engine=gen3 build=\"f42befb\" target=mac-O3-e0pre run_s=1.646 engine_s=1.553 driver_share=0.0569")
REAL_GAP = ("[EDMD3-GAP] contact audit over 1363030 events: max gap [px] pair 3.5e-12, wall 1.76e-12, divider 1.71e-12, piston 0; time "
            "quantum u_t 1.81899e-12 units; max gap / u_t [px/unit] pair 1.92, wall 0.969, divider 0.938, piston 0; v_ref 20 px/unit")


def gap_line(pair, wall, div, pis, ut):
    """the driver's printf (00ALLINONE.c, g3_run_extras) with the given printed values"""
    f = lambda s: float(s)
    q = lambda a: "%.3g" % (f(a) / f(ut)) if f(ut) not in (0.0,) else "inf"
    return (f"[EDMD3-GAP] contact audit over 1363030 events: max gap [px] pair {pair}, wall {wall}, divider {div}, piston {pis}; time quantum "
            f"u_t {ut} units; max gap / u_t [px/unit] pair {q(pair)}, wall {q(wall)}, divider {q(div)}, piston {q(pis)}; v_ref 20 px/unit")


def check(name, gap, base, fails):
    d = tempfile.mkdtemp(prefix=f"g3rate_{name}_", dir=base)
    with open(os.path.join(d, "run.log"), "w") as fh: fh.write("[EDMD3] built #1 at t = 0\n" + gap + "\n" + REAL_HEALTH + "\n")
    m = re.search(OLD, gap)
    try:
        max(float(x) for x in m.groups()); old = "parsed (no error)"
    except ValueError as e:
        old = f"ValueError: {e}"
    vals = re.search(r"pair (\S+?), wall (\S+?), divider (\S+?), piston (\S+?); time quantum u_t (\S+) units", gap).groups()
    want_gap, want_ut = max(float(v) for v in vals[:3]), float(vals[4])
    kv = R.record(d)
    same = lambda a, b: (math.isnan(a) and math.isnan(b)) or a == b
    ok_new = kv is not None and same(kv["gap"], want_gap) and same(kv["u_t"], want_ut) and kv["clean"] == "1" and kv["hash"] == "4428249ad3af787f" \
        and kv["events"] == 1097449 + 199212 + 681479 + 66369
    ok_old = old.startswith("ValueError") == fails and (not fails or old.endswith(f"'{vals[2]},'"))
    print(f"| {name} | divider '{vals[2]}' | old pattern: {old} | record(): gap {kv['gap'] if kv else None}, u_t {kv['u_t'] if kv else None}, "
          f"clean {kv['clean'] if kv else None}, events {kv['events'] if kv else None} | {'PASS' if ok_new and ok_old else '**FAIL**'} |")
    return ok_new and ok_old


def main():
    base = tempfile.mkdtemp(prefix="test_g3_rate_ld_record_")
    print(f"# record() of g3_rate_ld_koa.py (this commit) against the old pattern of a643bc7; run logs in {base}\n")
    print("| case | the divider's printed value | the old pattern (a643bc7) | record() of this commit | check |\n|---|---|---|---|---|")
    ok = check("real_line_stageA_acc4_afix", REAL_GAP, base, True)
    ok &= check("koa_15040861_divider_1.45e-12", gap_line("3.5e-12", "1.76e-12", "1.45e-12", "0", "1.81899e-12"), base, True)
    vals = ("0", "1e-05", "1.5e+03", "inf", "nan", "3.68e-12")
    for i, v in enumerate(vals):
        w = vals[(i + 1) % len(vals)]; x = vals[(i + 2) % len(vals)]; y = vals[(i + 3) % len(vals)]
        ok &= check(f"values_{i}", gap_line(v, w, x, y, "1.81899e-12"), base, True)
    ok &= check("u_t_long_double", gap_line("1e-15", "2e-15", "3e-15", "0", "8.88178e-16"), base, True)
    print(f"\n{'ALL PASS' if ok else '**FAILED**'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
