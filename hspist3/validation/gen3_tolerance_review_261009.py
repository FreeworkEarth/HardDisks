#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.18; stage E1 of the plan-author programme of sec. 4.7.12): the tolerance review of gen3. One table
of every tolerance and cut-off in edmd_core/edmd_gen3.c and in the gen3 driver path (00ALLINONE.c): its value, its scale, its
derivation, the largest measured value today and the margin (tolerance / measured). The values are evaluated at the production run
of amendment d (N = 400, eta 0.70, v_ref = 40 px/unit, box 538 px; sec. 4.7.14); the measured values are the maxima over every gen3
run record and harness output given on the command line (each with its source).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 validation/gen3_tolerance_review_261009.py <run logs and harness outputs ...> > gen3_tolerance_review_output.txt
"""
import glob, gzip, math, os, re, sys

U_T = 2.0 ** -39                     # ulp(2^13): the time quantum (internal units; 24 units = 1 sigma-time)
K = 4.0
D, R, W = 24.0, 12.0, 32.0           # disk diameter, radius, cell width [px]
V_REF, BOXW = 40.0, 538.0            # the production run (sec. 4.7.14): sqrt(2 E (1 + 1/m_min)) = sqrt(2 * 400 * 2)


def ulp(x): return math.ulp(x)


def scan(paths):
    """Maxima over the run records ([EDMD3-HEALTH], [EDMD3-GAP]) and the contact lines ([EDMD-CONTACT]) of every file given."""
    M = {}
    def upd(k, v, src):
        if v is None or not math.isfinite(v): return
        if k not in M or v > M[k][0]: M[k] = (v, src)
    def low(k, v, src):                                    # the smallest value (the margins)
        if v is None or not math.isfinite(v): return
        if k not in M or v < M[k][0]: M[k] = (v, src)
    nrec = nfile = 0
    for p in paths:
        op = gzip.open if p.endswith(".gz") else open
        try: t = op(p, "rt", errors="ignore").read()
        except OSError: continue
        if "[EDMD3-HEALTH]" not in t: continue           # gen3 run logs only (a gen2 log's contact line is not gen3's)
        nfile += 1; src = os.path.relpath(p)
        for rid, clean, rest in re.findall(r"^\[EDMD3-HEALTH\] (.*?): clean=(\S+) (.*)$", t, re.M):
            nrec += 1
            kv = dict(re.findall(r"(\w+)=(\S+)", rest))
            f = lambda k: float(kv[k]) if k in kv else None
            upd("local_worst_neg", -f("local_worst") if f("local_worst") is not None else None, src)
            upd("full_worst_neg", -f("full_worst") if f("full_worst") is not None else None, src)
            upd("contact_c_neg", -f("contact_c_min") if f("contact_c_min") is not None else None, src)
            upd("obj_gap_neg", -f("obj_contact_gap_min") if f("obj_contact_gap_min") is not None else None, src)
            upd("cross_residual", f("cross_residual_max"), src)
            upd("stagnation", f("stagnation"), src)
            for k in ("overlap_repair", "wall_overdue", "obj_overlap_repair", "cell_repair", "local_findings", "full_findings", "body_findings"):
                upd(k, f(k), src)
            upd("t_end_shifted", None, src)
        for m in re.finditer(r"\[EDMD-CONTACT\] executed events \d+; max abs\(contact distance\) \[px\]: disk-disk (\S+), outer walls (\S+), divider (\S+), pistons (\S+)", t):
            upd("gap_pair", float(m.group(1)), src); upd("gap_wall", float(m.group(2)), src); upd("gap_div", float(m.group(3)), src); upd("gap_pis", float(m.group(4)), src)
        for m in re.finditer(r"\[EDMD3-GAP\] contact audit over \d+ events: max gap \[px\] pair (\S+), wall (\S+), divider (\S+), piston (\S+);"
                             r".*?max gap / u_t \[px/unit\] pair (\S+), wall (\S+), divider (\S+), piston (\S+); v_ref (\S+) px/unit", t):
            gp_, gw_, gd_ = float(m.group(1)), float(m.group(2)), float(m.group(3))
            upd("gap_pair", gp_, src); upd("gap_wall", gw_, src); upd("gap_div", gd_, src)
            vr = float(m.group(9))                         # this run's v_ref: its own tolerances are K v_ref u_t (x 2d for c)
            for k, q in (("m_pair", float(m.group(5))), ("m_wall", float(m.group(6))), ("m_div", float(m.group(7)))):
                if q > 0: low(k, K * vr / q, src)
    return M, nrec, nfile


def main():
    paths = []
    for a in sys.argv[1:]: paths += sorted(glob.glob(a, recursive=True)) or [a]
    M, nrec, nfile = scan(paths)
    g = lambda k: M[k][0] if k in M else None
    s = lambda k: M[k][1] if k in M else "-"
    c_tol = K * 2 * D * V_REF * U_T; tol_face = K * V_REF * U_T + 8 * ulp(BOXW)
    tol_pair = max(1e-7, 1e-6 * D); tol_wall = max(1e-6, 1e-6 * max(1.0, R)); tol_cell = 1e-9
    print("# gen3 tolerance review (261012 sec. 4.7.18, stage E1), printed by validation/gen3_tolerance_review_261009.py\n")
    print(f"Inputs: {nfile} files, {nrec} gen3 run records. Values at the production run (v_ref = {V_REF:g} px/unit, box {BOXW:g} px); "
          f"u_t = 2^-39 = {U_T:.4g} units; K = {K:g}.\n")
    rows = []
    def row(name, where, value, scale, derivation, meas, msrc, kind="tolerance"):
        margin = (value / meas) if (meas is not None and meas > 0 and isinstance(value, float)) else None
        rows.append((name, where, value, scale, derivation, meas, msrc, margin, kind))
    gp = g("gap_pair")
    row("u_t (time quantum)", "edmd_gen3.c `S->u_time = ldexp(1.0, -39)`", U_T, "units", "ulp(2^13): the origin shift keeps now < 2^13 + one event, so every "
        "event time is resolved to u_t/2", None, "-", "constant")
    row("origin shift interval", "edmd_gen3.h `EDMD3_ORIGIN_SHIFT`", 8192.0, "units (341.3 sigma-time)", "2^13: the largest time with ulp = u_t", None, "-", "constant")
    row("K", "edmd_gen3.h `EDMD3_TOL_K`", K, "-", "the derived bound is 2.5 (sec. 4.7.6); K = 4 is the stated safety factor", None, "-", "constant")
    row("v_ref", "edmd_gen3.c `tol_update`", V_REF, "px/unit", "sqrt(2 E_bound (1 + 1/m_min)) + max abs(u) of mass-0 bodies: bounds every relative "
        "speed (amendment c, sec. 4.7.14)", None, "-", "scale")
    row("c_tol", "edmd_gen3.c `tol_update`: K 2 d v_ref u_t", c_tol, "px^2 (c = r^2 - d^2)", "a contact found with c in [-c_tol, 0] is within the time "
        "rounding (counted, executed); below it an overlap repair. c error <= 2 d x (gap error), gap error <= 2.5 v_ref u_t",
        (2 * D * gp) if gp is not None else None, f"2 d x the largest pair contact gap ({s('gap_pair')}); margin: the smallest over the runs of "
        f"K v_ref / (gap / u_t) at each run's own v_ref ({s('m_pair')})")
    rows[-1] = rows[-1][:7] + (g("m_pair"), rows[-1][8])
    row("tol_face", "edmd_gen3.c `tol_update`: K v_ref u_t + 8 ulp(boxW)", tol_face, "px", "a wall or body face gap in [-tol_face, 0] at a "
        "contact is within rounding; below it wall_overdue / obj_overlap_repair",
        max(x for x in (g("gap_wall"), g("gap_div")) if x is not None) if (g("gap_wall") is not None or g("gap_div") is not None) else None,
        f"the largest wall / divider contact gap ({s('gap_wall')}; {s('gap_div')}); margin: the smallest over the runs of K v_ref / (gap / u_t) "
        f"at each run's own v_ref, walls and divider (the 8 ulp(box) term only adds) ({s('m_wall')}; {s('m_div')})")
    rows[-1] = rows[-1][:7] + (min(x for x in (g("m_wall"), g("m_div")) if x is not None) if (g("m_wall") is not None or g("m_div") is not None) else None, rows[-1][8])
    row("tol_pair (validator)", "edmd_gen3.c `S->tol_pair = fmax(1e-7, 1e-6 d)`", tol_pair, "px", "experiment_validation.c's pair tolerance, so "
        "a local or full check finding means the same as the driver's validator", g("local_worst_neg"), f"-local_worst ({s('local_worst_neg')})")
    row("tol_wall (validator)", "edmd_gen3.c `S->tol_wall = fmax(1e-6, 1e-6 max(1, R))`", tol_wall, "px", "experiment_validation.c's wall "
        "tolerance; also the body checks of full_check", g("full_worst_neg"), f"-full_worst ({s('full_worst_neg')})")
    row("tol_cell", "edmd_gen3.c `S->tol_cell = 1e-9`", tol_cell, "px", "a local coordinate outside [-tol_cell, w + tol_cell] is a cell repair "
        "(bookkeeping failure); crossings leave a residual of a few ulp(w)", g("cross_residual"), f"cross_residual_max ({s('cross_residual')})")
    row("band margin", "edmd_gen3.c `S->band_margin = S->tol_wall`", tol_wall, "px", "the closed band test's margin; needs >= the membership error "
        "(tol_cell) + rounding of [lo, hi] and h (a few ulp(box)) + v_ref u_t, about 1.2e-9 px (the no-miss argument, sec. 4.7.14)",
        g("cross_residual"), f"the membership error measured as cross_residual_max ({s('cross_residual')})")
    row("check interval", "edmd_gen3.c `S->check_interval = 24.0`", 24.0, "units (1 sigma-time)", "the full overlap check's cadence (sec. 4.7.1 h)",
        None, "-", "cadence")
    row("same-time limit", "edmd_gen3.c `same_limit = max(5000, 4 N)`", 5000.0, "events at one time", "more events at one time than this stops "
        "the run (stagnation): a perfect lattice has at most ~4 N simultaneous events", g("stagnation"), f"stagnation stops ({s('stagnation')})", "cut-off")
    row("heap compaction", "edmd_gen3.c `compact_at = 64 N + 4096`", None, "heap entries", "performance only (stale entries removed); no effect "
        "on the schedule (the hash is unchanged by it)", None, "-", "performance")
    row("pair rule", "edmd_gen3.c `pair_rule`", None, "-", "b >= 0 or disc <= 0: no event; dt < 0 -> 0. NO time cut-off (gen2's `t <= 1e-12` is "
        "not in gen3)", None, "-", "rule")
    row("body rule (spring)", "edmd_gen3.c `body_rule`", None, "-", "the first downward zero of the gap, bracketed between the closed-form zeros of "
        "g', safeguarded Newton; the slow-approach O(1) jump with rounding guards (j -+ 1); no tolerance other than the bracket "
        "(amendment b: 0 failures in 8 x 2000 constructed cases)", None, "-", "rule")
    row("audit abs(dt)", "edmd_gen3.c `audit_cmp`: 1e-9 abs, 1e-10 of the horizon", 1e-9, "units", "diagnostic only (the schedule audit); not "
        "read by the dynamics", None, "-", "diagnostic")
    row("validator cadence", "00ALLINONE.c `G3_VALIDATOR_EVERY = 60`", 60.0, "steps (1 sigma-time)", "the driver's validator every 60 steps under "
        "gen3 (--validator-every); gen2 every step", None, "-", "cadence")
    row("psi6 sample time", "00ALLINONE.c `g3_psi6_series_sample`: t + 1e-9 < next", 1e-9, "sigma-time", "a reader's comparison; does not steer",
        None, "-", "reader")
    row("M4 fixed point", "00ALLINONE.c `g3_lattice_fit`: abs(g_new - g) <= 1e-15 (1 + g), 200 iterations", 1e-15, "relative", "the wall clearance "
        "g = s(g)/2 (a contraction); a positive gap is required after it", None, "-", "rule")
    print("| quantity | where | value | scale | derivation / role | largest measured today | source | margin (value / measured) | kind |\n"
          "|---|---|---|---|---|---|---|---|---|")
    for name, where, value, scale, der, meas, msrc, margin, kind in rows:
        vs = "-" if value is None else (f"{value:.4g}" if isinstance(value, float) else str(value))
        ms = "-" if meas is None else f"{meas:.3g}"
        print(f"| {name} | {where} | {vs} | {scale} | {der} | {ms} | {msrc} | {'-' if margin is None else f'{margin:.3g}'} | {kind} |")
    print("\nCounters that must be 0 (largest over every record read): " + ", ".join(
        f"{k} {g(k):g}" for k in ("overlap_repair", "wall_overdue", "obj_overlap_repair", "cell_repair", "local_findings", "full_findings",
                                   "body_findings", "stagnation") if k in M))


if __name__ == "__main__":
    main()
