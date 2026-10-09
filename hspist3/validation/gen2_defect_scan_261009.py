#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.10; plan-author decision 5 of 2026-10-09): READ-ONLY provenance scan for the gen2
pass-through defect. In edmd.c's collide_time_ab, `if(t<=1e-12) return 0;` drops a pair contact that is due exactly now (c = 0)
or within 1e-12 time units, so a disk whose velocity changes at the moment it touches a second disk can pass into it (261012
sec. 4.7.6, the cradle cells). gen2 is not patched (frozen builds). This script writes nothing.

1. BUILDS. Every `<= 1e-12` cut-off in edmd.c, quoted with line number and function, plus the pair safety net, for: the builds the
   plan author named (279282b, 73fc07f, 7b08827), the only two commits that ever changed collide_time_ab (3bb8c42, 20df8c1), and
   every build recorded in the scanned data. For each build: does 00ALLINONE.c print [EDMD-HEALTH] in speed-of-sound mode, in
   energy-transfer mode, and the contact audit line ([EDMD-CONTACT], HD_CONTACT_AUDIT only)?
2. CAMPAIGNS: the papers' data roots (the root list of experiments_loader_guard_261008/probe_data_trees.py, sec. 4.7.5), each
   marked if the 20 paper scripts read it directly (guard trace: experiments_gen2_defect_scan_261009/paper_inputs_by_campaign.txt),
   plus Test T and T-prime. Per campaign: runs; runs with a health record; runs with a health line; runs with overlap_repairs > 0
   and the largest count; validator findings by reason; runs WITHOUT a health record (a coverage gap: not counted as clean);
   collisions (measured by the contact audit where it ran, otherwise ESTIMATED); expected hits at 1e-12 per collision.

HEALTH RECORD (the definition used here): the run's stdout section is in a log; for speed-of-sound the section is complete
("Finished valid run"); and the build prints [EDMD-HEALTH] in that mode whenever a counter is non-zero (only then does a missing
line mean "all counters 0"). The build is the recorded one where the data record it (.build_git, the summaries' build_git column,
"00ALLINONE git" lines). Speed-of-sound campaigns made on the Mac record no build: for them the print is established BY DATE (the
line exists since the 2026-08-23 fix, 01_improvements_bugsfxed_dev/26_08_23_EDMD_PP_TUNNELLING_FIX_AND_PHASE_AWARE_EOS.md, first
committed in 20df8c1) and BY OBSERVATION (a campaign with at least one health line shows its build printed it); the basis is
printed per campaign. Energy-transfer runs have a health record only on a build with the energy-transfer print (9cafd7f and later
on engine-divider-resched); no paper campaign ran on one.
VALIDATOR: both modes run the fail-closed overlap validator (experiment_validation.c check_overlaps) at every recorded step.
Findings: speed-of-sound "INVALID speed-of-sound run" lines, the "N invalid" of the batch lines, speed_of_sound_failures.csv rows;
energy-transfer *.failures.csv rows and "ABORTING INVALID RUN" lines.
COLLISION ESTIMATE [DERIVATION]: pair collisions = N Gamma t / 2, Gamma = 4 (Z - 1) / sqrt(pi) per disk and sigma-time (2D
Enskog collision rate with the contact value from Z = 1 + 2 eta g(sigma); kT = m = sigma = 1), Z from Liu's global equation of
state (plot_speed_of_sound_edmd.Z_liu_global), t = hold + record per run. Walls, divider and pistons are not included (the audited
Test T and T-prime count them; their measured/estimated ratio is printed as the check).
usage (from hspist3/): python3 validation/gen2_defect_scan_261009.py
"""
import collections, csv, glob, math, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); REPO = os.path.dirname(HS)
sys.path.insert(0, HS)
import numpy as np
from plot_speed_of_sound_edmd import Z_liu_global

RATE = 1e-12                                     # the plan author's hit rate per collision (261012 sec. 4.7.6, report 2)
SOS = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN"
GATE2 = "experiments_resched_gate2_261007/" + SOS
TRACE = os.path.join(HS, "experiments_gen2_defect_scan_261009", "paper_inputs_by_campaign.txt")
NAMED = [("3bb8c42", "first commit of edmd.c"), ("20df8c1", "safety net and health line first committed"),
         ("279282b", "named: KOA ~/harddisks"), ("73fc07f", "named: KOA ~/harddisks_resched"), ("7b08827", "named: Test T and T-prime")]
FIX_DATE = "2026-08-23"


def git_show(rev, path):
    r = subprocess.run(["git", "-C", REPO, "show", f"{rev}:{path}"], capture_output=True, text=True, errors="replace")
    return r.stdout if r.returncode == 0 else None


def git_desc(rev):
    r = subprocess.run(["git", "-C", REPO, "log", "-1", "--format=%h %ad", "--date=format:%Y-%m-%d %H:%M", rev],
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else "(not in this repository)"


FUNC = re.compile(r"^static\s+[A-Za-z_][\w\s\*]*?\b([A-Za-z_]\w*)\s*\(")
CUT = re.compile(r"<=\s*1e-12\)")


def cutoffs(src):
    """(line, function, code) of every code line with a `<= 1e-12)` cut-off, and the pair safety net."""
    out, fn = [], "?"
    for i, line in enumerate(src.split("\n"), 1):
        m = FUNC.match(line)
        if m: fn = m.group(1)
        s = line.strip()
        if s.startswith(("/*", "*", "//")) or "`" in s: continue
        if CUT.search(s) or (fn == "collide_time_ab" and s.startswith("if(c<0.0)")):
            out.append((i, fn, s))
    return out


def build_facts(rev):
    e, a = git_show(rev, "hspist3/edmd_core/edmd.c"), git_show(rev, "hspist3/00ALLINONE.c")
    if e is None or a is None: return None
    return dict(cut=cutoffs(e), sos="[EDMD-HEALTH] L0=" in a, et="[EDMD-HEALTH] energy-transfer" in a,
                contact="[EDMD-CONTACT] executed events" in a)


# ---------------------------------------------------------------------------------------------------------- campaigns
def campaigns():
    """(label, path rel. to hspist3, kind, paper) -- the root list of probe_data_trees.py (sec. 4.7.5), T and T-prime."""
    out = [(x, f"{SOS}/{x}", "sos", "Paper 1") for x in
           ("A1v2_20260914", "famB_20260911", "A2_dilute50_20260917", "A2_dilute_20260916", "A2_long200_20260915",
            "A2_topup_20260912", "A2_alpha2_20260912", "campaign_r25_psi6_20260823", "routeA_lowdensity_20260912",
            "confinement_B_20261013", "confinement_pilot_20261013", "tests_20260913")]
    for x in sorted(os.listdir(os.path.join(HS, "experiments_energy_transfer"))):
        p = os.path.join(HS, "experiments_energy_transfer", x)
        if os.path.isdir(p) and x.startswith(("paper1_confinement_", "level")):
            out.append((x, f"experiments_energy_transfer/{x}", "et", "Paper 1 (confinement A)" if x.startswith("paper1") else "Paper 2"))
    for x in ("experiments_speed_of_sound", "experiments_energy_transfer"):
        if os.path.isdir(os.path.join(HS, "experiments_resched_gate_261005", x)):
            out.append((f"resched_gate_261005/{x.split('_', 1)[1]}", f"experiments_resched_gate_261005/{x}", "sos" if "sound" in x else "et", "gate (sec. 4.4)"))
    out += [("Test T (resched_testT_261007)", f"{GATE2}/resched_testT_261007", "sos", "gate (sec. 4.4.12)"),
            ("Test T-prime (resched_testTprime_261007)", f"{GATE2}/resched_testTprime_261007", "sos", "gate (sec. 4.4.15)")]
    return out


P1_SCRIPTS = {"populate", "figures", "damping", "massladder", "a2boxtrunc", "boxtrunc", "boxtrunc_tab", "conf_results", "conf_afix",
              "conf_heldwall", "conf_prereg", "canonical", "melting", "draft_audit", "roman"}
P2_SCRIPTS = {"p2_figures", "p2_geomfix", "p2_rampfast", "p2_level2Au"}


def traced():
    """campaign path -> the papers whose scripts read it ("Paper 1", "Paper 2", "gate script only")."""
    out = {}
    for line in open(TRACE):
        if line.startswith("#"): continue
        f = [x.strip() for x in line.split(" | ")]
        sc = set(f[3].split(","))
        who = [p for p, S in (("Paper 1", P1_SCRIPTS), ("Paper 2", P2_SCRIPTS)) if sc & S]
        out[f[0]] = ", ".join(who) if who else "gate script only"
    return out


# ---------------------------------------------------------------------------------------------------- speed-of-sound
RUN_HDR = re.compile(r"^##RUN (\d{4}-\d\d-\d\d)")
RUNNING = re.compile(r"Running: L0 = [-\d.]+, M = (\d+)\*m, run = (\d+), seed = (\d+)")
TARGET = re.compile(r"target=\d+ cycles\s+f_pred=\S+\s+samples=\d+\s+T=([-\d.eE+]+) sigma-time")
RELEASED = re.compile(r"Wall released \(EDMD\) at t = ([-\d.eE+]+)")
HEALTH_SOS = re.compile(r"\[EDMD-HEALTH\] L0=[-\d.]+ M=(\d+) run=(\d+) seed=(\d+): forced_advance=(\d+) wall_clamp_repairs=(\d+) "
                        r"overlap_repairs=(\d+) wall_overdue=(\d+)(?: past_events=(\d+))?")
CONTACT = re.compile(r"\[EDMD-CONTACT\] executed events (\d+)")
FINISHED = re.compile(r"Finished valid run (\d+)")
INVALID_SOS = re.compile(r"INVALID speed-of-sound run L0=\S+ M=(\d+) repeat=(\d+) \[([^\]]*)\]")
BATCH = re.compile(r"Speed-of-sound batch complete: (\d+) valid, (\d+) invalid, (\d+) requested")
KE_N = re.compile(r"\[KE-AUDIT\].*?\sN=(\d+)\s")
BUILD_LINE = re.compile(r"00ALLINONE\s+git (\S+)\s+target (\S+)")
TRACE_NAME = re.compile(r"^wall_x_positions_.*_run(\d+)\.csv$")
COUNTERS = ("forced_advance", "wall_clamp_repairs", "overlap_repairs", "wall_overdue", "past_events")


def eta_from_path(rel):
    for part in reversed(rel.split("/")):
        m = re.match(r"^eta_(\d+)p(\d+)$", part) or re.match(r"^e(\d+)p(\d+)_", part)
        if m: return float(f"{m.group(1)}.{m.group(2)}")
        if part.startswith(("epi8_", "koa_pi8_", "mac_pi8_")) or part == "pi8": return math.pi / 8
    return None


def first_row(path):
    try:
        with open(path, newline="") as f:
            r = csv.reader(f); h = next(r); v = next(r)
        return dict(zip(h, v))
    except Exception:
        return None


def scan_sos(root):
    """runs keyed by (directory, seed), from the logs; traces keyed the same way; validator findings; builds."""
    runs, unmatched_health, invalid, builds = {}, [], collections.Counter(), collections.Counter()
    batch_invalid = 0
    logs = [p for p in glob.glob(os.path.join(root, "**", "*.log"), recursive=True) if "/analysis/" not in p]
    logs += [p for p in glob.glob(os.path.join(root, "**", ".failed_run*", "*.log"), recursive=True)]
    for lg in sorted(set(logs)):
        d = os.path.dirname(lg); cur = None; date = None; lastN = None
        failed_dir = "/.failed_run" in lg
        for line in open(lg, errors="ignore"):
            m = RUN_HDR.match(line)
            if m: date = m.group(1); continue
            m = KE_N.search(line)
            if m: lastN = int(m.group(1)); continue
            m = RUNNING.search(line)
            if m:
                M, r, seed = int(m.group(1)), int(m.group(2)), int(m.group(3))
                key = (d.split("/.failed_run")[0], seed)
                cur = runs.get(key)
                if cur is None or (not cur["finished"] and not failed_dir):
                    cur = dict(dir=key[0], M=M, r=r, seed=seed, date=date, N=lastN, T=None, hold=None, contact=None,
                               finished=False, health=None, log=lg, failed_dir=failed_dir)
                    runs[key] = cur
                else:
                    cur = None                           # a second copy of a finished run: keep the first
                continue
            m = HEALTH_SOS.search(line)
            if m:
                h = dict(zip(COUNTERS, [int(x) if x is not None else 0 for x in m.groups()[3:]]))
                M, r, seed = int(m.group(1)), int(m.group(2)), int(m.group(3))
                if cur is not None and cur["seed"] == seed: cur["health"] = h
                else: unmatched_health.append((d, M, r, seed, h))
                continue
            m = BATCH.search(line)
            if m: batch_invalid += int(m.group(2)); continue
            m = INVALID_SOS.search(line)
            if m: invalid[m.group(3)] += 1; continue
            m = BUILD_LINE.search(line)
            if m: builds[m.group(1)] += 1; continue
            if cur is None: continue
            m = TARGET.search(line)
            if m: cur["T"] = float(m.group(1)); continue
            m = RELEASED.search(line)
            if m: cur["hold"] = float(m.group(1)); continue
            m = CONTACT.search(line)
            if m: cur["contact"] = int(m.group(1)); continue
            m = FINISHED.search(line)
            if m: cur["finished"] = True; continue
    for d, M, r, seed, h in unmatched_health:            # a health line printed outside its run's section: attach by seed
        if (d, seed) in runs: runs[(d, seed)]["health"] = h
    for f in glob.glob(os.path.join(root, "**", "speed_of_sound_failures.csv"), recursive=True):
        for row in csv.DictReader(open(f, newline="", errors="ignore")):
            if row.get("reason") and row["reason"] != "reason": invalid["ledger: " + row["reason"]] += 1
    for f in glob.glob(os.path.join(root, "**", ".build_git"), recursive=True):
        m = BUILD_LINE.search(open(f, errors="ignore").read())
        if m: builds[m.group(1)] += 1
    # traces on disk (Mac campaigns): the run universe beyond the logs (matched by seed), and N and eta per directory
    cell, ntr, tr_nolog = {}, 0, 0
    for f in glob.glob(os.path.join(root, "**", "wall_x_positions_*_run*.csv"), recursive=True):
        if "/analysis/" in f or "invalid" in os.path.basename(f).lower(): continue
        d = os.path.dirname(f); r0 = first_row(f) or {}; ntr += 1
        if d not in cell:
            try: cell[d] = dict(N=int(r0["Left_Count"]) + int(r0["Right_Count"]), eta=float(r0["eta"]))
            except Exception: cell[d] = dict(N=None, eta=None)
        try: seed = int(r0["Seed"])
        except Exception: seed = None
        if seed is None or (d, seed) not in runs: tr_nolog += 1
    return runs, invalid, batch_invalid, builds, cell, ntr, tr_nolog


# --------------------------------------------------------------------------------------------------- energy transfer
HEALTH_ET = re.compile(r"\[EDMD-HEALTH\] energy-transfer seed=(\d+): forced_advance=(\d+) wall_clamp_repairs=(\d+) "
                       r"overlap_repairs=(\d+) wall_overdue=(\d+)(?: past_events=(\d+))?")
ABORT = re.compile(r"ABORTING INVALID RUN \[([^\]]*)\]")
INIT = re.compile(r"\[EDMD\] initialized: N=(\d+)")


def fnum(v):
    try: return float(v)
    except Exception: return None


def scan_et(root):
    runs, invalid, builds, health = [], collections.Counter(), collections.Counter(), {}
    started = 0
    for f in glob.glob(os.path.join(root, "**", "*.csv"), recursive=True):
        b = os.path.basename(f)
        if b.endswith(".failures.csv"):
            for row in csv.DictReader(open(f, newline="", errors="ignore")):
                if row.get("reason") and row["reason"] != "reason": invalid["ledger: " + row["reason"]] += 1
            continue
        if not b.startswith("summary"): continue
        for row in csv.DictReader(open(f, newline="", errors="ignore")):
            if row.get("mode") not in ("energy_transfer", None) or row.get("timestamp") == "timestamp": continue
            bg = (row.get("build_git") or "").strip() or None
            builds[bg or "(not recorded)"] += 1
            N = fnum(row.get("particles_total")); eta = fnum(row.get("eta_particles_region")) or fnum(row.get("eta_nominal"))
            dt = fnum(row.get("dt_sigma")); st = fnum(row.get("steps_after_release")); hs = fnum(row.get("wall_hold_steps")) or 0.0
            t = (st + hs) * dt if (dt and st is not None) else None
            runs.append(dict(dir=os.path.dirname(f), seed=row.get("seed"), build=bg, N=N, eta=eta, t=t))
    for lg in glob.glob(os.path.join(root, "**", "*.log"), recursive=True):
        for line in open(lg, errors="ignore"):
            if INIT.search(line): started += 1; continue
            m = HEALTH_ET.search(line)
            if m: health[(os.path.dirname(lg), m.group(1))] = dict(zip(COUNTERS, [int(x) if x else 0 for x in m.groups()[1:]])); continue
            m = ABORT.search(line)
            if m: invalid["abort: " + m.group(1)] += 1
    return runs, invalid, builds, health, started


# ----------------------------------------------------------------------------------------------------------- helpers
def gamma(eta):
    """2D Enskog collision rate per disk and sigma-time, kT = m = sigma = 1: 4 (Z - 1) / sqrt(pi), Z from Liu's global EOS."""
    Z = float(np.asarray(Z_liu_global(np.array([eta])))[0])
    return 4.0 * (Z - 1.0) / math.sqrt(math.pi)


def est_pairs(N, eta, t):
    if not (N and eta and t) or eta <= 0: return None
    return N * gamma(eta) * t / 2.0


def g(x, fmt="{:.2e}"):
    return "-" if x is None else fmt.format(x)


def main():
    tr = traced(); facts = {}; data_builds = collections.Counter()
    rows, checks, details, kinds = [], [], [], collections.defaultdict(collections.Counter)
    for label, rel, kind, paper in campaigns():
        root = os.path.join(HS, rel)
        if not os.path.isdir(root): rows.append((label, paper, "missing", *["-"] * 11)); continue
        hits = sorted({v for t, v in tr.items() if t == rel or t.startswith(rel + "/")})
        read = "; ".join(hits) if hits else "no"
        if kind == "sos":
            runs, invalid, batch_inv, builds, cell, ntr, tr_nolog = scan_sos(root)
            for b, n in builds.items(): data_builds[b.split("-")[0]] += n
            nruns = len(runs) + tr_nolog
            hl = [r for r in runs.values() if r["health"]]
            ov = [r["health"]["overlap_repairs"] for r in hl if r["health"]["overlap_repairs"] > 0]
            other = sum(1 for r in hl if any(r["health"][k] for k in COUNTERS if k != "overlap_repairs"))
            for r in hl:
                for k in COUNTERS:
                    if r["health"][k]: kinds[label][k] += 1
                if r["health"]["overlap_repairs"] > 0: details.append((label, r))
            # the health print: the recorded build where there is one, else the date of the runs and the lines seen
            recb = {b.split("-")[0] for b in builds}
            for b in recb:
                if b not in facts: facts[b] = build_facts(b)
            if recb:
                prints = all(facts.get(b) and facts[b]["sos"] for b in recb); basis = "recorded build " + ",".join(sorted(recb))
            else:
                dates = sorted({r["date"] for r in runs.values() if r["date"]})
                mt = sorted({os.path.getmtime(r["log"]) for r in runs.values()})
                from datetime import datetime
                lo = dates[0] if dates else (datetime.fromtimestamp(mt[0]).strftime("%Y-%m-%d") if mt else None)
                prints = (lo is not None and lo >= FIX_DATE) or bool(hl)
                basis = (f"by date (runs from {lo})" if lo else "no date") + (", lines seen" if hl else "")
            have = [r for r in runs.values() if r["finished"] and prints]
            no_rec = nruns - len(have)
            meas = [r["contact"] for r in runs.values() if r["contact"] is not None]
            est = 0.0; n_est = 0
            for r in runs.values():
                c = cell.get(r["dir"], {})
                N = r["N"] or c.get("N"); eta = c.get("eta") or eta_from_path(os.path.relpath(r["dir"], HS))
                t = (r["T"] or 0.0) + (r["hold"] or 0.0) if r["T"] else None
                e = est_pairs(N, eta, t)
                if e: est += e; n_est += 1
                if r["contact"] is not None and e: checks.append((label, r["contact"], e))
            coll = (f"{sum(meas):.3e} measured ({len(meas)} runs)" if meas else "") + \
                   (("; " if meas else "") + f"{est:.3e} estimated ({n_est} runs)" if n_est else "")
            colln = sum(meas) if meas else est
            vf = sum(invalid.values()) + batch_inv
            vdet = ", ".join(f"{k} {v}" for k, v in invalid.most_common()) or ""
            vtxt = f"{vf}" + (f" ({vdet})" if vdet else "") + (f"; batch lines: {batch_inv} invalid" if batch_inv and not vdet else "")
            rows.append((label, paper, read, nruns, len(have), len(hl), len(ov), max(ov) if ov else 0, other, vtxt, no_rec,
                         basis, coll, colln))
        else:
            runs, invalid, builds, health, started = scan_et(root)
            for b, n in builds.items():
                if b != "(not recorded)": data_builds[b.split("-")[0]] += n
            recb = {b.split("-")[0] for b in builds if b != "(not recorded)"}
            for b in recb:
                if b not in facts: facts[b] = build_facts(b)
            et_print = {b for b in recb if facts.get(b) and facts[b]["et"]}
            have = [r for r in runs if r["build"] and r["build"].split("-")[0] in et_print]
            hl = list(health.values()); ov = [h["overlap_repairs"] for h in hl if h["overlap_repairs"] > 0]
            for h in hl:
                for k in COUNTERS:
                    if h[k]: kinds[label][k] += 1
            other = sum(1 for h in hl if any(h[k] for k in COUNTERS if k != "overlap_repairs"))
            est = sum(e for e in (est_pairs(r["N"], r["eta"], r["t"]) for r in runs) if e)
            n_est = sum(1 for r in runs if est_pairs(r["N"], r["eta"], r["t"]))
            vf = sum(invalid.values()); vdet = ", ".join(f"{k} {v}" for k, v in invalid.most_common())
            basis = ("builds " + ",".join(f"{b} x{n}" for b, n in sorted(builds.items())) if builds else "no summaries") + \
                    ("; none prints the energy-transfer health line" if not et_print else "") + f"; {started} runs started (logs)"
            rows.append((label, paper, read, len(runs), len(have), len(hl), len(ov), max(ov) if ov else 0, other,
                         f"{vf}" + (f" ({vdet})" if vdet else ""), len(runs) - len(have), basis,
                         f"{est:.3e} estimated ({n_est} runs)" if n_est else "-", est))
    # ------------------------------------------------------------------------------------------------ 1. builds
    print("# gen2 pass-through defect: read-only provenance scan (261012 sec. 4.7.10, decision 5)\n")
    print("## 1. The cut-offs per build (edmd.c, quoted), and what each build's 00ALLINONE.c prints\n")
    allb = [(b, why) for b, why in NAMED] + [(b, "recorded in the data") for b in sorted(data_builds) if b not in dict(NAMED)]
    for b, why in allb:
        f = facts.get(b) or build_facts(b)
        if f is None: print(f"### {b} ({why}): not in this repository\n"); continue
        print(f"### {b} ({why}; {git_desc(b)})  health line: speed-of-sound {'yes' if f['sos'] else 'NO'}, "
              f"energy-transfer {'yes' if f['et'] else 'NO'}; contact audit {'yes' if f['contact'] else 'no'}\n")
        print("```")
        for i, fn, s in f["cut"]: print(f"edmd.c:{i:<5d} {fn:24s} {s}")
        print("```\n")
    print("recorded builds in the scanned data (runs or records): " +
          ", ".join(f"{b} x{n}" for b, n in sorted(data_builds.items())) + "\n")
    # ---------------------------------------------------------------------------------------------- 2. campaigns
    print("## 2. Per campaign\n")
    print("| campaign | data root of | read directly by the scripts of | runs | with a health record | health lines | overlap_repairs > 0 "
          "| largest | other counters > 0 | validator findings | NO health record (gap) | health-print basis | collisions | "
          f"expected hits at {RATE:.0e} | observed (overlap_repairs > 0) |")
    print("|" + "---|" * 15)
    tot = collections.defaultdict(collections.Counter); totc = collections.Counter()
    for r in rows:
        if r[2] == "missing": print(f"| {r[0]} | {r[1]} | (missing) |" + " - |" * 12); continue
        label, paper, read, n, have, hl, ov, mx, other, vtxt, norec, basis, coll, colln = r
        print(f"| {label} | {paper} | {read} | {n} | {have} | {hl} | {ov} | {mx} | {other} | {vtxt} | {norec} | {basis} | {coll} | "
              f"{g(colln * RATE)} | {ov} |")
        grp = "gates (sec. 4.4)" if paper.startswith("gate") else paper
        for k, v in (("campaigns", 1), ("runs", n), ("have", have), ("norec", norec), ("hl", hl), ("ov", ov)): tot[grp][k] += v
        totc[grp] += colln
    print("\n### Totals\n")
    print("| group | campaigns | runs | with a health record | WITHOUT one (gap) | health lines | runs with overlap_repairs > 0 | "
          f"collisions (measured or estimated) | expected hits at {RATE:.0e} |\n|---|---|---|---|---|---|---|---|---|")
    for grp in sorted(tot):
        t = tot[grp]
        print(f"| {grp} | {t['campaigns']} | {t['runs']} | {t['have']} | {t['norec']} | {t['hl']} | {t['ov']} | {totc[grp]:.3e} | "
              f"{totc[grp] * RATE:.2e} |")
    print("\n### Every run with overlap_repairs > 0, and whether the paper loader uses it\n")
    for label, r in details:
        rel = os.path.relpath(r["dir"], HS)
        try:
            sys.path.insert(0, HERE); import tests_20260913 as TT
            used = {rr: d for rr, p, d in TT.cell_runs(r["dir"], r["M"])}
            verdict = ("DISCARDED by the health contract (tests_20260913.cell_runs)" if used.get(r["r"]) else
                       "USED by tests_20260913.cell_runs" if r["r"] in used else "not among the cell's traces")
        except Exception as e:
            verdict = f"loader check failed: {type(e).__name__}"
        print(f"- {label}: {rel}, M = {r['M']}, run {r['r']}, seed {r['seed']}, N = {r['N']}; record {r['T']} sigma-time; "
              + ", ".join(f"{k}={v}" for k, v in r["health"].items()) + f"; finished: {r['finished']}; {verdict}")
    if not details: print("(none)")
    print("\n### Health lines by counter (runs with that counter > 0)\n")
    for label in kinds:
        print(f"- {label}: " + ", ".join(f"{k} {v}" for k, v in kinds[label].most_common()))
    if not kinds: print("(none)")
    print("\n## 3. Check of the collision estimate on the audited runs (executed events, all kinds / estimated pair collisions)\n")
    for lab in sorted({c[0] for c in checks}):
        q = np.array([m / e for l, m, e in checks if l == lab])
        print(f"{lab}: {len(q)} runs, ratio median {np.median(q):.3f}, range {q.min():.3f}-{q.max():.3f}")


if __name__ == "__main__":
    main()
