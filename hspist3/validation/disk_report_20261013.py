#!/usr/bin/env python3
"""##CHRIS 2026-10-13: disk report -- the 20 largest campaign-level directories in the repo. READ-ONLY:
nothing is deleted, moved or copied. Writes hspist3/validation/261013_disk_report.txt and prints the same.

Campaign level = the children of the grouping directories below (a parent and its child are never both listed);
directories that are not split further are listed whole. "Summarised" = the directory's name appears in a
committed .md or .json (git grep, fixed string); "read by" = it appears in a committed .py or .sh.
The last column only marks candidates; Chris decides what is copied to an external drive.
"""
import os, re, subprocess, sys
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "hspist3", "validation", "261013_disk_report.txt")
SOS = "hspist3/experiments_speed_of_sound"
GROUPS = [SOS, f"{SOS}/EDMD", f"{SOS}/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN", f"{SOS}/EDMD/mode1_normalized_units",
          f"{SOS}/mode1_normalized_units", f"{SOS}/mode0_real_units", "hspist3/experiments_energy_transfer",
          "hspist3/experiments_szilard", "hspist3"]
WHOLE = [f"{SOS}/EDMD/mode0_real_units", "IB_Package", ".git"]
LABELS = [(r"experiments_szilard", "Szilard engine (not Paper 1/2)"), (r"A1v2_", "Paper 1 A1 v2 -- canonical c_s(eta), DROOT"), (r"A2_|famB|A3_", "Paper 1 A2/A3 finite-size ladders"),
          (r"campaign_r25|campaign_transition", "Paper 1 r25 campaign (Aug, superseded by A1 v2)"),
          (r"ladder_N|overnight_N|finitesize_aspect|routeA|routeB", "Paper 1 Aug size/aspect ladders"),
          (r"tests_2026|_orchestration|nofuse|validate_acc|debug_seed", "Paper 1 tests / checks"),
          (r"simulation_|number_and_L0|radius_sweep|eta_split|oscillation_design", "Paper 1 early sweeps (2026-02 to 08)"),
          (r"adiabatic_piston", "Paper 1/2 adiabatic piston check"), (r"^level0|level1|level2", "Paper 2 Levels 0-2"),
          (r"^level3", "Paper 2 Level 3"), (r"^level4", "Paper 2 Level 4"), (r"^level5", "Paper 2 efficiency map"),
          (r"experiments_szilard", "Szilard engine (not Paper 1/2)"), (r"mode0_real_units", "early real-units runs"),
          (r"SINGLE_analysis", "early single-run analyses"), (r"IB_Package", "IB_Package (not an experiment)"),
          (r"^\.git$", "git object store")]

def sh(*a): return subprocess.run(a, cwd=REPO, capture_output=True, text=True).stdout

def sizes():
    out, groups = {}, set(GROUPS)
    for g in GROUPS:
        for line in sh("du", "-k", "-d", "1", g).splitlines():
            kb, p = line.split("\t", 1)
            if p == g or p in groups or any(p == w or p.startswith(w + "/") for w in WHOLE): continue
            if any(o != g and p.startswith(o + "/") for o in GROUPS if o.startswith(g + "/")): continue
            if any(o.startswith(p + "/") for o in GROUPS + WHOLE): continue      # p is an ancestor of a listed directory
            out[p] = int(kb)
    for w in WHOLE:
        out[w] = int(sh("du", "-sk", w).split()[0])
    return out

def label(p):
    b = os.path.basename(p)
    for rx, lab in LABELS:
        if re.search(rx, p) if rx.startswith("experiments_") else (re.search(rx, b) or re.search(rx, p)): return lab
    return "unclassified"

def grep(name, globs):
    r = sh("git", "grep", "-l", "-F", name, "--", *globs).split()
    return r

def main():
    S = sizes(); top = sorted(S.items(), key=lambda kv: -kv[1])[:20]
    df = sh("df", "-k", ".").splitlines()[-1].split()
    log = os.path.join(REPO, "hspist3/experiments_energy_transfer/level5_effmap_20261010/RUN.log")
    prev = [l for l in open(log).read().splitlines() if "/System/Volumes/Data" in l] if os.path.exists(log) else []
    total = int(sh("du", "-sk", ".").split()[0])
    lines = ["Disk report 2026-10-13 (read-only; nothing deleted, moved or copied)", "",
             f"Data volume now: {int(df[1])/1048576:.0f} GiB total, {int(df[3])/1048576:.1f} GiB free ({df[4]} used)",
             f"Free space at the end of the efficiency map (RUN.log df line): {prev[-1].split()[3] if prev else 'n/a'}",
             f"Repository total: {total/1048576:.1f} GiB", "",
             "| # | path | GiB | experiment / level | summarised in committed md/json | read by committed scripts | archive candidate (Chris decides) |",
             "|---|---|---|---|---|---|---|"]
    for i, (p, kb) in enumerate(top, 1):
        name = os.path.basename(p) if p not in WHOLE else p
        md = grep(name, ["*.md", "*.json"]) if p != ".git" else []
        sc = grep(name, ["*.py", "*.sh"]) if p != ".git" else []
        if p == ".git": cand = "no (history)"
        elif md and not sc: cand = "YES -- summarised, no committed script reads it"
        elif md and sc: cand = f"with care -- summarised, but read by {len(sc)} script(s)"
        else: cand = "no -- not summarised in any committed md/json"
        lines.append(f"| {i} | `{p}` | {kb/1048576:.2f} | {label(p)} | {'yes (' + str(len(md)) + ' files, e.g. ' + md[0] + ')' if md else 'no'} | "
                     f"{', '.join(os.path.basename(x) for x in sc[:3]) + (' ...' if len(sc) > 3 else '') if sc else 'none'} | {cand} |")
    s = sum(kb for _, kb in top)
    lines += ["", f"Top 20 together: {s/1048576:.1f} GiB of the repository's {total/1048576:.1f} GiB."]
    txt = "\n".join(lines) + "\n"; open(OUT, "w").write(txt); print(txt)

if __name__ == "__main__":
    main()
