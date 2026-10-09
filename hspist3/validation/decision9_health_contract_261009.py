#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.13, plan-author decision 9 of 2026-10-09, read-only): are the speed-of-sound runs with a
health line (sec. 4.7.10) excluded by the loader's health contract, and do the paper tables show the reduced n?
  1. every run with an [EDMD-HEALTH] line in A1 v2 (2), A2_topup (1) and campaign_r25 (504): its cell, mass, run index, the
     non-zero counters, and the loader's verdict for it (tests_20260913.cell_runs: discarded or used);
  2. the paper tables that hold those cells: n used / discarded (A1 v2: 260919_A1v2_final_cs_vs_eta.csv; A2: n_runs of
     260919_A2_cs_per_mass.csv);
  3. what the paper scripts read from campaign_r25 (the guard trace of sec. 4.7.10, paper_inputs_by_campaign.txt and the
     per-path lists): which of its files, and whether any of the 504 runs is among them.
usage (from hspist3/): python3 validation/decision9_health_contract_261009.py
"""
import collections, csv, glob, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); REPO = os.path.dirname(HS)
sys.path.insert(0, HERE)
import tests_20260913 as T
SOS = os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN")
FINAL = os.path.join(REPO, "0000_PLAN_OVERALL/paper1_speedofsound/experiments/final")
LINE = re.compile(r"\[EDMD-HEALTH\] L0=([\d.]+) M=(\d+) run=(\d+) seed=(\d+): forced_advance=(\d+) wall_clamp_repairs=(\d+) "
                  r"overlap_repairs=(\d+) wall_overdue=(\d+)")


def health_runs(campaign):
    out = []
    for lg in sorted(glob.glob(os.path.join(SOS, campaign, "**", "run.log"), recursive=True)):
        if "/analysis/" in lg: continue
        for m in LINE.finditer(open(lg, errors="ignore").read()):
            out.append(dict(cell=os.path.dirname(lg), L0=float(m.group(1)), M=int(m.group(2)), r=int(m.group(3)), seed=int(m.group(4)),
                            counters={k: int(m.group(i)) for i, k in ((5, "forced_advance"), (6, "wall_clamp_repairs"),
                                                                       (7, "overlap_repairs"), (8, "wall_overdue")) if int(m.group(i))}))
    return out


def verdict(h):
    runs = T.cell_runs(h["cell"], h["M"])
    d = {r: disc for r, p, disc in runs}
    if h["r"] not in d: return "no trace of this run index in the cell"
    return "DISCARDED" if d[h["r"]] else "USED"


def main():
    print("# Decision 9: the loader's health contract and the paper tables (261012 sec. 4.7.13), printed by "
          "validation/decision9_health_contract_261009.py\n")
    print("## 1. The runs with a health line, and the loader's verdict (tests_20260913.cell_runs)\n")
    print("| campaign | cell | M | run | non-zero counters | loader |\n|---|---|---|---|---|---|")
    tally = collections.Counter()
    for camp in ("A1v2_20260914", "A2_topup_20260912"):
        for h in health_runs(camp):
            v = verdict(h); tally[(camp, v)] += 1
            print(f"| {camp} | {os.path.relpath(h['cell'], os.path.join(SOS, camp))} | {h['M']} | {h['r']} | "
                  f"{', '.join(f'{k}={c}' for k, c in h['counters'].items())} | {v} |")
    r25 = health_runs("campaign_r25_psi6_20260823")
    byv = collections.Counter(); bycell = collections.Counter()
    for h in r25:
        v = verdict(h); byv[v] += 1; bycell[(os.path.basename(h["cell"]), h["L0"], v)] += 1
    print(f"| campaign_r25_psi6_20260823 | {len(r25)} runs, by cell (below) | all | - | "
          f"{', '.join(sorted({k for h in r25 for k in h['counters']}))} | " + ", ".join(f"{v} {n}" for v, n in sorted(byv.items())) + " |")
    print("\ncampaign_r25 runs with a health line, by cell and loader verdict: " +
          "; ".join(f"{c} (L0 {L0:g}): {v} {n}" for (c, L0, v), n in sorted(bycell.items())))
    print("\n## 2. The paper tables\n")
    rows = {r["eta"]: r for r in csv.DictReader(open(os.path.join(FINAL, "260919_A1v2_final_cs_vs_eta.csv")))}
    print("A1 v2 (260919_A1v2_final_cs_vs_eta.csv, 9 masses x 25 trajectories per cell):")
    for e in ("0.009817", "0.006545"):
        r = rows.get(e)
        print(f"  eta {e}: trajectories_used {r['trajectories_used']}, trajectories_discarded {r['trajectories_discarded']}" if r else f"  eta {e}: not in the table")
    a2 = [r for r in csv.DictReader(open(os.path.join(FINAL, "260919_A2_cs_per_mass.csv")))
          if abs(float(r["eta"]) - 0.10) < 0.001 and r["N"] == "1600"]
    print("A2 (260919_A2_cs_per_mass.csv), eta 0.10, N = 1600: " + ", ".join(f"M {r['M']}: n_runs {r['n_runs']}" for r in a2))
    print("\n## 3. What the paper scripts read from campaign_r25 (the guard trace of sec. 4.7.10)\n")
    tr = os.path.join(HS, "experiments_gen2_defect_scan_261009", "paper_inputs_by_campaign.txt")
    for line in open(tr):
        if "campaign_r25" in line: print("  " + line.strip())
    read = set()
    SPT = sys.argv[1] if len(sys.argv) > 1 else None   # optional: the scratch folder with the per-script path lists
    if SPT and os.path.isdir(SPT):
        for f in glob.glob(os.path.join(SPT, "paths_*.txt")):
            for p in open(f):
                if "campaign_r25" in p: read.add(os.path.relpath(p.strip(), HS))
        kinds = collections.Counter(re.sub(r"L0_\d+|wallmassfactor_\d+|run\d+", "*", os.path.basename(p)) for p in read)
        print("  files read, by kind: " + ", ".join(f"{k} x{n}" for k, n in sorted(kinds.items())))
        hit = [h for h in r25 if any(os.path.dirname(os.path.join(HS, p)) == h["cell"] and f"wallmassfactor_{h['M']}_run{h['r']}.csv" in p for p in read)]
        print(f"  of the {len(r25)} campaign_r25 runs with a health line, traces among the files read: {len(hit)}"
              + (" (" + ", ".join(f"{os.path.basename(h['cell'])} M {h['M']} run {h['r']}" for h in hit) + ")" if hit else ""))
    print("  the paper scripts use campaign_r25 for metadata only: tests_20260913.a1_leaf_table reads the L0 strings and the first "
          "trace row (eta, Predicted_Frequency) of one M = 50 run per cell; no paper number is computed from its trajectories.")


if __name__ == "__main__":
    main()
