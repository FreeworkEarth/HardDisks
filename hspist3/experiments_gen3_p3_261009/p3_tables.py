#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.31; plan-author decision 13, item 5): the tables of PILOT P3, two starts -- a pilot, no paper
use, no verdict. Printed per cell and start (A = the M4 lattice at the target density, B = M4b, the compressed fluid), over the 6 seeds:
  - psi6 (global) and Z_global in the four quarters of the 2e4-sigma-time hold (quarter q = [5000 (q - 1), 5000 q) sigma-time): per seed
    the mean of the psi6(t) rows of the hold in the quarter, and the time average of Z_global over the quarter's 100 virial blocks; then
    the mean and SE (SD / sqrt(n)) over seeds;
  - A - B in the last quarter with its SE (sqrt(SE_A^2 + SE_B^2)), for psi6 and for Z, and the difference in units of its SE;
  - the per-seed psi6 at the end of the hold (its last psi6(t) row);
  - events per second (the run record's events / run_s, the run proper); for B also the compression's own wall time and its psi6 at the
    end of the compression ([EDMD3-M4B] lines);
  - the label per cell (information, the plan author's): "converged" if A - B is within 2 SE in both psi6 and Z in the last quarter.
The inventory (every run present, exit 0, one record clean=1, one build) is printed first; a missing or unclean run is listed and left
out of the means.
usage (from hspist3/ of the engine-gen3 worktree): python3 experiments_gen3_p3_261009/p3_tables.py --data <data root>
"""
import argparse, csv, glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p3_run as R


def mse(v):
    v = [x for x in v if x is not None and math.isfinite(x)]
    if len(v) < 2: return (v[0] if v else float("nan")), float("nan"), len(v)
    m = sum(v) / len(v); return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / math.sqrt(len(v)), len(v)


def run(d):
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    ec = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
    recs = re.findall(r"^\[EDMD3-HEALTH\] (.*)$", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
    kv = dict(re.findall(r" (\w+)=(\S+)", " " + recs[0].split(": ", 1)[-1])) if len(recs) == 1 else {}
    ok = ec == "0" and nb == 1 and kv.get("clean") == "1"
    v = dict(ok=ok, why="" if ok else f"exit {ec}, builds {nb}, records {len(recs)}, clean {kv.get('clean')}", build=(kv.get("build"), kv.get("target")))
    if not ok: return v
    ev = sum(int(kv[k]) for k in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band"))
    v["evs"] = ev / float(kv["run_s"])
    ps = [r for r in csv.DictReader(open(glob.glob(os.path.join(d, "psi6_t_*_run0.csv"))[0])) if r["phase"] == "hold"]
    vb = [r for r in csv.DictReader(open(glob.glob(os.path.join(d, "virial_blocks_*_run0.csv"))[0])) if r["phase"] == "hold"]
    v["psi_q"] = []; v["z_q"] = []
    for q in range(4):
        x = [float(r["psi6_global"]) for r in ps if 5000.0 * q <= float(r["t_sigma"]) < 5000.0 * (q + 1)]
        v["psi_q"].append(sum(x) / len(x) if x else float("nan"))
        b = [r for r in vb if 100 * q <= int(r["block"]) < 100 * (q + 1)]
        dt = [float(r["dt_sigma"]) for r in b]
        v["z_q"].append(sum(float(r["Z_global"]) * w for r, w in zip(b, dt)) / sum(dt) if len(b) == 100 else float("nan"))
    v["psi_end"] = float(ps[-1]["psi6_global"]); v["t_end"] = float(ps[-1]["t_sigma"])
    m = re.search(r"^\[EDMD3-M4B\] .*; ([0-9.]+) s$", lg, re.M); v["comp_s"] = float(m.group(1)) if m else None
    m = re.search(r"^\[EDMD3-M4B\] psi6 at the end of the compression \(the carried state\): global (\S+)", lg, re.M)
    v["psi_comp"] = float(m.group(1)) if m else None
    return v


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); a = ap.parse_args(); O = os.path.abspath(a.data)
    print("# Pilot P3, two starts (261012 sec. 4.7.31; decision 13, item 5), printed by experiments_gen3_p3_261009/p3_tables.py -- A PILOT, no paper use\n")
    print("A = the M4 lattice at the target density; B = M4b (the lattice at eta_0 = 0.60 in a longer box, 500 sigma-time held, the pistons at "
          "0.01 sigma per sigma-time, KE rescaled); then the divider held 2e4 sigma-time; 6 seeds per start and cell; SE = SD / sqrt(n) over seeds\n")
    V, builds = {}, set()
    print("## Inventory\n\n| cell | start | runs | clean | findings |\n|---|---|---|---|---|")
    for c in R.cells():
        for lab, _ in R.STARTS:
            rs = []
            for r in range(R.NSEED):
                v = run(os.path.join(O, c["name"], lab, f"seed{r}")); rs.append(v)
                if v["ok"]: builds.add(v["build"])
            V[(c["name"], lab)] = rs
            bad = [f"seed{r}: {v['why']}" for r, v in enumerate(rs) if not v["ok"]]
            print(f"| {c['name']} | {lab} | {R.NSEED} | {sum(v['ok'] for v in rs)} | {'; '.join(bad) if bad else '-'} |")
    print(f"\nbuilds: {sorted(builds)}")
    print("\n## psi6 (global) and Z_global in the four quarters of the hold (mean (SE) over seeds)\n")
    print("| cell | start | n | psi6 Q1 | psi6 Q2 | psi6 Q3 | psi6 Q4 | Z Q1 | Z Q2 | Z Q3 | Z Q4 | events/s (run) | B: compression wall time [s]; psi6 at its end (mean) |\n" + "|---" * 13 + "|")
    for c in R.cells():
        for lab, _ in R.STARTS:
            ok = [v for v in V[(c["name"], lab)] if v["ok"]]
            ps = [mse([v["psi_q"][q] for v in ok]) for q in range(4)]; zs = [mse([v["z_q"][q] for v in ok]) for q in range(4)]
            ev = mse([v["evs"] for v in ok])[0]
            comp = "-" if lab == "A" else f"{mse([v['comp_s'] for v in ok])[0]:.1f}; {mse([v['psi_comp'] for v in ok])[0]:.4f}"
            print(f"| {c['name']} | {lab} | {len(ok)} | " + " | ".join(f"{m:.4f} ({s:.4f})" for m, s, _ in ps) + " | "
                  + " | ".join(f"{m:.4f} ({s:.4f})" for m, s, _ in zs) + f" | {ev:.3g} | {comp} |")
    print("\n## A - B in the last quarter (Q4), and the label (information: 'converged' if within 2 SE in both psi6 and Z)\n")
    print("| cell | psi6 A - B (SE) | / SE | Z A - B (SE) | / SE | label |\n|---|---|---|---|---|---|")
    for c in R.cells():
        okA = [v for v in V[(c["name"], "A")] if v["ok"]]; okB = [v for v in V[(c["name"], "B")] if v["ok"]]
        pa, sa, _ = mse([v["psi_q"][3] for v in okA]); pb, sb, _ = mse([v["psi_q"][3] for v in okB])
        za, sza, _ = mse([v["z_q"][3] for v in okA]); zb, szb, _ = mse([v["z_q"][3] for v in okB])
        dp, sp = pa - pb, math.hypot(sa, sb); dz, sz = za - zb, math.hypot(sza, szb)
        lab = "converged" if abs(dp) < 2 * sp and abs(dz) < 2 * sz else "not converged"
        print(f"| {c['name']} | {dp:+.4f} ({sp:.4f}) | {dp / sp:+.2f} | {dz:+.4f} ({sz:.4f}) | {dz / sz:+.2f} | {lab} |")
    print("\n## The per-seed psi6 at the end of the hold (its last psi6(t) row)\n")
    print("| cell | start | seed 0 | seed 1 | seed 2 | seed 3 | seed 4 | seed 5 |\n|---|---|---|---|---|---|---|---|")
    for c in R.cells():
        for lab, _ in R.STARTS:
            print(f"| {c['name']} | {lab} | " + " | ".join(f"{v['psi_end']:.4f}" if v["ok"] else "-" for v in V[(c["name"], lab)]) + " |")
    print("\n(PILOT: no verdict, no paper use. The hold's last psi6(t) row is at t = "
          f"{next((v['t_end'] for vs in V.values() for v in vs if v['ok']), float('nan')):.3f} sigma-time.)")


if __name__ == "__main__":
    main()
