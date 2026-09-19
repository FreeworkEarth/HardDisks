#!/usr/bin/env python3
"""##CHRIS: run the damped-cosine frequency fit over every eligible route-A trace
and write nu_fitted + sigma_nu next to the binned nu. Read-only on the traces."""
import csv, glob, os, re, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit_nu_damped import fit_trace
C = sys.argv[1]; MAN = sys.argv[2]; OUT = sys.argv[3]
NW = int(sys.argv[4]) if len(sys.argv) > 4 else 3

man = list(csv.DictReader(open(MAN)))
idx = {}
for leaf in glob.glob(f"{C}/eta_*"):
    for f in glob.glob(f"{leaf}/wall_x_positions_*.csv"):
        m = re.search(r"wallmassfactor_(\d+)_run(\d+)\.csv$", f)
        if m: idx[(os.path.basename(leaf), int(m.group(1)), int(m.group(2)))] = f

def leaf_of(eta):
    s = f"{float(eta):.6f}".rstrip("0").rstrip(".")
    return "eta_" + s.replace(".", "p")

jobs = []
for i, r in enumerate(man):
    key = None
    for lk in {k[0] for k in idx}:
        if abs(float(lk.replace("eta_", "").replace("p", ".")) - float(r["eta"])) < 5e-4:
            key = (lk, int(r["M"]), int(r["run"])); break
    if key and key in idx: jobs.append((i, idx[key], float(r["nu"]) if r["nu"] not in ("", "nan") else None))

def work(a):
    i, path, hint = a
    try: res = fit_trace(path, hint)
    except Exception: res = None
    return (i, res)

if __name__ == "__main__":
    print(f"traces matched: {len(jobs)} of {len(man)}", flush=True)
    with Pool(NW) as p:
        out = p.map(work, jobs, chunksize=25)
    got = 0
    for r in man: r["nu_fitted"] = ""; r["sigma_nu_fitted"] = ""; r["gamma_fitted"] = ""; r["fit_rms"] = ""
    for i, res in out:
        if res:
            man[i]["nu_fitted"] = f"{res['nu']:.10g}"; man[i]["sigma_nu_fitted"] = f"{res['sigma_nu']:.6g}"
            man[i]["gamma_fitted"] = f"{res['gamma']:.6g}"; man[i]["fit_rms"] = f"{res['rms']:.6g}"; got += 1
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(man[0].keys())); w.writeheader(); w.writerows(man)
    print(f"fitted {got}/{len(jobs)}  -> {OUT}", flush=True)
