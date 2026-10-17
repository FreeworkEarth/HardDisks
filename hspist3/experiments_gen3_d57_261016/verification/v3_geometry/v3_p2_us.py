#!/usr/bin/env python3
"""Independent check: P2's measured per-CPU-model median microseconds per executed event. Read-only on the repo."""
import glob, os, re
import numpy as np, pandas as pd

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
OUT = os.path.dirname(os.path.abspath(__file__))
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016/P2")
TASKS = os.path.join(HS, "cluster/gen3_stage2_261016/tasks_p2_koa_261016.txt")

lines = [l.split() for l in open(TASKS) if l.strip()]
print("task lines:", len(lines))
pat = re.compile(r"^chunk (\d+) on (\S+) \((.*)\): group (\S+), lines (\d+)-(\d+);")
recs = []
for f in sorted(glob.glob(os.path.join(P2, "job_logs/gen3-p2_15083717_*.out"))):
    aid = int(re.search(r"_(\d+)\.out$", f).group(1))
    hits = [pat.match(l) for l in open(f, errors="replace")]
    hits = [h for h in hits if h]
    assert len(hits) == 1, (f, len(hits))
    h = hits[0]
    j, node, cpu, lo, hi = int(h.group(1)), h.group(2), h.group(3), int(h.group(5)), int(h.group(6))
    assert j == aid and lo == hi == j, (f, j, lo, hi)
    done = any(l.strip() == f"chunk {j} done" for l in open(f, errors="replace"))
    t = lines[j - 1]
    recs.append(dict(chunk=j, node=node, cpu=cpu, done=done, folder=t[1], M=int(t[2]), run=int(t[3]), seed=int(t[4])))
L = pd.DataFrame(recs)
def short(c):
    m = re.search(r"(EPYC \w+)", c) or re.search(r"(E5-\d+ v\d)", c) or re.search(r"(Gold \d+)", c)
    return m.group(1)
L["model"] = L.cpu.map(short)
L["cell"] = L.folder.str.split("/").str[2]
print("logs:", len(L), " all done:", bool(L.done.all()), " unique (cell,M,run):", L[["cell", "M", "run"]].drop_duplicates().shape[0])
print(L.model.value_counts().to_dict())

R = pd.concat([pd.read_csv(f, dtype={"pin": str}) for f in sorted(glob.glob(os.path.join(P2, "merged/stage2/gen3/P2_*/m_*/p2_runs.csv")))], ignore_index=True)
print("p2_runs rows:", len(R), " errors non-empty:", int(R["error"].notna().sum()))
J = R.merge(L, on=["cell", "M", "run"], how="outer", indicator=True, suffixes=("", "_log"))
print("join:", J["_merge"].value_counts().to_dict())
print("node column agrees with log node:", bool((J.node == J.node_log).all()), " seed agrees:", bool((J.seed == J.seed_log).all()))
J["us_per_event"] = J.seconds / J.events * 1e6
CC = {"EPYC 9654": 2.139, "Gold 6342": 2.209, "EPYC 7543P": 2.284, "Gold 6248": 2.764, "E5-2670 v3": 3.555, "E5-4640 v2": 3.935}
print("\nmodel        n   median us/event (mine)   CC value   |diff|   (min, max)")
out = []
for m, g in J.groupby("model"):
    med = float(np.median(g.us_per_event))
    out.append(dict(model=m, n=len(g), median=med, cc=CC[m], diff=med - CC[m], lo=g.us_per_event.min(), hi=g.us_per_event.max()))
    print(f"  {m:11s} {len(g):4d}   {med:.6f}   {CC[m]:.3f}   {abs(med - CC[m]):.2e}   ({g.us_per_event.min():.3f}, {g.us_per_event.max():.3f})")
pd.DataFrame(out).to_csv(os.path.join(OUT, "my_p2_us_per_model.csv"), index=False, float_format="%.8g")
# by pin, in case pin matters for the per-event cost
print("\nper model and pin: median us/event")
print(J.groupby(["model", "pin"]).us_per_event.median().unstack().round(4).to_string())
print("\nP2 total core-hours used (sum of seconds / 3600):", round(float(J.seconds.sum() / 3600), 4))
J[["cell", "M", "run", "seed", "node", "model", "seconds", "events", "us_per_event"]].to_csv(os.path.join(OUT, "my_p2_join.csv"), index=False)
