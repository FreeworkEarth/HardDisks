#!/usr/bin/env python3
"""Compare my per-cell cost numbers with the rounded values printed in CC's geometry_cost.md (all 14 cells). Read-only."""
import os, re
import pandas as pd
OUT = os.path.dirname(os.path.abspath(__file__))
MD = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_gen3_d57_261016/geometry_cost/geometry_cost.md"
P = pd.read_csv(os.path.join(OUT, "my_cost_percell.csv"), dtype={"pin": str})
txt = open(MD).read().split("## The cost")[1]
rows = [l for l in txt.splitlines() if re.match(r"^\| (900|1600) \|", l)]
nbad = 0
for l in rows:
    f = [x.strip() for x in l.strip("|").split("|")]
    N, eta, pin = int(f[0]), float(f[1]), f[2]
    zl, zu = map(float, f[3].split("/")); runs = int(f[4])
    e0, e1 = map(float, f[5].split(".."))
    med, mx = map(float, re.match(r"([\d.]+) \(([\d.]+)\)", f[6]).groups())
    v3, v2 = map(float, f[7].split("/")); chl, chu = map(float, f[8].split("/"))
    r = P[(P.N == N) & ((P.eta - eta).abs() < 1e-9) & (P.pin == pin)].iloc[0]
    mine = [round(r.Zl, 2), round(r.Zu, 2), r.runs, round(r.ev_min), round(r.ev_max), round(r.mix_med, 2), round(r.mix_max, 2), round(r.v3max, 2), round(r.v2max, 2), round(r.ch_lower, 1), round(r.ch_upper, 1)]
    cc = [zl, zu, runs, e0, e1, med, mx, v3, v2, chl, chu]
    ok = all(abs(a - b) < 1e-9 for a, b in zip(mine, cc))
    nbad += (not ok)
    print(f"{N} {eta:.2f} {pin:6s} {'AGREE' if ok else 'DIFFER'}  cc={cc}  mine={mine}")
print("cells:", len(rows), " differing:", nbad)
# per-model table
txt2 = open(MD).read().split("The minutes per run per CPU model")[1]
for l in [l for l in txt2.splitlines() if re.match(r"^\| (EPYC|Gold|E5)", l)]:
    print("  CC:", l)
