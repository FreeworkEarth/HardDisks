#!/usr/bin/env python3
"""Sensitivity (context only, not a recomputation of CC's numbers): core-hours with per-pin per-model medians; and with the 0.78 G1 max overshoot."""
import os
import numpy as np, pandas as pd
OUT = os.path.dirname(os.path.abspath(__file__))
C = pd.read_csv(os.path.join(OUT, "my_cost_runs.csv"), dtype={"pin": str})
J = pd.read_csv(os.path.join(OUT, "my_p2_join.csv"))
J["pin"] = J.cell.str.extract(r"_pin(\w+)$")[0]
W = {"EPYC 7543P": 32, "EPYC 9654": 59, "E5-2670 v3": 125, "E5-4640 v2": 10, "Gold 6248": 133, "Gold 6342": 41}
medp = J.groupby(["model", "pin"]).us_per_event.median().to_dict()
def us(m, N, pin):
    v = medp[(m, pin)]
    if N == 1600: v *= (4.580 / 4.242 if m == "E5-4640 v2" else 3.836 / 3.533)
    return v
for lab in ("lower", "upper"):
    g = C[C.model == lab]
    ch = sum(10 * r["events"] * sum(W[m] / 400 * us(m, r["N"], r["pin"]) for m in W) * 1e-6 / 3600 for r in g.to_dict("records"))
    print(f"core-hours {lab} with per-pin per-model medians: {ch:.2f}")
g = C[C.model == "upper"]
base = 10 * g.coresec_mix_per_run.sum() / 3600
print(f"upper x 1.018 everywhere (P2's largest executed/upper ratio): {base * 1.018:.2f}")
