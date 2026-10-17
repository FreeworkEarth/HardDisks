import math, os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
RU = pd.read_csv("my_runs.csv"); E1 = pd.read_csv(HS + "/experiments_gen3_p2_261016/report_261016/p2_cells_e1.csv").set_index("cell")
def cell(c):
    g = RU[RU.cell == c]; return dict(e=g.e_rec.mean(), se=g.e_rec.std(ddof=1)/math.sqrt(len(g)), c=E1.loc[c,"E1_c"], sc=E1.loc[c,"E1_sig"])
a, b = cell("P2_N0900_eta0.720000_G1_pin0"), cell("P2_N0900_eta0.720000_G1_pinturned")
d = a["e"] - b["e"]; sd = math.sqrt(a["se"]**2 + b["se"]**2)
print(f"G1 pair: e pin0 {100*a['e']:+.4f}%, turned {100*b['e']:+.4f}%, diff {100*d:.4f}% +- {100*sd:.4f}% = {d/sd:.2f} SE; c diff {a['c']-b['c']:+.4f} +- {math.hypot(a['sc'],b['sc']):.4f}")
def fit(e, se, y, sy):
    bb = 0.0
    for _ in range(200):
        w = 1/(sy**2 + bb**2*se**2); S, Sx, Sy, Sxx, Sxy = w.sum(), (w*e).sum(), (w*y).sum(), (w*e*e).sum(), (w*e*y).sum()
        D = S*Sxx - Sx**2; bn = (S*Sxy - Sx*Sy)/D; aa = (Sxx*Sy - Sx*Sxy)/D
        if abs(bn-bb) < 1e-13: bb = bn; break
        bb = bn
    w = 1/(sy**2 + bb**2*se**2); r = w*(y-aa-bb*e)**2
    return aa, bb, math.sqrt(S/D), float(r.sum()), r
w = np.array([1/a["sc"]**2, 1/b["sc"]**2]); cP1 = float((w*np.array([a["c"], b["c"]])).sum()/w.sum()); sP1 = w.sum()**-0.5
p2, p3 = cell("P2_N0900_eta0.720000_G2_pinturned"), cell("P2_N0900_eta0.720000_G2_pin0")
for lab, seP1 in (("SE_e(P1) as specified", math.sqrt(a["se"]**2+b["se"]**2)/2), ("SE_e(P1) = half the G1 difference", abs(d)/2),
                  ("SE_e(P1) = sd of the 100 G1 runs / sqrt(100)", None)):
    if seP1 is None:
        g = RU[RU.cell.isin(["P2_N0900_eta0.720000_G1_pin0","P2_N0900_eta0.720000_G1_pinturned"])]; seP1 = g.e_rec.std(ddof=1)/10
    e = np.array([(a["e"]+b["e"])/2, p2["e"], p3["e"]]); se = np.array([seP1, p2["se"], p3["se"]])
    y = np.log([cP1, p2["c"], p3["c"]]); sy = np.array([sP1/cP1, p2["sc"]/p2["c"], p3["sc"]/p3["c"]])
    aa, bb, sb, chi2, r = fit(e, se, y, sy)
    print(f"3-point, {lab} = {100*seP1:.4f}%: b {bb:.4f} +- {sb:.4f}, chi2 {chi2:.2f} (per point {np.round(r,2).tolist()})")
# four-cell fit chi2 per point
cs = [cell(c) for c in ["P2_N0900_eta0.720000_G1_pin0","P2_N0900_eta0.720000_G1_pinturned","P2_N0900_eta0.720000_G2_pin0","P2_N0900_eta0.720000_G2_pinturned"]]
e = np.array([q["e"] for q in cs]); se = np.array([q["se"] for q in cs]); y = np.log([q["c"] for q in cs]); sy = np.array([q["sc"]/q["c"] for q in cs])
aa, bb, sb, chi2, r = fit(e, se, y, sy); print(f"4-cell fit: b {bb:.6f} chi2 {chi2:.4f}, per point (G1 pin0, G1 turned, G2 pin0, G2 turned) {np.round(r,2).tolist()}")
# the pilot's two runs' own c (CC's md: 16.006 and 15.913)
cr = pd.read_csv(HS + "/experiments_gen3_d56_261015/c_analysis_b/c_runs.csv"); cr = cr[cr.cell == 76]
cc = pd.read_csv(HS + "/experiments_gen3_d56_261015/c_analysis_b/c_cells.csv"); x = float(cc[cc.cell == 76].x.iloc[0])
print("pilot runs own c:", [round(math.sqrt(n**2 + (1/(2*math.pi*t))**2)/x, 4) for n, t in zip(cr.nu_d, cr.tau_r)], "Z mean", cr.Z.mean(),
      "psi6", ((cr.psi6_rec_0 + cr.psi6_rec_1)/2).mean())
