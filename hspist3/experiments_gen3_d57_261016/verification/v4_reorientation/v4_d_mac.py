#!/usr/bin/env python3
"""Independent verification of Report D57-3 task 4 (d): the Mac run (read-only on the repo)."""
import glob, gzip, hashlib, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/gen3_stage2_261016")
import p2_lattice as LT

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
MAC = os.path.join(HS, "experiments_gen3_d57_261016/data/stage2/gen3/P2_N0900_eta0.720000_G1_pin0/m_4500")
CCF = os.path.join(HS, "experiments_gen3_d57_261016/reorientation/reorientation_mac_frames.csv")
KOA = os.path.join(HS, "experiments_gen3_stage2_koa_261016/P2/merged/stage2/gen3/P2_N0900_eta0.720000_G1_pin0/m_4500")
OUT = os.path.dirname(os.path.abspath(__file__))

# checksums of the uncompressed files against .sha256_uncompressed
want = dict(reversed(l.split()) for l in open(os.path.join(MAC, ".sha256_uncompressed")).read().strip().split("\n"))
for name, h in want.items():
    p = os.path.join(MAC, name + ".gz")
    got = hashlib.sha256(gzip.open(p, "rb").read()).hexdigest()
    print(f"sha256 {name}: {'OK' if got == h else 'MISMATCH ' + got}")

snp = glob.glob(os.path.join(MAC, "snapshots_L0_*_run0.bin.gz")); assert len(snp) == 1
S = LT.read_snapshots(snp[0])
print(f"snapshots: N {S['N']}, W {S['W']:.6f}, H {S['H']:.6f}, R {S['R']}, every {S['every']}, frames {len(S['frames'])}, "
      f"t first {S['frames'][0]['t']:.4f} last {S['frames'][-1]['t']:.4f}")
rows = []
for fr in S["frames"]:
    if fr["t"] > 600.5: break
    res = LT.config_strain(fr["xy"], fr["cen"], S["W"], S["H"])
    row = dict(t=fr["t"], phase=fr["phase"], cen=fr["cen"])
    for k, s in ((0, "left"), (1, "right")):
        th = res[k]["theta0"]
        row.update({f"theta0_{s}": th, f"line_{s}": res[k]["line"], f"ordered_bulk_{s}": res[k]["ordered_bulk"], f"e_{s}": res[k]["e"],
                    f"ytheta_{s}": float(abs(np.mod(th, 60.0) - 30.0) <= 15.0)})
    rows.append(row)
F = pd.DataFrame(rows); F.to_csv(os.path.join(OUT, "v4_d_mac_frames.csv"), index=False, float_format="%.9g")
dts = np.diff(F["t"].to_numpy()); print(f"frames to 600: {len(F)}, dt min {dts.min():.5f} max {dts.max():.5f}; phases {F['phase'].unique()}")
for s in ("left", "right"):
    agree = (F[f"line_{s}"] == 90).to_numpy() == (F[f"ytheta_{s}"] == 1).to_numpy()
    y = F[F[f"line_{s}"] == 90]["t"]
    t1 = y.iloc[0]
    back = F[(F["t"] > t1) & (F[f"line_{s}"] == 0)]
    nans = F[f"line_{s}"].isna().sum()
    print(f"{s}: first config_strain line-y snapshot t = {t1:.5f}; later x snapshots up to 600: {len(back)}; NaN lines {nans}; "
          f"config_strain line vs theta0-mod-60 rule agree on {agree.sum()} of {len(F)} frames")

ps = pd.read_csv(glob.glob(os.path.join(MAC, "psi6_t_L0_*_run0.csv.gz"))[0])
print(f"psi6_t: {len(ps)} samples, dt median {np.median(np.diff(ps['t_sigma'])):.5f}")
for k, s in ((0, "left"), (1, "right")):
    o = np.mod(ps[f"psi6_orient_deg_{k}"].to_numpy(float) + 15.0, 60.0) - 15.0
    i = np.where(np.abs(o - 30.0) <= 7.5)[0][0]
    print(f"{s}: first psi6_t sample within 7.5 deg of 30: t = {ps['t_sigma'].iloc[i]:.6f} (orientation {o[i]:.3f}); previous sample t {ps['t_sigma'].iloc[i-1]:.4f} o {o[i-1]:.3f}")

# against CC's frame CSV
C = pd.read_csv(CCF)
print(f"\nCC frames rows {len(C)}, mine {len(F)}; t max |diff| {np.abs(C['t'].to_numpy() - F['t'].to_numpy()).max():.3e}")
for col in C.columns:
    if col in ("t",): continue
    a, b = F[col].to_numpy(float), C[col].to_numpy(float)
    if col.startswith("line") or col == "phase":
        print(f"  {col}: equal {np.array_equal(a, b)}")
    else:
        rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-12)
        print(f"  {col}: max abs diff {np.abs(a - b).max():.3e}, max rel diff {rel.max():.3e}")

# run.log: pin and seed
L = open(os.path.join(MAC, "run.log")).read()
print("\nrun.log [EDMD3-M4C]:", [l for l in L.split("\n") if "[EDMD3-M4C]" in l])
print("run.log seeds:", sorted(set(re.findall(r"seed[ =]+(\d+)", L))))
R = pd.read_csv(os.path.join(KOA, "p2_runs.csv")); r0 = R[R["run"] == 0].iloc[0]
print(f"KOA P2 (0.72, G1, pin0, M 4500, run 0): seed {r0['seed']}, m4c '{r0['m4c']}', arrangement {r0['arrangement']}, snapshots {r0['snapshots']}")
kl = open(os.path.join(KOA, "run.log")).read()
print("KOA run.log seeds (run 0 lines):", sorted(set(re.findall(r"run=0 seed=(\d+)", kl))))
