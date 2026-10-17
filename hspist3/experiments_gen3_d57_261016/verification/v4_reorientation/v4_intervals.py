import os, numpy as np, pandas as pd
P2 = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_gen3_stage2_koa_261016/P2/merged/stage2/gen3"
for cell in ["P2_N0900_eta0.720000_G1_pin0", "P2_N0900_eta0.720000_G2_pin0", "P2_N0900_eta0.720000_G1_pinturned", "P2_N0900_eta0.780000_G1_pin0"]:
    for M in [900, 1800, 4500, 9000, 18000]:
        f = pd.read_csv(os.path.join(P2, cell, f"m_{M}", "p2_frames.csv"))
        dts = np.concatenate([np.diff(g.sort_values("t")["t"].to_numpy()) for _, g in f.groupby("run")])
        print(cell, M, f"diff min {dts.min():.4f} median {np.median(dts):.4f} max {dts.max():.4f}; frames/run {f.groupby('run').size().unique()}; t_last {f['t'].max():.1f}; hold frames phase0 {f[f.phase==0].groupby('run').size().unique()}")
