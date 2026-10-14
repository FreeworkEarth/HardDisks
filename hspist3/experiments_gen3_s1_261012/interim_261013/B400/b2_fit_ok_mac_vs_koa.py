import sys, os, math, time, warnings
sys.path.insert(0, "validation"); sys.path.insert(0, ".")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import resched_testT_followup_261007 as FU
d = sys.argv[1]
red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
z = np.load(os.path.join(d, "acf_runs.npz")); nl = min(len(z[f"run{r}"]) for r in red["run"])
dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
koa = pd.read_csv(os.path.join(d, "nud_runs.csv")).set_index("run")
for r in red["run"]:
    if int(koa.loc[r, "fit_ok"]) == 0: continue
    t0 = time.time()
    try:
        FU.fit_acf(np.asarray(z[f"run{r}"][:nl], dtype=float), dt, nu_ref); ok = 1; msg = ""
    except Exception as e:
        ok = 0; msg = f"{type(e).__name__}: {str(e).strip()}"
    if not ok: print(f"{os.path.basename(os.path.dirname(d))}/{os.path.basename(d)} run {int(r)} seed {int(koa.loc[r, 'seed'])}: KOA fit_ok 1 (nu_d {koa.loc[r, 'nu_d']:.10g}), the Mac raises {msg} ({time.time() - t0:.0f} s)", flush=True)
print(f"{d}: done", flush=True)
