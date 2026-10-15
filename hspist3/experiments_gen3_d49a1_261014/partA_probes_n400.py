#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 49 addendum 1 item 1; 261012 sec. 4.7.111): code-fact probes of Part A's held path at N = 400,
the way Report S1-0's AF1-AF6 probed it (gate G0 f; experiments_gen3_s10_261012/s10_runs.py, af_case): the worker's AF command
(cluster/confinement_20261013/conf_worker.sh, mode AF) under --engine=gen3, on the Mac, with the clean Mac build of f72d764
(experiments_gen3_d18_261010/bin_f72d764; its SHA-256 is checked against its build record first). Test seeds, a few hundred held steps.
NOTHING HERE IS PART A DATA: no stencil, no record; the facts are which cells initialise, where the divider sits, which start the held
path uses at N = 400, and the event log's rate and divider share.
  Q1  eta 0.30, L0 = N pi / (8 H eta) exactly (off the 1/48 grid), divider at L0                  -> expected: the box is cut, the run stops
  Q2  eta pi/8, L0 = 20 (on the grid), divider at L0 + sigma_x/2 (off the 1/24 grid)              -> where the divider sits at t = 0
  Q3  L0 = 546/48 (eta 0.69046, next to the 0.69 cell), 3000 held steps, without --quiet       -> the seeder's lines; rows per sigma-time
  Q4  L0 = 7540/48 (eta 0.049999, next to the 0.05 cell), 3000 held steps                        -> rows per sigma-time at the dilute end
usage (from hspist3/): python3 experiments_gen3_d49a1_261014/partA_probes_n400.py
"""
import hashlib, math, os, re, subprocess, sys, time
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
BIN = os.path.join(HS, "experiments_gen3_d18_261010", "bin_f72d764", "00ALLINONE")
OUT = os.path.join(HERE, "probes")
ENV = dict(os.environ, SDL_VIDEODRIVER="dummy")
N, H, NS = 400, 20.0, 200


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def case(name, L0, xw, hold=600, post=60, quiet=True):
    d = os.path.join(OUT, name); seed = 7740100 + int(name[1])
    argv = [BIN, "--mode=edmd", "--experiment=energy_transfer", "--headless"] + (["--quiet"] if quiet else []) + [
            "--edmd-acc=0", "--seed-drift-order=drift-first",
            f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv", "--trace-every=60", f"--particles={N}",
            f"--particles-boxes={NS},{NS}", "--particle-radius=0.5", f"--l0={L0:.17g}", f"--height={H:g}", "--num-walls=1",
            f"--wall-positions={xw:.17g}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
            "--eff-output=wall-ke", f"--wall-hold-steps={hold}", f"--steps={post}", "--fixed-dt=0.4", "--kbt1", f"--seed={seed}",
            "--engine=gen3"]
    return name, d, argv, dict(ENV, HD_PISTON_EVENTS=os.path.join(d, "ev.csv")), dict(L0=L0, xw=xw, hold=hold, post=post)


def main():
    rec = open(os.path.join(os.path.dirname(BIN), "build_record.txt")).read()
    s = sha(BIN); ok = s in rec
    print(f"# Part A code-fact probes at N = 400 (partA_probes_n400.py); binary {os.path.relpath(BIN, HS)} SHA-256 {s[:16]}...: "
          f"{'in its build record' if ok else 'NOT in its build record'}\n")
    if not ok: sys.exit("STOP: the binary is not the recorded clean build")
    import importlib.util
    spec = importlib.util.spec_from_file_location("prereg", os.path.join(HS, "validation", "paper1_confinement_prereg_20261012.py"))
    P = importlib.util.module_from_spec(spec); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS); spec.loader.exec_module(P)
    L0_30 = N * math.pi / (8 * H * 0.30); L0_p8 = N * math.pi / (8 * H * (math.pi / 8))
    e8 = P.eta_of(NS, H, L0_p8); sx8 = math.sqrt(P.leff(L0_p8) ** 2 / (2 * NS * P.cs2(e8)))
    cases = [case("Q1", L0_30, L0_30), case("Q2", L0_p8, L0_p8 + sx8 / 2),
             case("Q3", 546 / 48, 546 / 48, hold=3000, quiet=False), case("Q4", 7540 / 48, 7540 / 48, hold=3000)]
    print("| probe | L0 | eta exact | divider requested | exit | seconds | the lines that decide |\n|---|---|---|---|---|---|---|")
    facts = []
    for name, d, argv, env, meta in cases:
        if os.path.exists(d): sys.exit(f"STOP: {d} exists (a new name; nothing is deleted)")
        os.makedirs(d)
        with open(os.path.join(d, "command.txt"), "w") as fh:
            fh.write(" ".join(f"{k}={v}" for k, v in env.items() if k.startswith("HD_")) + " " + " ".join(argv) + "\n")
        t0 = time.time(); p = subprocess.run(argv, env=env, cwd=d, capture_output=True, text=True); dt = time.time() - t0
        open(os.path.join(d, "stdout.txt"), "w").write(p.stdout); open(os.path.join(d, "stderr.txt"), "w").write(p.stderr)
        lines = [l for l in (p.stdout + p.stderr).splitlines()
                 if re.search(r"WARNING: L_0|ABORTING|STOP|EDMD3-HEALTH\] energy|\[SEED\]|Particles placed|fallback", l)]
        print(f"| {name} | {meta['L0']:.9f} | {P.eta_of(NS, H, meta['L0']):.6f} | {meta['xw']:.9f} | {p.returncode} | {dt:.1f} | "
              + "; ".join(l.strip()[:220] for l in lines[:4]).replace("|", "/") + " |")
        facts.append((name, d, meta, p.returncode))
    print()
    for name, d, meta, rc in facts:
        tr = os.path.join(d, "tr.csv"); ev = os.path.join(d, "ev.csv")
        if rc == 0 and os.path.exists(tr):
            t = pd.read_csv(tr, low_memory=False); col = [c for c in t.columns if c.startswith("W0_x")]
            if col: print(f"{name}: the divider at the first trace row (after the hold) at {t[col[0]].iloc[0]:.9f} sigma; requested {meta['xw']:.9f} "
                          f"(difference {t[col[0]].iloc[0] - meta['xw']:+.2e})")
        if rc == 0 and os.path.exists(ev):
            e = pd.read_csv(ev, usecols=["t_sigma", "kind"]); held = e[e["t_sigma"] < meta["hold"] * 0.4 / 24.0]
            T = meta["hold"] * 0.4 / 24.0; nd = int(held["kind"].str.startswith("D").sum())
            print(f"{name}: event log over the {T:.2f} held sigma-time: {len(held)} rows, {len(held) / T:.1f} rows per sigma-time, "
                  f"divider rows {100 * nd / max(1, len(held)):.1f} %, {os.path.getsize(ev) / max(1, len(e)):.1f} bytes per row; kinds "
                  + ", ".join(f"{k} {v}" for k, v in held["kind"].value_counts().items()))


if __name__ == "__main__":
    main()
