#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (decision 38 item 5; Report S1-0): gate G0 of the recorded sec. 4.8.10, items (a)-(h), printed from the code
(quoted lines, read from git at the named commits) and from the runs of s10_runs.py (data_f72d764/). Read-only.
  (a) the campaign commit against f72d764; the KOA binary of a643bc7 as the campaign binary?
  (b) the flags of sec. 4.8.3.2 together in the speed-of-sound loop; the combinations that never ran together
  (c) virial blocks and psi6 in the record phase as well as in the hold
  (d) nu_pred in the trace header under the exact box; valid for every cell
  (e) the nu_d fit and the averaged-ACF fit: window, start values, failure criterion
  (f) Part A: the start, the exact box and off-grid positions in the held path under gen3; the event-log size at N = 1600;
      the window parameters of reduce_AF.py
  (g) MaxArraySize on KOA
  (h) the step limits of the B command against the longest record of the grid
usage (from hspist3/): python3 experiments_gen3_s10_261012/s10_tables.py > experiments_gen3_s10_261012/s10_tables_output.txt"""
import collections, glob, gzip, math, os, re, struct, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); ROOT = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS); sys.path.insert(0, HERE)
import tests_20260913 as T
import plot_speed_of_sound_edmd as SOS
import s10_runs as RUNS
DATA = os.path.join(HERE, "data_f72d764")
MAIN = os.path.join(os.path.dirname(ROOT), "HardDisks")                      # the main worktree (Test G's Mac data)
DTF = float(np.float32(0.4)) / 24.0                                          # the binary's sample_dt_sigma: (double) 0.4f / 24


def git(*a):
    return subprocess.run(["git", "-C", ROOT] + list(a), capture_output=True, text=True, check=True).stdout


def src(commit, path):
    return git("show", f"{commit}:{path}").split("\n")


def find(lines, pat, start=0):
    for i in range(start, len(lines)):
        if pat in lines[i]: return i
    return None


def quote(label, lines, i, n=1):
    for k in range(i, i + n): print(f"    {label}:{k + 1}: {lines[k].rstrip()}")


# ---------------------------------------------------------------------------------------------------------------- (a)
def item_a():
    print("## (a) The campaign commit against f72d764; the KOA binary of a643bc7\n")
    A, F = src("a643bc7", "hspist3/00ALLINONE.c"), src("f72d764", "hspist3/00ALLINONE.c")
    print("The registered flags (sec. 4.8.3.2 and block S) in the argument parser, line of the first parse in each source:\n")
    print("| flag | a643bc7 | f72d764 |\n|---|---|---|")
    flags = ["--engine", "--gen3-seeding", "--gen3-exact-box", "--gen3-psi6-compartments", "--psi6-every", "--gen3-virial-blocks",
             "--wall-hold-steps", "--gen3-snapshots", "--gen3-lattice-orientation", "--target-oscillations", "--oscillation-safety",
             "--oscillation-min-steps", "--oscillation-max-steps", "--speed-sound-log-stride", "--speed-sound-exact-seed", "--lengths",
             "--wall-masses", "--fixed-dt"]
    miss = []
    for f in flags:
        pa = re.compile(r'(strn?cmp\(arg, *"' + re.escape(f) + r'[="]|"' + re.escape(f) + r'=?", *strlen)')
        la = next((i + 1 for i, s in enumerate(A) if pa.search(s)), None)
        lf = next((i + 1 for i, s in enumerate(F) if pa.search(s)), None)
        if la is None: miss.append(f)
        print(f"| `{f}` | {la if la else '**absent**'} | {lf if lf else '**absent**'} |")
    print(f"\nabsent in a643bc7: {', '.join(miss) if miss else 'none'}")
    i = find(A, "Unknown option '%s'.")
    if i is not None: print("an unknown option stops the a643bc7 binary:"); quote("a643bc7 00ALLINONE.c", A, i, 1)
    anc = subprocess.run(["git", "-C", ROOT, "merge-base", "--is-ancestor", "a643bc7", "f72d764"]).returncode == 0
    print(f"\na643bc7 is an ancestor of f72d764: {'yes' if anc else 'no'}")
    print("\nThe build inputs (the files the build compiles: 00ALLINONE.c, experiment_validation.c, edmd_core/*.c, their headers; the "
          "kissfft submodule commit):\n")
    print("| from | to | files that differ (git diff --numstat) |\n|---|---|---|")
    paths = ["hspist3/00ALLINONE.c", "hspist3/experiment_validation.c", "hspist3/experiment_validation.h", "hspist3/edmd_core",
             "hspist3/kissfft", "hspist3/Makefile"]
    for a, b in (("a643bc7", "f72d764"), ("f72d764", "engine-gen3")):
        ns = git("diff", "--numstat", a, b, "--", *paths).strip().split("\n")
        ns = [x for x in ns if x.strip()]
        print(f"| {a} | {b} ({git('rev-parse', '--short', b).strip()}) | " +
              ("; ".join(f"{x.split()[2]} (+{x.split()[0]} -{x.split()[1]})" for x in ns) if ns else "none") + " |")
    print("\nThe commits between a643bc7 and f72d764 that change 00ALLINONE.c:")
    for ln in git("log", "--format=%h %s", "a643bc7..f72d764", "--", "hspist3/00ALLINONE.c").strip().split("\n"):
        print(f"- {ln[:150]}")
    print("\nThe commits after f72d764 on engine-gen3 that change a build input:")
    for ln in git("log", "--format=%h %s", "f72d764..engine-gen3", "--", *paths).strip().split("\n"):
        if ln.strip(): print(f"- {ln[:150]}")
    sha = open(os.path.join(HS, "cluster", "gen3_koa_261009", "k2_binaries_a643bc7.sha256")).readline().split()[0]
    print(f"\nthe KOA binary of a643bc7 (cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256): {sha}")


# ---------------------------------------------------------------------------------------------------------------- runs
def run_log(case):
    return open(os.path.join(DATA, case, "run.log"), errors="replace").read()


def snapshots(path):
    b = open(path, "rb").read(); assert b[:8] == b"G3SNAP1\n"
    N, nd = struct.unpack_from("<ii", b, 8); off = 16 + 4 * 8
    per = collections.Counter()
    while off < len(b):
        t, ph, ndv = struct.unpack_from("<dii", b, off); off += 16 + 8 * ndv + 16 * N; per[ph] += 1
    return N, per


def item_bc():
    print("\n## (b) and (c) The flags of sec. 4.8.3.2 together, on the Mac (s10_runs.py, the Mac build of f72d764)\n")
    print("Commands: data_f72d764/<case>/command.txt (Test G's B command with the changes of sec. 4.8.3.2; HD_KE_TRACE=1 and "
          "HD_CONTACT_AUDIT=1; snapshots on).\n")
    print("| case | cell | exit | clean | exact box, eta (run record) | start (run record) | f_pred (binary / replication) | "
          "record: planned steps, sigma-time | trace rows; HD_KE_TRACE's [KE-AUDIT], [EDMD-ENERGY] lines; Box_Width_sigma | "
          "psi6 rows hold / record; per-compartment columns; "
          "spacing | virial blocks hold / record | snapshots hold / record | contact audit |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    ok_all = True
    for case in ("B1", "S0", "ST"):
        meta = eval(open(os.path.join(DATA, case, "meta.txt")).read())
        L = run_log(case); rc = open(os.path.join(DATA, case, "exit_code.txt")).read().strip()
        clean = re.search(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", L, re.M)
        rr = re.search(r"^\[EDMD3-RUN\] .*?exact_box=(\d) eta=(\S+) .*?seeding=(\S+) left=\[orient (\S+ \S+ \S+ n_s \d+ n_r \d+)", L, re.M)
        m4c = re.search(r"^\[EDMD3-M4C\] lattice orientation requested: (\S+) .*?\); (.*)$", L, re.M)
        start = f"{rr.group(3)}, " + (f"requested {m4c.group(1)}: {m4c.group(2)}" if m4c else f"default: {rr.group(4)}") if rr else "-"
        nke, nen = len(re.findall(r"^\[KE-AUDIT\] ", L, re.M)), len(re.findall(r"^\[EDMD-ENERGY\] ", L, re.M))
        fp = re.search(r"f_pred=(\S+)\s+samples=(\d+)\s+T=(\S+) sigma-time", L)
        ca = re.search(r"^\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: (.*)$", L, re.M)
        tr = glob.glob(os.path.join(DATA, case, "wall_x_positions_L0_*_run0.csv"))[0]
        hdr = open(tr).readline().strip().split(","); rows = sum(1 for _ in open(tr)) - 1
        first = open(tr).readlines()[1].strip().split(",")
        bw = float(first[hdr.index("Box_Width_sigma")]); tgt = first[hdr.index("Target_Oscillations")]
        ps = glob.glob(os.path.join(DATA, case, "psi6_t_L0_*_run0.csv"))[0]
        ph = open(ps).readline().strip().split(","); pc = collections.Counter(); tt = []
        for ln in open(ps).readlines()[1:]:
            v = ln.split(","); pc[v[1]] += 1; tt.append(float(v[0]))
        sp = np.diff(tt); comp_cols = ph[5:]                      # after t_sigma, phase, psi6_global, psi6_local, mean_neighbors
        vb = glob.glob(os.path.join(DATA, case, "virial_blocks_L0_*_run0.csv"))[0]
        vc = collections.Counter(ln.split(",")[1] for ln in open(vb).readlines()[1:])
        Nsn, sc = snapshots(glob.glob(os.path.join(DATA, case, "snapshots_L0_*_run0.bin"))[0])
        nup_rep = RUNS.nu_pred_binary(meta["N"], meta["H"], meta["L0"], meta["M"])
        steps_rep = math.ceil((200.0 / nup_rep) / DTF)
        good = (rc == "0" and clean and clean.group(1) == "1" and rr and rr.group(1) == "1" and fp and pc["hold"] > 0 and pc["record"] > 0
                and vc["hold"] > 0 and vc["record"] > 0 and sc[0] > 0 and sc[1] > 0 and ca and nke > 0 and nen > 0
                and abs(bw - 2 * meta["L0"]) < 1e-5 and tgt == "200" and comp_cols and int(fp.group(2)) == steps_rep)
        ok_all &= bool(good)
        print(f"| {case} | N {meta['N']}, eta {meta['eta']}, alpha {meta['alpha']}, L0 {meta['L0']:.17g} | {rc} | {clean.group(1) if clean else '-'} | "
              f"{rr.group(1)}, {rr.group(2)} | {start} | {fp.group(1)} / {nup_rep:.8g} | {fp.group(2)} (replication "
              f"{steps_rep}), {fp.group(3)} | {rows}; {nke}, {nen}; {bw:.6f} (2 L0 = {2 * meta['L0']:.6f}) | "
              f"{pc['hold']} / {pc['record']}; {len(comp_cols)} ({', '.join(comp_cols[:4])}{'...' if len(comp_cols) > 4 else ''}); "
              f"{sp.min():.4f}-{sp.max():.4f} (Delta_psi {meta['dpsi']:.4f}) | {vc['hold']} / {vc['record']} | {sc[0]} / {sc[1]} (N {Nsn}) | "
              f"{'yes: ' + ca.group(2) if ca else 'NO'} |")
    print(f"\n(b), (c) on the Mac: {'every check above holds' if ok_all else 'A CHECK FAILED (see the table)'}")
    print("\nThe loop's sampling calls (f72d764), the same in the hold and in the record:")
    F = src("f72d764", "hspist3/00ALLINONE.c")
    for pat in ('g3_psi6_series_sample(psi6_series, &psi6_next, "hold");   /* ##CHRIS (M3)', 'g3_virial_block_sample("hold"); g3_snapshot_sample(0);',
                'g3_psi6_series_sample(psi6_series, &psi6_next, "record");', 'g3_virial_block_sample("record"); g3_snapshot_sample(1);'):
        i = find(F, pat); quote("f72d764 00ALLINONE.c", F, i, 1)
    print("\nThe recorded gen3 speed-of-sound commands (every committed command.txt on engine-gen3; environment variables are not in "
          "command.txt), by the registered flags they carry:\n")
    KEYS = ["--gen3-seeding=lattice", "--gen3-exact-box", "--gen3-psi6-compartments", "--psi6-every", "--gen3-virial-blocks",
            "--gen3-snapshots", "--gen3-lattice-orientation", "--target-oscillations", "--record-sigma-time"]
    files = [f for f in git("ls-files", "hspist3").split("\n") if f.endswith("/command.txt") and "experiments_gen3_s10_261012" not in f]
    combos = collections.Counter(); ex = {}
    for f in files:
        s = open(os.path.join(ROOT, f), errors="replace").read()
        if "--engine=gen3" not in s or "speed_of_sound" not in s: continue
        have = tuple(k for k in KEYS if k in s); combos[have] += 1; ex.setdefault(have, f)
    print("| runs | flags | example |\n|---|---|---|")
    for have, n in combos.most_common():
        print(f"| {n} | {' '.join(k[2:] for k in have) or '(none of them)'} | {ex[have].replace('hspist3/', '')} |")
    both = [h for h in combos if "--target-oscillations" in h]
    with_any = set(k for h in both for k in h)
    print(f"\nrecorded with --target-oscillations: {', '.join(k[2:] for k in KEYS if k in with_any and k != '--target-oscillations')}; "
          f"never with --target-oscillations: {', '.join(k[2:] for k in KEYS if k not in with_any and k not in ('--target-oscillations', '--record-sigma-time'))}")


# ---------------------------------------------------------------------------------------------------------------- grid
ETA_B = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, math.pi / 8, 0.45, 0.50, 0.55, 0.60, 0.63, 0.66, 0.68, 0.69]
ETA_D = [0.02, 0.035, 0.075]; ETA_S = [0.60, 0.63, 0.66, 0.68, 0.69]
ALPHA = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20]


def grid():
    out = []
    for blk, Ns_, etas, alphas in (("B", (100, 400, 900, 1600), ETA_B, ALPHA), ("D", (100, 400), ETA_D, ALPHA), ("S", (100, 400, 900, 1600), ETA_S, (1, 10))):
        for N in Ns_:
            for eta in etas:
                H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta)
                for a in alphas:
                    M = int(round(a * N)); nup = RUNS.nu_pred_binary(N, H, L0, M)
                    nukr = T.kr_cs(eta) * T.k_root(M / N) / (2 * math.pi * (L0 - 1.025))
                    raw = math.ceil((200.0 / nup) / DTF)
                    out.append(dict(blk=blk, N=N, eta=eta, a=a, M=M, H=H, L0=L0, nup=nup, nukr=nukr, raw=raw,
                                    steps=min(max(raw, 10000), 400000000), Teq=1000.0 if eta <= 0.60 + 1e-12 else 5000.0))
    return out


def item_d(G):
    print("\n## (d) nu_pred in the trace header under the exact box\n")
    F = src("f72d764", "hspist3/00ALLINONE.c")
    i = find(F, "const double eta_nominal_const ="); quote("f72d764 00ALLINONE.c", F, i, 3)
    i = find(F, "const double predicted_frequency =", i); quote("f72d764 00ALLINONE.c", F, i, 6)
    i = find(F, "static double speed_sound_henderson_frequency_estimate("); quote("f72d764 00ALLINONE.c", F, i - 3, 3)
    i = find(F, "const double L_eff = L0_sigma - diameter_sigma;"); quote("f72d764 00ALLINONE.c", F, i, 1)
    i = find(F, "L0_UNITS = L0;"); quote("f72d764 00ALLINONE.c", F, i, 1)
    i = find(F, '"Target_Oscillations,Predicted_Frequency,Planned_Steps,Planned_Duration");'); quote("f72d764 00ALLINONE.c", F, i - 1, 2)
    R = [x["nup"] / x["nukr"] for x in G]
    lo = min(G, key=lambda x: x["nup"] / x["nukr"]); hi = max(G, key=lambda x: x["nup"] / x["nukr"])
    print(f"\nover the {len(G)} (cell, mass) pairs of blocks B, D and S: nu_pred finite and > 0 in all: "
          f"{all(math.isfinite(x['nup']) and x['nup'] > 0 for x in G)}; nu_pred / nu_KR (KR c_s, L_eff = L0 - 1.025) from "
          f"{min(R):.4f} (N {lo['N']}, eta {lo['eta']:.4g}, alpha {lo['a']}) to {max(R):.4f} (N {hi['N']}, eta {hi['eta']:.4g}, alpha {hi['a']})")
    print("reduce_B.py (engine-gen3) takes nu_pred from the trace (tests_20260913._load: the first row's Predicted_Frequency) and uses it for "
          "the analysed prefix (TD = 200 predicted periods), the argmax floor (bin TD / X_EDGE = 80, i.e. 0.4 nu_pred) and the ACF "
          "length (20 predicted periods); the record is ceil(200 / nu_pred / dt) steps, so the prefix is the whole record:")
    B = src("engine-gen3", "hspist3/cluster/confinement_20261013/reduce_B.py")
    for pat in ("t, x, nup = T._load(p)", "n = T._prefix(t, nup, TD)", "nu = (k + int(np.argmax(P[k:]))) * df", "y = x[:n] - x[:n].mean(); nl ="):
        quote("engine-gen3 cluster/confinement_20261013/reduce_B.py", B, find(B, pat), 1)
    V = src("engine-gen3", "hspist3/validation/tests_20260913.py")
    quote("engine-gen3 validation/tests_20260913.py", V, find(V, 'return d["Time"].to_numpy(float), d["Displacement(σ)"].to_numpy(float), float(r0["Predicted_Frequency"])'), 1)
    i = find(V, "def _prefix(t, nup, N):"); quote("engine-gen3 validation/tests_20260913.py", V, i, 5)
    print("the trace's Time starts at the release (B1: first row "
          f"{open(glob.glob(os.path.join(DATA, 'B1', 'wall_x_positions_L0_*_run0.csv'))[0]).readlines()[1].split(',')[0]} sigma-time after a "
          "1000 sigma-time hold), so the prefix does not depend on the hold.")


def item_h(G):
    print("\n## (h) The step limits of the B command against the longest record\n")
    F = src("f72d764", "hspist3/00ALLINONE.c")
    i = find(F, "static int speed_sound_steps_for_target_cycles("); quote("f72d764 00ALLINONE.c", F, i, 16)
    W = src("engine-gen3", "hspist3/cluster/gen3_koa_261009/testG_worker_koa.sh")
    quote("engine-gen3 cluster/gen3_koa_261009/testG_worker_koa.sh", W, find(W, "--oscillation-min-steps=10000 --oscillation-max-steps=400000000"), 1)
    mx = max(G, key=lambda x: x["raw"]); mn = min(G, key=lambda x: x["raw"])
    print(f"\nplanned record steps over the {len(G)} (cell, mass) pairs (safety 1.0, target 200, dt = (double) 0.4f / 24):")
    print(f"- longest: {mx['raw']:,} steps = {mx['raw'] * DTF:,.0f} sigma-time (block {mx['blk']}, N {mx['N']}, eta {mx['eta']:.4g}, alpha {mx['a']}); "
          f"the maximum is 400,000,000 (factor {4e8 / mx['raw']:.2f}); INT_MAX 2,147,483,647")
    print(f"- shortest: {mn['raw']:,} steps (block {mn['blk']}, N {mn['N']}, eta {mn['eta']:.4g}, alpha {mn['a']}); the minimum is 10,000")
    print(f"- records cut short by the maximum: {sum(x['raw'] > 4e8 for x in G)}; lengthened by the minimum: {sum(x['raw'] < 10000 for x in G)}")
    print(f"- hold steps: 60,000 (T_eq 1000) and 300,000 (T_eq 5000), parsed as int (--wall-hold-steps)")


# ---------------------------------------------------------------------------------------------------------------- (e)
def item_e():
    print("\n## (e) The nu_d fit and the averaged-ACF fit (notes sec. 4.4.13 item 4a, as Test G and methods sec. 16 used them)\n")
    Fp = "hspist3/validation/resched_testT_followup_261007.py"; Fl = src("engine-gen3", Fp)
    print("The model, the start values, the bounds and the window (the whole stored ACF):")
    quote(f"engine-gen3 {Fp[8:]}", Fl, find(Fl, "def _g(t, A, tT, B, tr, om, ph):"), 2)
    quote(f"engine-gen3 {Fp[8:]}", Fl, find(Fl, "def fit_acf(a, dt, nu_ref):"), 5)
    print("Per trajectory (nu_d), the reference frequency and the failure criterion (any exception of curve_fit: no convergence within "
          "maxfev = 200000 evaluations; a fit that ends on a bound is not a failure):")
    i = find(Fl, "ref = {M: 0.5 * (nu[\"minimal\"][M].mean() + nu[\"legacy\"][M].mean()) for M in Ms}"); quote(f"engine-gen3 {Fp[8:]}", Fl, i, 1)
    i = find(Fl, "try: v.append(fit_acf(a, e[\"dt\"], ref[M])[4] / (2 * math.pi))"); quote(f"engine-gen3 {Fp[8:]}", Fl, i, 2)
    print("The seed-averaged ACF (tau_r, hence nu_0 and Q): the same model and start values on the mean of the trajectories' ACFs cut to "
          "the shortest; no failure handling there (an exception stops the script):")
    i = find(Fl, "p = fit_acf(np.mean([a[:min(len(b) for b in e[\"acf\"])] for a in e[\"acf\"]], axis=0), e[\"dt\"], ref[M])"); quote(f"engine-gen3 {Fp[8:]}", Fl, i, 3)
    Gp = "hspist3/validation/testG_verdict_261009.py"; Gl = src("engine-gen3", Gp)
    print("Test G's nu_d (the G3 target values 0.0370095 and 0.0172665), the same fit with nu_ref = the mean of the two engines' argmax means:")
    i = find(Gl, "ref = 0.5 * (D[(\"gen3\", M)][\"red\"][\"nu\"].mean() + D[(\"gen2\", M)][\"red\"][\"nu\"].mean()); v = {}"); quote(f"engine-gen3 {Gp[8:]}", Gl, i, 1)
    i = find(Gl, "try: w.append(F.fit_acf(z[f\"run{r}\"].astype(float), dt, ref)[4] / (2 * math.pi))"); quote(f"engine-gen3 {Gp[8:]}", Gl, i, 2)
    print("The ACF itself (reduce_B.py): the mean-removed position over the prefix, unbiased normalisation, lags to 20 predicted periods:")
    B = src("engine-gen3", "hspist3/cluster/confinement_20261013/reduce_B.py")
    i = find(B, "y = x[:n] - x[:n].mean(); nl ="); quote("engine-gen3 cluster/confinement_20261013/reduce_B.py", B, i, 2)


# ---------------------------------------------------------------------------------------------------------------- (f)
def zkr(eta):
    return float(SOS.Z_kolafa_rottner_2006(np.array([eta]))[0])


def bulk_rows(N, eta, L0):
    H = 10.0 * math.sqrt(N / 100.0); bp = 4 * eta / math.pi * zkr(eta)
    return bp / math.sqrt(2 * math.pi) * (4 * H + 4 * L0)       # rows per sigma-time: outer walls 2 (2 L0) + 2 H, divider faces 2 H


def ev_stats(lines):
    c = collections.Counter(); n = b = 0; t0 = t1 = None
    for k, ln in enumerate(lines):
        if k == 0: hb = len(ln); continue
        v = ln.split(",", 2); c[v[1]] += 1; n += 1; b += len(ln); t = float(v[0]); t0 = t if t0 is None else t0; t1 = t
    return c, n, b + hb, t0, t1


def item_f():
    print("\n## (f) Part A: the held path under gen3\n")
    F = src("f72d764", "hspist3/00ALLINONE.c")
    print("The start and the exact box: gen3 refuses both outside the speed-of-sound experiment (the AF command is --experiment=energy_transfer):")
    i = find(F, "else if ((cli_gen3_seeding || cli_gen3_exact_box || cli_gen3_dump_initial) && !sos)"); quote("f72d764 00ALLINONE.c", F, i, 2)
    W = src("engine-gen3", "hspist3/cluster/confinement_20261013/conf_worker.sh")
    i = find(W, 'HD_PISTON_EVENTS="$d/ev_${seed}.csv" "$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet', find(W, 'elif [ "$mode" = AF ]; then'))
    quote("engine-gen3 cluster/confinement_20261013/conf_worker.sh", W, i, 6)
    print("The divider's position: a fraction of the requested 2 L0, applied to the box width on the 1/48 grid; the validator then "
          "compares it with the request to 1e-5 sigma:")
    i = find(F, "float frac = pos_sigma / fmaxf(1e-6f, span_sigma);"); quote("f72d764 00ALLINONE.c", F, i, 1)
    i = find(F, "all_wall_positions[w] = XW1 + frac * width_px;"); quote("f72d764 00ALLINONE.c", F, i, 1)
    i = find(F, "const double position_tolerance_sigma = 1e-5;"); quote("f72d764 00ALLINONE.c", F, i, 1)
    A = src("a643bc7", "hspist3/00ALLINONE.c")
    print("The same three rules in a643bc7 (Part A has the same answer with either binary):")
    for pat in ("else if ((cli_gen3_seeding || cli_gen3_exact_box || cli_gen3_dump_initial) && !sos)",
                "all_wall_positions[w] = XW1 + frac * width_px;", "const double position_tolerance_sigma = 1e-5;"):
        quote("a643bc7 00ALLINONE.c", A, find(A, pat), 1)
    print("\n| probe | what | exit | the lines that decide |\n|---|---|---|---|")
    what = {"AF1": "eta 0.30 cell, --gen3-seeding=lattice", "AF2": "eta 0.30 cell, --gen3-exact-box",
            "AF3": "eta 0.30 cell (L0 off the 1/48 grid), legacy start, divider at L0", "AF4": "L0 = 628/48 (on the grid), divider at L0 + 0.0123",
            "AF5": "N 100, L0 = 273/48 (eta 0.6905), legacy start", "AF6": "N 1600, L0 = 1092/48 (eta 0.6905), 50 sigma-time held",
            "AF5v": "AF5 without --quiet (the seeder's lines)", "AF6v": "AF6 without --quiet (the seeder's lines)"}
    pat = re.compile(r"(STOP \(--engine=gen3\).*|ABORTING INVALID RUN.*|WARNING: L_0 not on the 1/48.*|\[SEED\].*|✅ Particles placed.*|"
                     r"\[EDMD3-HEALTH\] \S+ seed=\d+: clean=\S+)")
    for c in what:
        L = run_log(c); rc = open(os.path.join(DATA, c, "exit_code.txt")).read().strip()
        hits = list(dict.fromkeys(m.group(1) for m in pat.finditer(L)))
        print(f"| {c} | {what[c]} | {rc} | {'; '.join(hits)} |")
    tr = open(os.path.join(DATA, "AF4", "tr.csv")).readlines(); h = [x.strip() for x in tr[0].split(",")]; r1 = [x.strip() for x in tr[1].split(",")]
    print(f"\nAF4: the divider at t = 0 at {r1[h.index('W0_x_sigma')]} sigma (requested {628 / 48 + 0.0123:.6f}): an off-grid position is taken "
          "when the box is on the grid.")
    reg = [x for x in grid() if x["blk"] == "B" and x["a"] == 1]
    on = [x for x in reg if abs(48 * x["L0"] - round(48 * x["L0"])) < 1e-9]
    print(f"cells of block B whose L0 is on the 1/48 grid (they would pass the validator): {len(on)} of {len(reg)}: "
          + ", ".join(f"N {x['N']} eta {x['eta']:.4g} (L0 {x['L0']:.6g})" for x in on))
    print("\nThe event log (HD_PISTON_EVENTS; under gen3 one row per outer-wall, divider-face and piston collision):")
    E = src("f72d764", "hspist3/edmd_core/edmd_gen3.c")
    i = find(E, "one row per outer-wall, divider-face and piston collision."); quote("f72d764 edmd_core/edmd_gen3.c", E, i - 1, 2)
    meas = []
    tg = os.path.join(MAIN, "hspist3", "experiments_gen3_gate_261009", "testG", "gen3", "epi8_H_H10_L10", "x_0", "ev_974863060.csv.gz")
    with gzip.open(tg, "rt") as fh: c, n, b, t0, t1 = ev_stats(fh.readlines())
    meas.append(("Test G (Mac), gen3, N 100, eta pi/8, held, ev_974863060", 100, math.pi / 8, 10.0, c, n, b, t0, t1))
    with open(os.path.join(DATA, "AF6", "ev.csv")) as fh: c, n, b, t0, t1 = ev_stats(fh.readlines())
    meas.append(("AF6, N 1600, eta 0.6905, held 50 sigma-time", 1600, 1600 * math.pi / 4 / (2 * 1092 / 48 * 40), 1092 / 48, c, n, b, t0, t1))
    print("\n| log | rows | bytes | bytes per row | sigma-time | rows per sigma-time | of them divider (D0) | bulk estimate | measured / bulk |\n|---|---|---|---|---|---|---|---|---|")
    cal = {}
    for lab, N, eta, L0, c, n, b, t0, t1 in meas:
        rate = n / (t1 - t0); be = bulk_rows(N, eta, L0); cal[N] = (rate / be, b / n)
        print(f"| {lab} | {n:,} | {b:,} | {b / n:.1f} | {t1 - t0:.2f} | {rate:.1f} | {c['D0'] / n:.1%} | {be:.1f} | {rate / be:.3f} |")
    print("(bulk estimate: beta P / sqrt(2 pi) per unit wall length with beta P = rho Z_KR, the walls 4 H + 4 L0 long)")
    f16, bpr = cal[1600]
    print(f"\nPart A's event log at N = 1600 per (seed, position): the bulk estimate times {f16:.3f} (AF6), {bpr:.1f} bytes per row, "
          "held T_eq + 5000 sigma-time (+ the 20 sigma-time tail):\n")
    print("| eta | L0 | rows per sigma-time | held sigma-time | rows | size |\n|---|---|---|---|---|---|")
    for eta in ETA_B:
        L0 = 1600 * math.pi / (8 * 40 * eta); Teq = 1000 if eta <= 0.60 + 1e-12 else 5000; r = bulk_rows(1600, eta, L0) * f16
        nn = r * (Teq + 5000 + 20)
        print(f"| {eta:.4g} | {L0:.4f} | {r:.1f} | {Teq + 5000} | {nn:,.0f} | {nn * bpr / 1e6:,.0f} MB |")
    print("\nreduce_AF.py's window: [t0, t1) from its arguments; the worker passes t0 = 200 and t1 = the release:")
    R = src("engine-gen3", "hspist3/cluster/confinement_20261013/reduce_AF.py")
    quote("engine-gen3 cluster/confinement_20261013/reduce_AF.py", R, find(R, "ev, tr, out, t0, t1 = sys.argv[1]"), 1)
    quote("engine-gen3 cluster/confinement_20261013/reduce_AF.py", R, find(R, "e = pd.read_csv(ev, usecols="), 1)
    quote("engine-gen3 cluster/confinement_20261013/reduce_AF.py", R, find(R, 'd = e[(e["kind"] == "D0")'), 1)
    quote("engine-gen3 cluster/confinement_20261013/conf_worker.sh", W, find(W, 't1=$(awk -v h="$hold" \'BEGIN{printf "%.9f", h * 0.4 / 24.0}\')'), 2)
    K = src("engine-gen3", "hspist3/cluster/gen3_gate_261009/af_stream_reduce.py")
    quote("engine-gen3 cluster/gen3_gate_261009/af_stream_reduce.py", K, find(K, "data = fh.read()"), 1)


# ---------------------------------------------------------------------------------------------------------------- (g)
def item_g(G):
    print("\n## (g) MaxArraySize on KOA\n")
    S = src("engine-gen3", "hspist3/cluster/KOA_RUNSHEET_261002.md")
    i = find(S, "The cap is our own choice. KOA's limits (printed by Chris on 2026-10-08):"); quote("engine-gen3 cluster/KOA_RUNSHEET_261002.md", S, i, 5)
    n = collections.Counter(x["blk"] for x in G)
    print(f"\nregistered trajectories: B {n['B'] * 25:,} (25 seeds), D {n['D'] * 25:,}, S {n['S'] * 10 * 2:,} (10 seeds, two starts); "
          f"total {n['B'] * 25 + n['D'] * 25 + n['S'] * 20:,}: below 25001 even at one trajectory per array task")


def main():
    print("# Report S1-0, gate G0 (sec. 4.8.10): the code facts, printed by experiments_gen3_s10_261012/s10_tables.py\n")
    G = grid()
    item_a(); item_bc(); item_d(G); item_e(); item_f(); item_g(G); item_h(G)


if __name__ == "__main__":
    main()
