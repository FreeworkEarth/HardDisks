#!/usr/bin/env python3
"""##CHRIS 2026-10-02 (Task R1): gate 4 of 261012 sec. 1.7 applied to the method-A pilot at pi/8, on the Mac.

Gate 4, verbatim (261012 sec. 1.7): "Pilot at pi/8: 4 seeds per position at the anchor. This measures the noise
coefficient, fixes the record, and checks that the steps-to-sigma-time conversion gives 5000 sigma-time per seed (read
from the event-log time range, not computed from --steps)."  Sec. 1.8: "Only the pilot of gate 4: it may change the
record length and the cost of (A) at pi/8, by the rule already written."
Noise model (sec. 1.4): sigma_F/F = eps0 (100/N_s)^(1/2) / sqrt(T); the rule (paper1_confinement_prereg_20261012.py:163-165):
    Tpos = ((sqrt(130)/12) eps FoverK / (dL NOISE_MAX))^2 / 2,  eps = eps0 sqrt(100/N_s),  seeds/position = ceil(Tpos/5000)
so Tpos is proportional to eps0^2 [DERIVATION]: Tpos_new = Tpos_plan (eps0_pilot/eps0_plan)^2.

Input: the pilot's per-seed reductions red_<seed>.csv (written ON KOA by conf_worker.sh mode A via reduce_A.py) and
run_<seed>.log, copied to the Mac (runsheet step 8). Per seed: F_L, F_R from the divider's event log over the window
[200, t_last]; window = t_last - 200 is read from the EVENT LOG (reduce_A.py), which is the conversion check.
eps measured per (position, face) = SD(F over the 4 seeds)/mean(F) * sqrt(window), pooled in quadrature over the ten
(position, face) groups; each face has N_s = 50 behind it, so eps0_pilot = eps_pool / sqrt(100/50).
Checks: 20 seeds present, health lines 0, |window - 5000|/5000 <= 1 % for every seed (the 1 % tolerance is this script's
operational reading of "gives 5000 sigma-time"; stated, not pre-registered), and -- before anything else -- that the
planning eps0 reproduces the committed task files' seeds per position exactly.
usage: python3 hspist3/cluster/gate4_pilot_261002.py [pilot_dir]      (--dry: only the reproduction gate, no pilot data)

##CHRIS 2026-10-02 (Task U1): also printed now -- the error of the pooled eps0 (chi-square: 10 groups x 3 degrees of
freedom, relative SE of an SD = 1/sqrt(2 dof)), a Bartlett test that the ten groups share one relative variance
(information only), the seeds per position at eps0 + 1 sigma (information only; the verdict uses eps0 itself, by the
pre-registered rule), and the conf_A_0.39 core-hour totals (Mac cost model, and at the measured KOA speed of
round_plan_261002.koa_speed()). On PASS the result is recorded in cluster/confinement_20261013/gate4_pi8_result.txt,
which gen_confinement_sbatch.py reads to write the conf_A_0.39 task files.
"""
import contextlib, glob, io, math, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import paper1_confinement_prereg_20261012 as PR
REL = "experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10"
CONF = os.path.join(HERE, "confinement_20261013")
HEALTH = re.compile(r"EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue")

def cid(c): return f"e{'0p10' if c['eta_lab'] == '0.10' else 'pi8'}_{c['scan']}_H{c['H']:g}_L{c['L0']:g}"   # as gen_confinement_sbatch

def planned():
    """Table A of the pre-registration, recomputed by its own code; deduplicated like the generator."""
    with contextlib.redirect_stdout(io.StringIO()):
        C, A = PR.main(); eps0_plan = PR.noise_eps0()[0]
    seen, out = set(), []
    for c in C:
        key = (c["eta_lab"], round(c["H"], 6), round(c["L0"], 6))
        if key not in seen: seen.add(key); out.append(c)
    return out, eps0_plan

def main():
    cells, eps0_plan = planned()
    print(f"planning eps0 (Level 3 c0, eta = 0.10005): {eps0_plan:.4f}")
    # ##CHRIS 2026-10-02 (Task U3): the task files now carry the gate-4 seeds, so the reproduction gate reads them as
    # committed BEFORE gate 4, at the planning commit PLAN_REF (git show), not from the working tree.
    import subprocess
    PLAN_REF = "70b2069"
    print(f"\n### Reproduction gate: the pre-registered rule with the planning eps0 against the task files at git {PLAN_REF}\n")
    print("| cell | seeds/position (rule) | seeds/position (tasks file) | reproduced |\n|---|---|---|---|")
    ok_all = True
    for c in [c for c in cells if c["eta_lab"] == "0.39"]:
        n_tasks = len(subprocess.run(["git", "show", f"{PLAN_REF}:hspist3/cluster/confinement_20261013/tasks_A_{cid(c)}.txt"],
                                     capture_output=True, text=True, cwd=HS, check=True).stdout.splitlines()) // 5
        ok = n_tasks == c["nseed"]; ok_all &= ok
        print(f"| {cid(c)} | {c['nseed']} | {n_tasks} | {'yes' if ok else '**NO**'} |")
    print(f"\nreproduction gate: {'PASS' if ok_all else 'FAIL -- STOP'}")
    if not ok_all: return 1
    if "--dry" in sys.argv: return 0
    d = next((a for a in sys.argv[1:] if not a.startswith("--")), os.path.join(HS, REL))
    print(f"\n### Pilot: {d}\n")
    rows, health, missing = [], 0, []
    import pandas as pd
    for pos in ("m2", "m1", "0", "p1", "p2"):
        for seed in range(9700, 9704):
            f = os.path.join(d, f"x_{pos}", f"red_{seed}.csv"); lg = os.path.join(d, f"x_{pos}", f"run_{seed}.log")
            if not os.path.exists(f): missing.append(f"x_{pos}/red_{seed}.csv"); continue
            r = pd.read_csv(f).iloc[0]
            if os.path.exists(lg): health += len(HEALTH.findall(open(lg, errors="ignore").read()))
            else: missing.append(f"x_{pos}/run_{seed}.log")
            rows.append(dict(pos=pos, seed=seed, w=float(r["window"]), FL=float(r["F_L"]), FR=float(r["F_R"])))
    print(f"seeds found: {len(rows)} of 20; missing files: {missing or 'none'}; health lines: {health}")
    # ##CHRIS 2026-10-02 (Task U2): when was each file written (mtimes preserved by rsync -a; UTC), and is any run log different
    import datetime, hashlib
    fs = sorted(glob.glob(os.path.join(d, "x_*", "*")))
    utc = lambda f: datetime.datetime.fromtimestamp(os.path.getmtime(f), datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    for kind in ("red_", "run_", "summary_"):
        m = sorted(utc(f) for f in fs if os.path.basename(f).startswith(kind))
        print(f"  {kind}*: {len(m)} files, mtime (UTC) {m[0]} .. {m[-1]}" if m else f"  {kind}*: none")
    logs = {hashlib.md5(open(f, "rb").read()).hexdigest() for f in fs if os.path.basename(f).startswith("run_")}
    print(f"  run_*.log: {len(logs)} distinct content(s); all files: {len(fs)}, mtime (UTC) {min(map(utc, fs))} .. {max(map(utc, fs))}")
    wins = np.array([q["w"] for q in rows]); conv_ok = bool(len(wins)) and bool(np.all(np.abs(wins - 5000) / 5000 <= 0.01))
    print(f"window per seed from the event log: {wins.min():.1f} .. {wins.max():.1f} sigma-time "
          f"-> conversion check (5000 +- 1 %): {'PASS' if conv_ok else 'FAIL'}")
    print("\n| position | face | F mean | SD over 4 seeds | eps = SD/mean x sqrt(window) |\n|---|---|---|---|---|")
    eps2 = []
    for pos in ("m2", "m1", "0", "p1", "p2"):
        sel = [q for q in rows if q["pos"] == pos]
        for face in ("FL", "FR"):
            v = np.array([q[face] for q in sel]); w = np.mean([q["w"] for q in sel])
            if len(v) < 2: continue
            e = v.std(ddof=1) / v.mean() * math.sqrt(w); eps2.append(e * e)
            print(f"| x_{pos} | {face[1]} | {v.mean():.5f} | {v.std(ddof=1):.5f} | {e:.4f} |")
    eps_pool = math.sqrt(np.mean(eps2)); eps0_pilot = eps_pool / math.sqrt(100 / 50)
    dof = sum(len([q for q in rows if q["pos"] == p]) - 1 for p in ("m2", "m1", "0", "p1", "p2")) * 2
    err = eps0_pilot / math.sqrt(2 * dof)
    print(f"\npooled eps at N_s = 50: {eps_pool:.4f} -> eps0 (pi/8 pilot) = {eps0_pilot:.4f} +- {err:.4f} ({dof} degrees of freedom)"
          f"  (planning value {eps0_plan:.4f}, ratio {eps0_pilot/eps0_plan:.3f})")
    from scipy.stats import bartlett                               # ##CHRIS 2026-10-02 (Task U1): information only
    grp = [np.array([q[f] for q in rows if q["pos"] == p]) for p in ("m2", "m1", "0", "p1", "p2") for f in ("FL", "FR")]
    print(f"Bartlett test, one relative variance across the 10 (position, face) groups: p = {bartlett(*[g / g.mean() for g in grp]).pvalue:.3f}")
    # ##CHRIS 2026-10-02 (Task U1), information only [INFERENCE]: divider impacts per face per sigma-time, H n Z / sqrt(2 pi)
    # (contact theorem, kT = m = 1, n = 4 eta / pi), at the planning anchor (Level 3 c0) and at the pilot; H = 10 in both
    rate = lambda e: 10 * (4 * e / math.pi) * PR.Z(e) / math.sqrt(2 * math.pi)
    r10, r39 = rate(0.10005), rate(math.pi / 8)
    print(f"impacts per face per sigma-time: {r10:.3f} (eta 0.10005) vs {r39:.3f} (pi/8), ratio {r39 / r10:.2f}; pure shot noise "
          f"would scale eps by sqrt(1/ratio) = {math.sqrt(r10 / r39):.3f}; measured eps(pi/8, N_s 50)/eps0(plan) = {eps_pool / eps0_plan:.3f}")
    import round_plan_261002 as RP
    speed = RP.koa_speed()
    print("\n### Seeds per position for conf_A_0.39, by the pre-registered rule with the pilot's eps0\n")
    print("| cell | T per position (plan) | seeds/position (plan) | T per position (pilot eps0) | seeds/position (new) | "
          "seeds/position at eps0 + 1 sigma (= amendment C3, USED) | core-h plan (Mac model) | core-h new (Mac model) |\n|---|---|---|---|---|---|---|---|")
    ch_plan = ch_new = ch_c3 = 0.0
    for c in [c for c in cells if c["eta_lab"] == "0.39"]:
        Tn = c["Tpos"] * (eps0_pilot / eps0_plan) ** 2; nn = math.ceil(Tn / PR.T_SEED_A)
        nhi = math.ceil(c["Tpos"] * ((eps0_pilot + err) / eps0_plan) ** 2 / PR.T_SEED_A)
        cn = c["cpuA"] * nn / c["nseed"]; ch_plan += c["cpuA"]; ch_new += cn; ch_c3 += c["cpuA"] * nhi / c["nseed"]
        print(f"| {cid(c)} | {c['Tpos']:.3g} | {c['nseed']} | {Tn:.3g} | {nn} | {nhi} | {c['cpuA']:.2f} | {cn:.2f} |")
    print(f"\nconf_A_0.39 core-h: plan {ch_plan:.1f} -> new {ch_new:.1f} (Mac cost model); at the measured KOA speed x{speed:.3f}: "
          f"plan {ch_plan * speed:.1f} -> new {ch_new * speed:.1f}")
    # ##CHRIS 2026-10-02 (Task V2): amendment C3 uses the eps0 + 1 sigma column
    print(f"amendment C3 (eps0 + 1 sigma = {eps0_pilot + err:.4f}): conf_A_0.39 core-h {ch_c3:.1f} (Mac cost model), "
          f"{ch_c3 * speed:.1f} at KOA speed (+{(ch_c3 - ch_new) * speed:.1f} over the gate-4 seeds)")
    ok = len(rows) == 20 and health == 0 and conv_ok
    if ok:
        res = os.path.join(CONF, "gate4_pi8_result.txt")
        open(res, "w").write(f"# gate 4 (261012 sec. 1.7 item 4), written by cluster/gate4_pilot_261002.py; pilot job 14966594 (KOA, 70b2069)\n"
                             f"eps0_pilot {float(eps0_pilot)!r}\neps0_pilot_err {float(err)!r}\neps0_plan {float(eps0_plan)!r}\nverdict PASS\n"
                             # ##CHRIS 2026-10-02 (Task V2): amendment C3 (261012 sec. 1.9) -- the seeds come from eps0 + 1 sigma
                             f"eps0_c3_upper {float(eps0_pilot + err)!r}\n")
        print(f"recorded: {os.path.relpath(res, HS)}")
    print(f"\n**GATE 4: {'PASS' if ok else 'FAIL'}** -- conf_A_0.39 uses the new seeds per position "
          f"(its tasks files must be regenerated before Round 2). conf_A_0.10 is NOT affected: its eps0 was measured at "
          f"eta = 0.10005 (Level 3 c0), and sec. 1.8 allows the pilot to change only (A) at pi/8.")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
