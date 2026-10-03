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
    print("\n### Reproduction gate: the pre-registered rule with the planning eps0 against the committed task files\n")
    print("| cell | seeds/position (rule) | seeds/position (tasks file) | reproduced |\n|---|---|---|---|")
    ok_all = True
    for c in [c for c in cells if c["eta_lab"] == "0.39"]:
        n_tasks = sum(1 for _ in open(os.path.join(CONF, f"tasks_A_{cid(c)}.txt"))) // 5
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
    print(f"\npooled eps at N_s = 50: {eps_pool:.4f} -> eps0 (pi/8 pilot) = {eps0_pilot:.4f}  (planning value {eps0_plan:.4f}, "
          f"ratio {eps0_pilot/eps0_plan:.3f})")
    print("\n### Seeds per position for conf_A_0.39, by the pre-registered rule with the pilot's eps0\n")
    print("| cell | T per position (plan) | seeds/position (plan) | T per position (pilot eps0) | seeds/position (new) |\n|---|---|---|---|---|")
    for c in [c for c in cells if c["eta_lab"] == "0.39"]:
        Tn = c["Tpos"] * (eps0_pilot / eps0_plan) ** 2
        print(f"| {cid(c)} | {c['Tpos']:.3g} | {c['nseed']} | {Tn:.3g} | {math.ceil(Tn / PR.T_SEED_A)} |")
    ok = len(rows) == 20 and health == 0 and conv_ok
    print(f"\n**GATE 4: {'PASS' if ok else 'FAIL'}** -- conf_A_0.39 uses the new seeds per position "
          f"(its tasks files must be regenerated before Round 2). conf_A_0.10 is NOT affected: its eps0 was measured at "
          f"eta = 0.10005 (Level 3 c0), and sec. 1.8 allows the pilot to change only (A) at pi/8.")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
