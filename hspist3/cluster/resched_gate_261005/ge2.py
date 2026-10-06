#!/usr/bin/env python3
"""##CHRIS 2026-10-05 (261012 sec. 4.4, engine gate G-E2): the minimal divider rescheduling (default) against the legacy full
reschedule (--legacy-resched), SAME binary, SAME seed, on two cells:

  (a) smoke: the determinism trajectory of cluster/confinement_pilot.py -- eta = pi/8, H = L_0 = 10, N_s = 50, M = 50
      (alpha = 0.5), seed run_seed(20261013, 0, 0, 0), 25 oscillations; speed-of-sound mode: divider held for 2000 steps,
      then released (free divider). Outputs: the wall trace and the psi6 file.
  (b) afix:  the A-fixed cell epi8_H_H10_L10 at x_0 (the centre), the campaign's first x_0 seed (9700), the campaign
      protocol (conf_worker.sh AF: hold 312000 steps = divider held over [0, 5200), then 1200 steps released), event log on;
      reduced by reduce_AF.py over the held window [200, 5200) as in the campaign.
Every run has HD_KE_TRACE=1, so the binary prints [EDMD-ENERGY] lines (E_gas + E_divider at %.17g) at 0 = state loaded,
1 = release, 2 = end.

Registered before any KOA run (261012 sec. 4.4):
  byte identity  cmp of every output pair (a: trace, psi6; b: ev, tr, red). Not expected: the minimal path keeps pair events
                 that the legacy path recomputes from drifted positions, which changes event times in the last bit, and the
                 dynamics amplifies that. If not identical, the first differing trace row / event-log line is printed.
  energy         per run and phase (hold 0 -> 1, record 1 -> 2): |dE/E| of E_tot. PASS if the minimal run's value is
                 <= max(10 x the legacy run's value in the same phase, 1e-12) in every phase ("within legacy tolerance").
  divider ledger (b) u_wall_max = 0 and W_div = 0 in the held window, both runs. Forces F_L, F_R and the event counts are
                 printed new/legacy for information only (independent realisations once the runs diverge).
  health         no [EDMD-HEALTH] line in any run (the new build prints it in both modes, with past_events).
  policy         each run printed the [EDMD-RESCHED] line of its policy.
  279282b        with --bin-old (the recorded 279282b binary of ~/harddisks): the same two runs with the OLD binary, and the
                 legacy path's outputs compared with them byte for byte. Registered expectation [INFERENCE]: IDENTICAL -- in
                 the legacy path the new code only adds integer counters, an Event.cb that is never read for divider
                 events there, and print-only lines; the floating-point operations and their order are unchanged and the
                 build has -ffp-contract=off. Reported, not part of the plan author's verdict rule.
usage (from hspist3/):  python3 cluster/resched_gate_261005/ge2.py run --bin ./00ALLINONE --out <dir> [--bin-old <279282b binary>]
                        python3 cluster/resched_gate_261005/ge2.py compare --out <dir>
"""
import argparse, glob, math, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, CL); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)

AF_TASKS = os.path.join(CL, "confinement_20261013", "tasks_AF_epi8_H_H10_L10.txt")
POLICIES = (("minimal", []), ("legacy", ["--legacy-resched"]))
EN = re.compile(r"\[EDMD-ENERGY\] (\d) .*?E_tot=(\S+) resched=(\w+)")


def af_task():
    for line in open(AF_TASKS):
        f = line.split()
        if f[0] == "AF" and f[1].endswith("/x_0"):
            return dict(xw=f[2], seed=f[3], L0=f[4], H=f[5], NS=int(f[6]), hold=int(f[7]), post=int(f[8]), every=f[9])
    sys.exit("no x_0 line in " + AF_TASKS)


def run(a):
    import confinement_pilot as CP
    a.bin, a.out = os.path.abspath(a.bin), os.path.abspath(a.out)
    if os.path.exists(a.out): sys.exit(f"{a.out} exists -- not overwriting")
    t = af_task(); procs = []
    env = dict(os.environ, HD_KE_TRACE="1")
    runs = [(pol, a.bin, extra) for pol, extra in POLICIES]
    if a.bin_old: runs.append(("279282b", os.path.abspath(a.bin_old), []))
    for pol, binp, extra in runs:
        d = os.path.join(a.out, f"smoke_{pol}"); os.makedirs(d)
        c = CP.cmd(binp, 50, 0, d, target=25) + extra
        open(os.path.join(d, "command.txt"), "w").write("HD_KE_TRACE=1 " + " ".join(c) + "\n")
        procs.append((d, subprocess.Popen(c, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT, env=env, cwd=d)))
        d = os.path.join(a.out, f"afix_{pol}"); os.makedirs(d); s = t["seed"]
        c = [binp, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0",
             "--seed-drift-order=drift-first", f"--energy-transfer-summary={d}/summary_{s}.csv",
             f"--energy-transfer-trace={d}/tr_{s}.csv", f"--trace-every={t['every']}", f"--particles={2 * t['NS']}",
             f"--particles-boxes={t['NS']},{t['NS']}", "--particle-radius=0.5", f"--l0={t['L0']}", f"--height={t['H']}",
             "--num-walls=1", f"--wall-positions={t['xw']}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05",
             "--wall-thickness-vis=0.05", "--eff-output=wall-ke", f"--wall-hold-steps={t['hold']}", f"--steps={t['post']}",
             "--fixed-dt=0.4", "--kbt1", f"--seed={s}"] + extra
        open(os.path.join(d, "command.txt"), "w").write(f"HD_KE_TRACE=1 HD_PISTON_EVENTS={d}/ev_{s}.csv " + " ".join(c) + "\n")
        procs.append((d, subprocess.Popen(c, stdout=open(os.path.join(d, f"run_{s}.log"), "w"), stderr=subprocess.STDOUT,
                                          env=dict(env, HD_PISTON_EVENTS=f"{d}/ev_{s}.csv"), cwd=d)))
    rc = [(d, p.wait()) for d, p in procs]
    for d, r in rc: print(f"{os.path.basename(d)}: exit {r}")
    t1 = t["hold"] * 0.4 / 24.0
    for pol, _, _ in runs:
        d = os.path.join(a.out, f"afix_{pol}"); s = t["seed"]
        r = subprocess.run([sys.executable, os.path.join(CL, "confinement_20261013", "reduce_AF.py"), f"{d}/ev_{s}.csv",
                            f"{d}/tr_{s}.csv", f"{d}/red_{s}.csv", "200", f"{t1:.9f}"])
        print(f"reduce_AF afix_{pol}: exit {r.returncode}")
    return 0 if all(r == 0 for _, r in rc) else 1


def first_diff(p, q):
    """(identical?, 1-based line of the first difference, the two lines, line counts)"""
    A = open(p, errors="ignore").read().split("\n"); B = open(q, errors="ignore").read().split("\n")
    for i, (x, y) in enumerate(zip(A, B)):
        if x != y: return False, i + 1, x, y, len(A), len(B)
    return (len(A) == len(B)), min(len(A), len(B)) + 1, "", "", len(A), len(B)


def energies(log):
    e = {}; pol = set()
    for m in EN.finditer(open(log, errors="ignore").read()):
        e[int(m.group(1))] = float(m.group(2)); pol.add(m.group(3))
    return e, pol


def compare(a):
    a.out = os.path.abspath(a.out); s = af_task()["seed"]; ok = dict(energy=True, ledger=True, health=True, policy=True)
    print(f"## G-E2 -- minimal vs legacy rescheduling, same binary, same seed ({a.out})\n")
    print("### Byte identity\n\n| cell | file | identical | lines (min / leg) | first differing line | minimal | legacy |\n|---|---|---|---|---|---|---|")
    pairs = [("smoke", os.path.basename(sorted(glob.glob(os.path.join(a.out, "smoke_minimal", "wall_x_positions_*_run0.csv")))[0])),
             ("smoke", "speed_of_sound_psi6.csv"), ("afix", f"ev_{s}.csv"), ("afix", f"tr_{s}.csv"), ("afix", f"red_{s}.csv")]
    for cell, fn in pairs:
        p, q = os.path.join(a.out, f"{cell}_minimal", fn), os.path.join(a.out, f"{cell}_legacy", fn)
        same, i, x, y, na, nb = first_diff(p, q)
        cut = lambda z: z[:70] + (" ..." if len(z) > 70 else "")
        print(f"| {cell} | {fn} | {'IDENTICAL' if same else 'no'} | {na} / {nb} | {'-' if same else i} | "
              f"{'' if same else '`' + cut(x) + '`'} | {'' if same else '`' + cut(y) + '`'} |")
    print("\n### Energy (E_tot = E_gas + E_divider, [EDMD-ENERGY] lines)\n")
    print("| cell | policy | E_0 | E_1 | E_2 | hold abs(dE/E) 0->1 | record abs(dE/E) 1->2 |\n|---|---|---|---|---|---|---|")
    rel = {}
    for cell, log in (("smoke", "run.log"), ("afix", f"run_{s}.log")):
        for pol, _ in POLICIES:
            e, pols = energies(os.path.join(a.out, f"{cell}_{pol}", log))
            if pols != {pol}: ok["policy"] = False
            h = abs(e[1] / e[0] - 1); r = abs(e[2] / e[1] - 1); rel[(cell, pol)] = (h, r)
            print(f"| {cell} | {pol} | {e[0]:.17g} | {e[1]:.17g} | {e[2]:.17g} | {h:.3e} | {r:.3e} |")
    print("\n| cell | phase | minimal | legacy | limit max(10 x legacy, 1e-12) | PASS |\n|---|---|---|---|---|---|")
    for cell in ("smoke", "afix"):
        for k, ph in ((0, "hold"), (1, "record")):
            m, l = rel[(cell, "minimal")][k], rel[(cell, "legacy")][k]; lim = max(10 * l, 1e-12); good = m <= lim
            ok["energy"] &= good
            print(f"| {cell} | {ph} | {m:.3e} | {l:.3e} | {lim:.3e} | {'yes' if good else '**NO**'} |")
    import pandas as pd
    print("\n### Divider ledger and forces, A-fixed cell (held window [200, 5200))\n")
    print("| policy | F_L | F_R | n_L | n_R | T_L | T_R | u_wall_max | W_div | t_last |\n|---|---|---|---|---|---|---|---|---|---|")
    for pol, _ in POLICIES:
        r = pd.read_csv(os.path.join(a.out, f"afix_{pol}", f"red_{s}.csv")).iloc[0]
        ok["ledger"] &= (r["u_wall_max"] == 0.0 and r["W_div"] == 0.0)
        print(f"| {pol} | {r['F_L']:.6f} | {r['F_R']:.6f} | {int(r['n_L'])} | {int(r['n_R'])} | {r['T_L']:.9f} | {r['T_R']:.9f} | "
              f"{r['u_wall_max']:g} | {r['W_div']:g} | {r['t_last']:.3f} |")
    print("\n### Health and policy lines\n\n| run | [EDMD-HEALTH] lines | [EDMD-RESCHED] |\n|---|---|---|")
    for cell, log in (("smoke", "run.log"), ("afix", f"run_{s}.log")):
        for pol, _ in POLICIES:
            t = open(os.path.join(a.out, f"{cell}_{pol}", log), errors="ignore").read()
            nh = t.count("[EDMD-HEALTH]"); ok["health"] &= nh == 0
            rs = re.findall(r"\[EDMD-RESCHED\] divider events: (\w+)", t)
            ok["policy"] &= rs == [pol]
            print(f"| {cell}_{pol} | {nh} | {','.join(rs) or 'none'} |")
    if os.path.isdir(os.path.join(a.out, "smoke_279282b")):
        print("\n### Legacy path of the new binary vs the 279282b binary (same seed; expectation IDENTICAL, reported)\n")
        print("| cell | file | identical | first differing line |\n|---|---|---|---|")
        for cell, fn in pairs:
            same, i, x, y, na, nb = first_diff(os.path.join(a.out, f"{cell}_legacy", fn), os.path.join(a.out, f"{cell}_279282b", fn))
            print(f"| {cell} | {fn} | {'IDENTICAL' if same else '**no**'} | {'-' if same else i} |")
    print("\nG-E2: " + "; ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in ok.items())
          + " (byte identity is reported, not gated: see the table)")
    return 0 if all(ok.values()) else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "compare"])
    ap.add_argument("--bin"); ap.add_argument("--bin-old"); ap.add_argument("--out", required=True); a = ap.parse_args()
    sys.exit({"run": run, "compare": compare}[a.what](a))
