#!/bin/bash
# ##CHRIS 2026-10-03 (Task W2): READ-ONLY check of every confinement cell on scratch against its task file.
# Per cell: trajectories expected (task file) vs present (non-empty output), the missing list, the .build_git record(s),
# zero-size / duplicate / unexpected outputs, leftovers of killed runs (.run<r>, .stale_*, partial seeds), held
# .guard.lock directories, a timing line for cells with leftovers or missing work, and a one-line verdict.
# It writes nothing, anywhere (os.scandir / os.stat / open-for-read only).
#   expected output: B  <cell>/m_<M>/wall_x_positions_L0_*_wallmassfactor_<M>_run<r>.csv   (conf_worker.sh mode B)
#                    A  <cell>/x_<pos>/red_<seed>.csv, non-empty, 2 lines                    (conf_worker.sh mode A)
# usage (KOA, INSIDE a sandbox session, from ~/harddisks/hspist3):  bash cluster/check_cells.sh [data_root]
set -uo pipefail
cd "$(dirname "$0")/.." || exit 2
case "$(hostname)" in login*) echo "STOP: run inside a sandbox session (srun -p sandbox ... --pty /bin/bash), not on the login node"; exit 2;; esac
ROOT="${1:-/mnt/lustre/koa/scratch/charing/harddisks/hspist3}"
PY="$HOME/envs/hd/bin/python3"; [ -x "$PY" ] || PY=python3
echo "check_cells: data root $ROOT; repo $(git log --oneline -1 2>/dev/null | cut -c1-60); $(date -u +%Y-%m-%dT%H:%M:%SZ)"
exec "$PY" - "$ROOT" <<'PYEOF'
import os, re, sys, time
from collections import defaultdict
ROOT = sys.argv[1]; CONF = "cluster/confinement_20261013"
GROUPS = ["B_0.10", "B_0.39", "A_0.10", "A_0.39", "A_pilot", "Afix_0.10", "Afix_0.39", "Afix_pilot"]   # ##CHRIS 2026-10-04: + A-fixed (sec. 3)

def ls(d):
    try: return {e.name: e for e in os.scandir(d)}
    except FileNotFoundError: return None

def ranges(xs):
    xs = sorted(xs); out = []; i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[j] + 1: j += 1
        out.append(str(xs[i]) if i == j else f"{xs[i]}-{xs[j]}"); i = j + 1
    return ",".join(out)

def build_records(dirs):
    rec, norec = defaultdict(int), 0
    for d, ents in dirs.items():
        if ents is None: continue
        if ".build_git" in ents:
            rec[open(os.path.join(d, ".build_git")).read().strip()] += 1
        else: norec += 1
    return dict(rec), norec

totals = defaultdict(int)
for g in GROUPS:
    tsv = os.path.join(CONF, f"cells_{g}.tsv")
    if not os.path.exists(tsv): continue
    print(f"\n== {g}")
    for task, cell in enumerate(open(tsv).read().split(), 1):
        pre = "AF" if g.startswith("Afix") else ("A" if g[0] == "A" else "B")      # ##CHRIS 2026-10-04: AF lines = A layout
        lines = [l.split() for l in open(os.path.join(CONF, f"tasks_{pre}_{cell}.txt"))]
        exp = defaultdict(list)                                  # dir -> expected ids
        for f in lines:
            exp[os.path.join(ROOT, f[1])].append((int(f[2]), int(f[3])) if f[0] == "B" else int(f[3]))
        dirs = {d: ls(d) for d in exp}
        n_exp, missing, zero, dup, extra, malformed, partial, leftovers, locks = sum(len(v) for v in exp.values()), [], 0, 0, 0, 0, 0, [], 0
        mt, rec_t, run_secs = [], [], []
        for d, ids in exp.items():
            ents = dirs[d]
            if ents is None: missing += [(d, i) for i in ids]; continue
            if ".guard.lock" in ents: locks += 1
            if ".build_git" in ents: rec_t.append(ents[".build_git"].stat().st_mtime)
            leftovers += [n for n in ents if (n.startswith(".run") and n != ".runlog.lock") or n.startswith(".stale") or n.startswith(".failed")]
            if lines[0][0] == "B":
                found = defaultdict(list)
                for n, e in ents.items():
                    m = re.match(r"wall_x_positions_L0_\d+_wallmassfactor_(\d+)_run(\d+)\.csv$", n)
                    if m: found[(int(m.group(1)), int(m.group(2)))].append(e)
                for i in ids:
                    fs = found.get(i, [])
                    if len(fs) > 1: dup += 1
                    if not fs: missing.append((d, i))
                    elif all(e.stat().st_size == 0 for e in fs): zero += 1; missing.append((d, i))
                    else: mt.append(max(e.stat().st_mtime for e in fs))
                extra += len(set(found) - set(ids))
                if "run.log" in ents:
                    run_secs += [int(x) for x in re.findall(r"##RUN .*?\((\d+) s\)", open(os.path.join(d, "run.log"), errors="ignore").read())]
            else:
                reds = {int(m.group(1)): e for n, e in ents.items() for m in [re.match(r"red_(\d+)\.csv$", n)] if m}
                for s in ids:
                    e = reds.get(s)
                    if e is None or e.stat().st_size == 0:
                        if e is not None: zero += 1
                        missing.append((d, s))
                        if any(n.endswith(f"_{s}.csv") or n.endswith(f"_{s}.log") for n in ents): partial += 1
                    else:
                        if sum(1 for _ in open(os.path.join(d, f"red_{s}.csv"))) != 2: malformed += 1
                        mt.append(e.stat().st_mtime)
                extra += len(set(reds) - set(ids))
        rec, norec = build_records(dirs)
        n_ok = n_exp - len(missing)
        if n_ok == 0 and not any(dirs.values()): verdict = "NOT STARTED"
        elif zero or dup or malformed or locks or len(rec) > 1 or (norec and n_ok): verdict = "PROBLEM"
        elif missing: verdict = f"INCOMPLETE ({len(missing)} missing)"
        else: verdict = "COMPLETE"
        totals[verdict.split(" (")[0]] += 1
        miss_txt = ""
        if missing and n_ok:
            byd = defaultdict(list)
            for d, i in missing: byd[os.path.basename(d)].append(i[1] if isinstance(i, tuple) else i)
            miss_txt = "; missing " + " ".join(f"{k}:{ranges(v)}" for k, v in sorted(byd.items()))[:300]
        print(f"[{task}] {cell}: {n_ok}/{n_exp} present; build {rec or '-'}{'; dirs without record: ' + str(norec) if norec and n_ok else ''}"
              f"{'; zero-size ' + str(zero) if zero else ''}{'; duplicates ' + str(dup) if dup else ''}{'; unexpected ' + str(extra) if extra else ''}"
              f"{'; malformed red ' + str(malformed) if malformed else ''}{'; partial seeds ' + str(partial) if partial else ''}"
              f"{'; leftovers ' + str(len(leftovers)) if leftovers else ''}{'; HELD .guard.lock ' + str(locks) if locks else ''}{miss_txt}")
        if (missing or leftovers) and mt:
            span_h = (max(mt) - min(rec_t or mt)) / 3600
            extra_t = f", B run.log: {len(run_secs)} runs, mean {sum(run_secs) / len(run_secs):.0f} s, max {max(run_secs)} s" if run_secs else ""
            print(f"      timing: {len(mt)} outputs between {time.strftime('%m-%d %H:%M', time.gmtime(min(rec_t or mt)))} and "
                  f"{time.strftime('%m-%d %H:%M', time.gmtime(max(mt)))} UTC ({span_h:.2f} h){extra_t}")
        print(f"      verdict: {verdict}")
print("\nSUMMARY: " + ", ".join(f"{k} {v}" for k, v in sorted(totals.items())))
PYEOF
