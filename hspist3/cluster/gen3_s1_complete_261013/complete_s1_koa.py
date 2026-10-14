#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 47 item 3; the contingency for D's merge, job 15055588): THE REDUCTION COMPLETION OF A STAGE 1
MERGE STOPPED DURING ITS REDUCTION (killed at its time limit after merged/.merge_complete). For the (cell, mass) folders of
<root>/merged/stage1/gen3 that the merge did not finish it does the merge's own post-move steps (cluster/gen3_stage1_261012/merge_s1_koa.py
execute(), after .merge_complete), in the merge's order, with the package's own scripts and commands:
  per cell (sorted):   reduce_B.py <cell> <the masses whose red_nu.csv or acf_runs.npz is missing or unreadable>
  per m_<M> (sorted):  reduce_s1_koa.py <m_M>, if its tables are missing or incomplete
  per m_<M> (sorted):  the SHA-256 lines of its uncompressed traces appended to .sha256_uncompressed (those not there yet), then gzip -9
                       of each uncompressed trace
WHAT IS FINISHED is read from the files, never from the merge log: the merge's python prints are block-buffered (cluster/koa_env.sh sets
no PYTHONUNBUFFERED), so a killed merge's log lacks its last lines. A folder is
  reduced by reduce_B  when red_nu.csv reads with one complete row per trace (no empty field) and acf_runs.npz loads with run<r> for each;
  reduced by reduce_s1 when nud_runs.csv, runs_health.csv, failed_runs.csv, acf_mean_fit.csv and s1_runs.csv read, nud_runs.csv and
                       s1_runs.csv with one complete row per red_nu.csv row, and snap_runs.csv with one row per snapshot file when the
                       folder holds snapshot files;
  compressed           when no trace is left uncompressed, every trace has its .csv.gz and its line in .sha256_uncompressed.
IT CHANGES NO FINISHED FILE: a finished folder is not touched; a finished step of an unfinished folder is not redone (its files are not
rewritten); reduce_B runs only for the masses whose own files are missing or unreadable. A trace whose gzip was cut (its .csv and a .csv.gz
both there) has the partial .csv.gz renamed to <name>.csv.gz.partial_<job id> (kept, never deleted) before gzip -9 runs again; the trace's
SHA-256 line must already be there and equal the .csv's SHA-256 (the merge writes every line of a folder before its first gzip).
IT REFUSES (exit 2, nothing changed) without merged/.merge_complete, with a state it does not recognise (a .sha256_uncompressed that holds
some but not all lines of a folder's uncompressed traces, a line that differs from its trace, a last line without its newline, a .csv.gz
without its line, a compressed trace before reduce_s1 ran ...), and (exit 3) when merged/.reduction_completed exists (it ran before to the
end). Every line it prints is flushed at once. At the end, when every step exited 0 and something was done, it writes
merged/.reduction_completed (a new file: the time and the counts).
usage: python3 cluster/gen3_s1_complete_261013/complete_s1_koa.py <root of one group: .../stage1/<group>> [--plan-only]
exit: 0 completed, or nothing left to do; 2 refused, nothing changed; 3 refused, it ran before; 1 a step failed (no .reduction_completed).
"""
import glob, hashlib, os, re, subprocess, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE)
REDUCE_B = os.path.join(CL, "confinement_20261013", "reduce_B.py")
REDUCE_S1 = os.path.join(CL, "gen3_stage1_261012", "reduce_s1_koa.py")
TRACE = re.compile(r"^wall_x_positions_.*_run(\d+)\.csv$")
SNAP = re.compile(r"^snapshots_.*_run\d+\.bin(\.gz)?$")
SHA_LINE = re.compile(r"^([0-9a-f]{64})  (\S+)$")
S1_TABLES = ("nud_runs.csv", "runs_health.csv", "failed_runs.csv", "acf_mean_fit.csv", "s1_runs.csv")


def say(*a): print(*a, flush=True)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def traces(md):
    """{trace name (.csv): (has .csv, has .csv.gz)}"""
    out = {}
    for f in os.listdir(md):
        b = f[:-3] if f.endswith(".csv.gz") else f
        if TRACE.match(b): out.setdefault(b, [False, False])[1 if f.endswith(".gz") else 0] = True
    return {k: tuple(v) for k, v in sorted(out.items())}


def rows_ok(p, n, cols, runs=None):
    """the table reads, has n rows, no empty field in cols, and (if runs) its run column equals runs"""
    try:
        t = pd.read_csv(p)
    except Exception:
        return False
    if len(t) != n or any(c not in t.columns for c in cols) or t[list(cols)].isna().any().any(): return False
    return runs is None or sorted(t["run"].astype(int)) == sorted(runs)


def b_done(md, tr):
    runs = sorted(int(TRACE.match(t).group(1)) for t in tr)
    if not os.path.exists(os.path.join(md, "red_nu.csv")) or not os.path.exists(os.path.join(md, "acf_runs.npz")): return False, runs
    if not rows_ok(os.path.join(md, "red_nu.csv"), len(runs), ("run", "seed", "nu", "nu_pred", "dt", "n"), runs): return False, runs
    try:
        z = np.load(os.path.join(md, "acf_runs.npz"))
        if sorted(z.files) != sorted(f"run{r}" for r in runs): return False, runs
        for k in z.files: np.asarray(z[k])
    except Exception:
        return False, runs
    return True, runs


def s1_done(md, runs):
    if not all(os.path.exists(os.path.join(md, f)) for f in S1_TABLES): return False
    n = len(runs)
    if not rows_ok(os.path.join(md, "nud_runs.csv"), n, ("run", "seed", "nu_arg", "fit_ok"), runs): return False
    if not rows_ok(os.path.join(md, "s1_runs.csv"), n, ("run", "seed", "trace", "psi6", "virial", "snapshot"), runs): return False
    if not rows_ok(os.path.join(md, "acf_mean_fit.csv"), 1, ("n", "fit_ok")): return False
    for f in ("runs_health.csv", "failed_runs.csv"):
        try: pd.read_csv(os.path.join(md, f))
        except Exception: return False
    ns = sum(1 for f in os.listdir(md) if SNAP.match(f))
    if ns: return rows_ok(os.path.join(md, "snap_runs.csv"), ns, ("run", "file"))
    return True


def sha_lines(md):
    """{name: sha} of the complete lines of .sha256_uncompressed; None if its last line has no newline or a line does not parse"""
    p = os.path.join(md, ".sha256_uncompressed")
    if not os.path.exists(p): return {}
    b = open(p, "rb").read()
    if b and not b.endswith(b"\n"): return None
    out = {}
    for l in b.decode().splitlines():
        m = SHA_LINE.match(l)
        if not m: return None
        out[m.group(2)] = m.group(1)
    return out


def plan_folder(md, rel):
    """(state, details): state in finished / todo; the steps to do"""
    tr = traces(md)
    if not tr: return dict(rel=rel, tr=tr, runs=[], b=True, s1=True, comp=True, note="no traces: nothing to do")
    bd, runs = b_done(md, tr)
    sd = bd and s1_done(md, runs)
    unc = [t for t, (c, z) in tr.items() if c]; gz_only = [t for t, (c, z) in tr.items() if z and not c]; both = [t for t, (c, z) in tr.items() if c and z]
    L = sha_lines(md)
    if L is None: return dict(rel=rel, refuse=f"{rel}/.sha256_uncompressed: a last line without its newline, or a line that does not parse")
    have = [t for t in tr if t in L]
    if (gz_only or both) and not sd: return dict(rel=rel, refuse=f"{rel}: compressed traces before its reduction is complete")
    for t in gz_only:
        if t not in L: return dict(rel=rel, refuse=f"{rel}/{t}.gz: no line in .sha256_uncompressed")
    if have and len([t for t in unc if t in L]) not in (0, len(unc)):
        return dict(rel=rel, refuse=f"{rel}: .sha256_uncompressed holds {len([t for t in unc if t in L])} of the {len(unc)} uncompressed traces' lines")
    if unc and not [t for t in unc if t in L] and (gz_only or both):
        return dict(rel=rel, refuse=f"{rel}: traces compressed but the uncompressed ones have no line")
    comp = not unc
    return dict(rel=rel, tr=tr, runs=runs, b=bd, s1=sd, comp=comp, unc=unc, both=both, lines_needed=bool(unc) and not [t for t in unc if t in L])


def main():
    av = sys.argv[1:]; plan_only = "--plan-only" in av; args = [a for a in av if not a.startswith("--")]
    if len(args) != 1: say("usage: complete_s1_koa.py <root of one group> [--plan-only]"); sys.exit(2)
    root = os.path.abspath(args[0]); M = os.path.join(root, "merged"); base = os.path.join(M, "stage1", "gen3")
    job = os.environ.get("SLURM_JOB_ID", "nojob")
    rel = lambda p: os.path.relpath(p, M)
    if not os.path.exists(os.path.join(M, ".merge_complete")):
        say(f"STOP: {M} has no .merge_complete: the merge did not finish its moves; this job completes only a reduction; NOTHING was changed"); sys.exit(2)
    if os.path.exists(os.path.join(M, ".reduction_completed")):
        say(f"STOP: {M}/.reduction_completed exists: the completion ran before to the end; NOTHING was changed"); sys.exit(3)
    cells = [c for c in sorted(glob.glob(os.path.join(base, "*"))) if glob.glob(os.path.join(c, "m_*"))]
    plans = {c: [plan_folder(md, rel(md)) for md in sorted(glob.glob(os.path.join(c, "m_*")))] for c in cells}
    bad = [p["refuse"] for ps in plans.values() for p in ps if "refuse" in p]
    say(f"## the plan ({'--plan-only' if plan_only else 'job ' + job}): {len(cells)} cells, {sum(len(v) for v in plans.values())} (cell, mass) folders under {rel(base)}\n")
    say("| folder | traces | reduce_B | reduce_s1 | compression | to do |\n|---|---|---|---|---|---|")
    todo = 0
    for c, ps in plans.items():
        for p in ps:
            if "refuse" in p: say(f"| {p['rel']} | - | - | - | - | REFUSED: {p['refuse']} |"); continue
            steps = ([] if p["b"] else ["reduce_B"]) + ([] if p["s1"] else ["reduce_s1"]) + ([] if p["comp"] else
                     [("the SHA-256 lines and " if p["lines_needed"] else "") + f"gzip of {len(p['unc'])}" + (f" ({len(p['both'])} partial .gz moved aside)" if p["both"] else "")])
            todo += bool(steps)
            nz = sum(1 for c_, z in p["tr"].values() if z); ntr = len(p["tr"])
            cs = "done" if p["comp"] else f"{nz} of {ntr} compressed"
            say(f"| {p['rel']} | {ntr} | {'done' if p['b'] else 'TO DO'} | {'done' if p['s1'] else 'TO DO'} | {cs} | {', '.join(steps) or 'nothing (finished)'} |")
    if bad:
        say(f"\nSTOP: {len(bad)} folder(s) in a state this job does not recognise (listed above); NOTHING was changed"); sys.exit(2)
    say(f"\nfolders finished: {sum(len(v) for v in plans.values()) - todo}; folders to complete: {todo}")
    if plan_only: say("--plan-only: nothing was changed"); return
    if not todo: say("nothing to complete: every folder is finished; NOTHING was changed"); return
    n = dict(b_cells=0, b_masses=0, s1=0, comp=0, traces=0, partial=0, failed=0)
    for c, ps in plans.items():
        if all(p["b"] and p["s1"] and p["comp"] for p in ps): continue
        need_b = [os.path.basename(os.path.join(c, p["rel"]))[2:] for p in ps if not p["b"]]
        if need_b:
            r = subprocess.run([sys.executable, REDUCE_B, c] + need_b, capture_output=True, text=True)
            say(f"reduce_B {rel(c)} (masses {' '.join(need_b)}): exit {r.returncode} {r.stdout.strip()} {r.stderr.strip()[-300:]}")
            n["b_cells"] += 1; n["b_masses"] += len(need_b)
            if r.returncode != 0: n["failed"] += 1; continue                       # as the merge: no reduce_s1, no compression for this cell
        for p in ps:
            if p["s1"]: continue
            md = os.path.join(M, p["rel"])
            q = subprocess.run([sys.executable, REDUCE_S1, md], capture_output=True, text=True)
            say(f"reduce_s1 {p['rel']}: exit {q.returncode} {q.stdout.strip()} {q.stderr.strip()[-300:]}")
            n["s1"] += 1; n["failed"] += q.returncode != 0
        for p in ps:
            if p["comp"]: continue
            md = os.path.join(M, p["rel"]); tr = traces(md)
            unc = sorted(t for t, (cc, z) in tr.items() if cc)
            L = sha_lines(md) or {}
            if not [t for t in unc if t in L]:                                    # the merge's sha step had not run: run it as the merge does
                with open(os.path.join(md, ".sha256_uncompressed"), "a") as fh:
                    for x in unc: fh.write(f"{sha(os.path.join(md, x))}  {x}\n")
            else:
                for x in unc:
                    if L.get(x) != sha(os.path.join(md, x)): say(f"STOP: {p['rel']}/{x}: its SHA-256 differs from its line; the rest is NOT done"); sys.exit(1)
            for x in unc:
                if os.path.exists(os.path.join(md, x + ".gz")):
                    os.rename(os.path.join(md, x + ".gz"), os.path.join(md, f"{x}.gz.partial_{job}")); n["partial"] += 1
                    say(f"  {p['rel']}/{x}.gz: a cut gzip, renamed to {x}.gz.partial_{job} (kept)")
                subprocess.run(["gzip", "-9", os.path.join(md, x)], check=True)
            say(f"  {p['rel']}: {len(unc)} traces compressed")
            n["comp"] += 1; n["traces"] += len(unc)
    line = (f"reduction completed by job {job}: reduce_B in {n['b_cells']} cells ({n['b_masses']} masses), reduce_s1 in {n['s1']} folders, "
            f"{n['comp']} folders compressed ({n['traces']} traces), {n['partial']} cut gzip files kept aside, {n['failed']} steps with an exit "
            f"other than 0")
    say(line)
    if n["failed"]: say("STOP: a step failed; merged/.reduction_completed is NOT written; inspect by hand"); sys.exit(1)
    with open(os.path.join(M, ".reduction_completed"), "w") as fh: fh.write(time.strftime("%Y-%m-%d %H:%M:%S %Z") + "\n" + line + "\n")


if __name__ == "__main__":
    main()
