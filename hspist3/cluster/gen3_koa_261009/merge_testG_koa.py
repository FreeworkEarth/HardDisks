#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): the merge of TEST G ON KOA's chunks, one process (no two writers): every chunk
folder chunks/chunk_<j> holds the Mac layout for its own trajectories (testG/<engine>/<cell>/m_<M>/ and .../x_0/). The merge MOVES each
chunk's files into merged/testG/..., appends each chunk's run.log ##RUN sections to the merged cell's run.log in chunk order, keeps
every .failed_* folder (renamed with its chunk), checks that every chunk's .build_git is the same line, then reduces every B cell with
reduce_B.py (unchanged) and compresses its traces (gzip -9, the SHA-256 of each uncompressed trace in .sha256_uncompressed first), as
run_testG_mac.py does on the Mac. Nothing is deleted: a chunk folder that is empty after the merge stays.
##CHRIS 2026-10-10 (decision 20; sec. 4.7.47): the A-fixed runs' 00_COMMAND.md and run_params.json are kept per chunk
(<name>_chunk<j>), the fault the Mac rehearsal found.
##CHRIS 2026-10-10 (plan-author decision 21, item 2; 261012 sec. 4.7.49): THE MERGE PLANS BEFORE IT MOVES.
  1. SCAN (reads only): every chunk, every cell folder, every name, against the RULES below. The complete plan (treatment, chunk,
     source, destination) is written to <root>/merge_plan_<time>.tsv and summarised: each destination folder with the chunks that
     deliver into it, every name that occurs in more than one chunk for one destination with its treatment, and the top-level names
     inside the moved run folders (.done_runs/run<r>; one line: the rehearsal of 882dba4 listed every file there, and the driver's
     experiments_speed_of_sound/.../simulation_<time>/ folders made that line grow with the number of runs).
  2. CHECK: a name without a rule, an entry outside the layout chunks/chunk_<j>/testG/<gen2|gen3>/<cell>/<m_M | x_0>/, two plan entries
     with one destination, a cell folder without .build_git, and .build_git files that differ are ALL printed; then the merge stops
     with exit 2 and NOTHING moved (merged/ is not created), and it may be started again once the cause is removed.
  3. EXECUTE: only a plan without problems; merged/ is created, the plan is carried out in chunk order, and merged/.merge_complete
     records the end of the moves. A start is refused (exit 3) when merged/ exists: with .merge_complete the merge was completed; without
     it a merge stopped during its execution and is inspected by hand. reduce_B and the compression follow, as before.
RULES (per cell folder testG/<gen2|gen3>/<cell>/<m_M | x_0>; B = the speed-of-sound folders m_<M>, AF = the A-fixed folder x_0):
  EQUAL       .build_git -- every chunk's copy identical; the merged cell gets that line once
  APPEND      run.log (B), .sha256_uncompressed (B, AF) -- the chunks' files appended to the merged file in chunk order, byte for byte
  MOVE        B: wall_x_positions_*_run<r>.csv[.gz], psi6_t_*_run<r>.csv[.gz]; AF: ev_, red_, run_, summary_, tr_, stream_<seed>.csv /
              .csv.gz / .log -- the run index or the seed is in the name; the same name from two chunks stops the merge
  PER CHUNK   B: .done_runs/run<r> -> .done_runs/run<r>_chunk<j> (each run folder holds that run's 00_COMMAND.md, 01_PLOT_COMMANDS.md,
              00_COMMAND_AND_PLOT_COMMAND.md, speed_of_sound_*, stdout.log ...: the binary runs inside it); .failed_*, .stale_* ->
              <name>_chunk<j>; AF: 00_COMMAND.md, run_params.json (every A-fixed run of a chunk rewrites them in the shared x_0, so each
              chunk holds its last run's copy) -> <base>_chunk<j><ext>; .failed_*, .stale_*, .fifo_used -> <name>_chunk<j>
  INCOMPLETE  B: .run<r> (a run folder an interrupted worker left) -> .incomplete_run<r>_chunk<j>, printed
  LEFT        .guard.lock, .runlog.lockdir, .sha.lockdir, .build_git.tmp<pid> (an interrupted worker's lock or temp file) -- left in
              the chunk folder, printed
usage (inside merge_testG_k2.sbatch): python3 cluster/gen3_koa_261009/merge_testG_koa.py <root: .../testG_koa> [--plan-only] [--no-reduce]
  --plan-only   scan and check, print and write the plan, move nothing (exit 0 for a clean plan, 2 for a plan with problems)
  --no-reduce   skip reduce_B and the compression (the unit test's synthetic chunks hold no real traces)
exit: 0 merged; 2 the plan has problems, NOTHING moved; 3 refused at the start (merged/ exists); 1 an error during the execution.
"""
import glob, hashlib, os, re, shutil, subprocess, sys, time, traceback
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE)
B_MOVE = re.compile(r"^(wall_x_positions_.*_run\d+\.csv|psi6_t_.*_run\d+\.csv)(\.gz)?$")
AF_MOVE = re.compile(r"^(ev|red|run|summary|tr|stream)_\d+\.(csv|csv\.gz|log)$")
LEFT = re.compile(r"^(\.guard\.lock|\.runlog\.lockdir|\.sha\.lockdir|\.build_git\.tmp\d+)$")
WHAT = {"EQUAL": "keep one after an equality check", "APPEND": "append in chunk order", "MOVE": "move (unique name)",
        "PER CHUNK": "rename per chunk", "INCOMPLETE": "rename per chunk, reported", "LEFT": "left in the chunk folder, reported"}


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def rule(kind, name, isdir):
    """(treatment, destination name with {j} for the chunk) for one name in a cell folder of kind 'B' or 'AF'; (None, why) without a rule"""
    if name == ".build_git" and not isdir: return "EQUAL", name
    if LEFT.match(name): return "LEFT", None
    if name == ".sha256_uncompressed" and not isdir: return "APPEND", name
    if name.startswith((".failed_", ".stale_")) and isdir: return "PER CHUNK", name + "_chunk{j}"
    if kind == "B":
        if name == "run.log" and not isdir: return "APPEND", name
        if name == ".done_runs" and isdir: return "DONE RUNS", None
        if re.match(r"^\.run\d+$", name) and isdir: return "INCOMPLETE", ".incomplete_" + name[1:] + "_chunk{j}"
        if B_MOVE.match(name) and not isdir: return "MOVE", name
    else:
        if name in ("00_COMMAND.md", "run_params.json") and not isdir:
            b, e = os.path.splitext(name); return "PER CHUNK", b + "_chunk{j}" + e
        if name == ".fifo_used" and isdir: return "PER CHUNK", name + "_chunk{j}"
        if AF_MOVE.match(name) and not isdir: return "MOVE", name
    return None, f"no rule for a {'folder' if isdir else 'file'} named '{name}' in {'a speed-of-sound (B)' if kind == 'B' else 'the A-fixed'} cell folder"


def scan(root):
    """the plan, the problems, the build lines and the names inside the run folders; reads only"""
    plan, problems, builds, inside = [], [], defaultdict(list), defaultdict(int)
    R = lambda p: os.path.relpath(p, root)
    kind_of = lambda p: "folder" if os.path.isdir(p) else "file"
    cdir = os.path.join(root, "chunks")
    if not os.path.isdir(cdir): return [], plan, [f"no folder {cdir}"], builds, inside
    chunks = []
    for n in sorted(os.listdir(cdir)):
        p = os.path.join(cdir, n)
        if re.match(r"^chunk_\d+$", n) and os.path.isdir(p): chunks.append(p)
        else: problems.append(f"unexpected {kind_of(p)} {R(p)} (the layout is chunks/chunk_<j>/)")
    chunks.sort(key=lambda p: int(p.rsplit("_", 1)[1]))
    for ch in chunks:
        j = ch.rsplit("_", 1)[1]
        for n in sorted(os.listdir(ch)):
            if n != "testG" or not os.path.isdir(os.path.join(ch, n)): problems.append(f"unexpected {kind_of(os.path.join(ch, n))} {R(os.path.join(ch, n))}")
        tg = os.path.join(ch, "testG")
        cells = []
        for e in (sorted(os.listdir(tg)) if os.path.isdir(tg) else []):
            pe = os.path.join(tg, e)
            if e not in ("gen2", "gen3") or not os.path.isdir(pe): problems.append(f"unexpected {kind_of(pe)} {R(pe)} (engines: gen2, gen3)"); continue
            for c in sorted(os.listdir(pe)):
                pc = os.path.join(pe, c)
                if not os.path.isdir(pc): problems.append(f"unexpected file {R(pc)}"); continue
                for k in sorted(os.listdir(pc)):
                    pk = os.path.join(pc, k)
                    if os.path.isdir(pk) and (re.match(r"^m_\d+$", k) or k == "x_0"): cells.append(pk)
                    else: problems.append(f"unexpected {kind_of(pk)} {R(pk)} (cell folders: m_<M>, x_0)")
        for cell in cells:
            rel = os.path.relpath(cell, ch); kind = "AF" if os.path.basename(cell) == "x_0" else "B"
            names = sorted(os.listdir(cell))
            if ".build_git" not in names and names: problems.append(f"{R(cell)}: entries without a .build_git")
            for name in names:
                src = os.path.join(cell, name); isdir = os.path.isdir(src) and not os.path.islink(src)
                t, dname = rule(kind, name, isdir)
                if t is None: problems.append(f"{R(src)}: {dname}"); continue
                if t == "EQUAL":
                    builds[open(src, errors="replace").read().strip()].append(R(src)); plan.append((t, src, os.path.join(rel, name), j, name, rel)); continue
                if t == "LEFT": plan.append((t, src, None, j, name, rel)); continue
                if t == "DONE RUNS":
                    for r in sorted(os.listdir(src)):
                        pr = os.path.join(src, r)
                        if not re.match(r"^run\d+$", r) or not os.path.isdir(pr):
                            problems.append(f"{R(pr)}: no rule (an entry of .done_runs must be a folder run<r>)"); continue
                        for e in os.listdir(pr):                        # the top level only: deeper names carry run times
                            inside[e + ("/" if os.path.isdir(os.path.join(pr, e)) else "")] += 1
                        plan.append(("PER CHUNK", pr, os.path.join(rel, ".done_runs", f"{r}_chunk{j}"), j, ".done_runs/run<r>", rel))
                    continue
                plan.append((t, src, os.path.join(rel, dname.replace("{j}", j)), j, name, rel))
    if len(builds) > 1:
        for line, where in sorted(builds.items()):
            problems.append(f".build_git differs: '{line}' in {len(where)} cell folder(s), e.g. {where[0]}")
    if not builds and plan: problems.append("no .build_git in any cell folder")
    seen = {}
    for t, src, dst, j, name, cell in plan:                        # one destination per entry (EQUAL and APPEND share theirs by rule)
        if dst is None or t in ("EQUAL", "APPEND"): continue
        if dst in seen: problems.append(f"two entries for one destination {dst}: {R(seen[dst])} and {R(src)}")
        else: seen[dst] = src
    return chunks, plan, problems, builds, inside


def summary(plan, inside):
    by_dst, names = defaultdict(set), defaultdict(lambda: defaultdict(set))
    for t, src, dst, j, name, cell in plan:
        by_dst[cell].add(int(j)); names[cell][(name, t)].add(int(j))
    gen = lambda n: re.sub(r"_run\d+", "_run<r>", re.sub(r"_\d+(\.csv|\.log)", r"_<seed>\1", n))
    print("## the merge plan: every destination folder, the chunks that deliver into it, and the names found in more than one chunk "
          "for that destination\n")
    print("| destination (under merged/) | chunks | names in more than one chunk: treatment |\n|---|---|---|")
    for cell in sorted(by_dst):
        js = sorted(by_dst[cell])
        span = f"{js[0]}-{js[-1]} ({len(js)})" if len(js) > 3 and js == list(range(js[0], js[-1] + 1)) else ", ".join(map(str, js))
        multi = sorted({f"{gen(n)}: {WHAT[t]}" for (n, t), cj in names[cell].items() if len(cj) > 1})
        print(f"| {cell} | {span} | {'; '.join(multi) if multi else '-'} |")
    cnt = defaultdict(int)
    for t, *_ in plan: cnt[t] += 1
    print("\nplan entries by treatment: " + ", ".join(f"{k} {v}" for k, v in sorted(cnt.items())))
    if inside:
        print("names inside the moved run folders .done_runs/run<r> (each run's own folder; they never meet at the cell level): "
              + ", ".join(f"{n} x{c}" for n, c in sorted(inside.items())))


def execute(root, plan, builds, no_reduce):
    M = os.path.join(root, "merged"); os.makedirs(M)
    n = defaultdict(int)
    for t, src, dst, j, name, cell in plan:
        n[t] += 1
        if t in ("LEFT", "EQUAL"): continue
        out = os.path.join(M, dst); os.makedirs(os.path.dirname(out), exist_ok=True)
        if t == "APPEND":
            with open(out, "ab") as fo, open(src, "rb") as fi: shutil.copyfileobj(fi, fo)
            continue
        if os.path.lexists(out): raise RuntimeError(f"{out} exists although the plan had one entry for it")
        shutil.move(src, out)
    b = next(iter(builds)) if builds else ""
    for cell in glob.glob(os.path.join(M, "testG", "*", "*", "*")):
        with open(os.path.join(cell, ".build_git"), "w") as fh: fh.write(b + "\n")
    failed = sum(1 for e in plan if e[4].startswith(".failed_"))
    line = (f"merged {len(set(e[3] for e in plan))} chunks: {n['MOVE']} files moved, {n['APPEND']} files appended, {n['PER CHUNK']} kept per "
            f"chunk ({failed} failed-run folders), {n['INCOMPLETE']} incomplete run folders, {n['LEFT']} locks or temp files left in place, "
            f"{n['EQUAL']} build files equal; build '{b}'")
    with open(os.path.join(M, ".merge_complete"), "w") as fh: fh.write(time.strftime("%Y-%m-%d %H:%M:%S %Z") + "\n" + line + "\n")
    print(line)
    if no_reduce: return
    for c in sorted(glob.glob(os.path.join(M, "testG", "*", "*"))):
        if not glob.glob(os.path.join(c, "m_*")): continue
        r = subprocess.run([sys.executable, os.path.join(CL, "confinement_20261013", "reduce_B.py"), c], capture_output=True, text=True)
        print(f"reduce_B {os.path.relpath(c, M)}: exit {r.returncode} {r.stdout.strip()} {r.stderr.strip()[-300:]}")
        if r.returncode != 0: continue
        for md in sorted(glob.glob(os.path.join(c, "m_*"))):
            tr = sorted(x for x in os.listdir(md) if x.startswith("wall_x_positions_") and x.endswith(".csv"))
            with open(os.path.join(md, ".sha256_uncompressed"), "a") as fh:
                for x in tr: fh.write(f"{sha(os.path.join(md, x))}  {x}\n")
            for x in tr: subprocess.run(["gzip", "-9", os.path.join(md, x)], check=True)
            print(f"  {os.path.relpath(md, M)}: {len(tr)} traces compressed")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1: print("usage: merge_testG_koa.py <root> [--plan-only] [--no-reduce]"); sys.exit(3)
    root = os.path.abspath(args[0]); M = os.path.join(root, "merged")
    if os.path.exists(os.path.join(M, ".merge_complete")): print(f"STOP: {M} holds a completed merge (.merge_complete); a merge is never repeated"); sys.exit(3)
    if os.path.lexists(M): print(f"STOP: {M} exists without .merge_complete: a merge stopped during its execution; inspect by hand"); sys.exit(3)
    try:
        chunks, plan, problems, builds, inside = scan(root)
        stamp = time.strftime("%Y%m%d_%H%M%S")
        with open(os.path.join(root, f"merge_plan_{stamp}.tsv"), "w") as fh:
            fh.write("treatment\tchunk\tsource (under the root)\tdestination (under merged/)\n")
            for t, src, dst, j, name, cell in plan: fh.write(f"{t}\t{j}\t{os.path.relpath(src, root)}\t{dst or '(left in place)'}\n")
    except Exception:
        traceback.print_exc(); print("STOP: the scan failed; NOTHING moved"); sys.exit(2)
    print(f"scan: {len(chunks)} chunks, {len(plan)} plan entries, {len(problems)} problem(s); plan written to {root}/merge_plan_{stamp}.tsv\n")
    summary(plan, inside)
    if problems:
        print(f"\n## {len(problems)} problem(s): NOTHING MOVED\n")
        for p in problems: print(f"  - {p}")
        print("STOP: the merge plan has problems; nothing was moved (remove the cause, then start the merge again)")
        sys.exit(2)
    odd = [e for e in plan if e[0] in ("LEFT", "INCOMPLETE")]
    if odd: print(f"\nnote, {len(odd)} lock / temp / incomplete entries: " + "; ".join(f"{e[0]} {os.path.relpath(e[1], root)}" for e in odd[:20])
                  + (f"; ... and {len(odd) - 20} more (all in the plan file)" if len(odd) > 20 else ""))
    if "--plan-only" in sys.argv: print("\n--plan-only: the plan has no problems; nothing moved"); return
    print()
    try:
        execute(root, plan, builds, "--no-reduce" in sys.argv)
    except Exception:
        traceback.print_exc(); print(f"STOP: the merge failed during its execution; {M} has no .merge_complete; inspect by hand"); sys.exit(1)


if __name__ == "__main__":
    main()
