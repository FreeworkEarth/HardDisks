#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. The merge of TEST G ON KOA's chunks, one process (no two
writers): every chunk folder chunks/chunk_<j> holds the Mac layout for its own trajectories (testG/<engine>/<cell>/m_<M>/ and .../x_0/).
The merge MOVES each chunk's files into merged/testG/... (names are unique: the run index or the seed is in every name), appends each
chunk's run.log ##RUN sections to the merged cell's run.log in chunk order, keeps every .failed_* folder (renamed with its chunk),
checks that every chunk's .build_git is the same line, then reduces every B cell with reduce_B.py (unchanged) and compresses its
traces (gzip -9, the SHA-256 of each uncompressed trace in .sha256_uncompressed first), as run_testG_mac.py does on the Mac.
Nothing is deleted: a chunk folder that is empty after the merge stays; a name collision stops the merge.
usage (inside merge_testG_koa.sbatch): python3 cluster/gen3_koa_261009/merge_testG_koa.py <root: .../testG_koa>
"""
import glob, hashlib, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    root = os.path.abspath(sys.argv[1]); M = os.path.join(root, "merged")
    chunks = sorted(glob.glob(os.path.join(root, "chunks", "chunk_*")), key=lambda p: int(p.rsplit("_", 1)[1]))
    builds, moved, logs, failed = set(), 0, 0, 0
    for ch in chunks:
        j = ch.rsplit("_", 1)[1]
        for cell in sorted(glob.glob(os.path.join(ch, "testG", "*", "*", "*"))):
            rel = os.path.relpath(cell, ch); dst = os.path.join(M, rel); os.makedirs(dst, exist_ok=True)
            bg = os.path.join(cell, ".build_git")
            if os.path.exists(bg): builds.add(open(bg).read().strip())
            for name in sorted(os.listdir(cell)):
                src = os.path.join(cell, name)
                if name in (".build_git", ".guard.lock", ".runlog.lockdir", ".sha.lockdir"): continue
                if name == "run.log":
                    with open(os.path.join(dst, "run.log"), "a") as out, open(src) as inp: out.write(inp.read())
                    logs += 1; continue
                if name in (".sha256_uncompressed",):
                    with open(os.path.join(dst, name), "a") as out, open(src) as inp: out.write(inp.read())
                    continue
                if name.startswith(".failed_"):
                    shutil.move(src, os.path.join(dst, f"{name}_chunk{j}")); failed += 1; continue
                if name in (".done_runs",):
                    os.makedirs(os.path.join(dst, ".done_runs"), exist_ok=True)
                    for r in os.listdir(src): shutil.move(os.path.join(src, r), os.path.join(dst, ".done_runs", f"{r}_chunk{j}"))
                    continue
                if os.path.exists(os.path.join(dst, name)): sys.exit(f"STOP: name collision {os.path.join(rel, name)} (chunk {j})")
                shutil.move(src, os.path.join(dst, name)); moved += 1
    if len(builds) != 1: sys.exit(f"STOP: the chunks were written by {len(builds)} builds: {sorted(builds)}")
    b = sorted(builds)[0] if builds else ""
    for cell in glob.glob(os.path.join(M, "testG", "*", "*", "*")):
        with open(os.path.join(cell, ".build_git"), "w") as fh: fh.write(b + "\n")
    print(f"merged {len(chunks)} chunks: {moved} files moved, {logs} run logs appended, {failed} failed-run folders kept; build '{b}'")
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


if __name__ == "__main__":
    main()
