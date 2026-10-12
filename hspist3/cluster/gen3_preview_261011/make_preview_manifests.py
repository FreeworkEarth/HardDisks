#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): the two input manifests of the PREVIEW jobs. THIS FILE IS
cluster/gen3_testR_261011/make_testR_manifests.py (b5378f4) with the preview's file names; the merge files are Test R's, run unchanged.
Format: sha256sum -c, paths relative to hspist3/. They are written before each package commit, ON THE DELIVERED TREE (the engine-gen3
worktree; Report R0's lesson: cluster/confinement_20261013/reduce_B.py differs between main and engine-gen3), and the jobs check them
at their start (no git at run time):
  preview_inputs.sha256        what every chunk uses: the task list, its SHA-256 file, the chunk size, the build record (Test R's
                               cluster/gen3_testR_261011/build_record_a643bc7.txt), the worker cluster/gen3_koa_261009/testG_worker_koa.sh
                               (unchanged since 1cdb7ff6) and cluster/koa_env.sh;
  preview_merge_inputs.sha256  what the merge job runs: Test R's chunk_logs_koa.sh, merge_koa.py and reduce_testR_koa.py, reduce_B.py, and
                               every module under hspist3/ that reduce_B.py and reduce_testR_koa.py import (found by importing them here).
usage (from hspist3/ of the engine-gen3 worktree): python3 cluster/gen3_preview_261011/make_preview_manifests.py
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
P = os.path.relpath(HERE, HS); R = "cluster/gen3_testR_261011"


def sha(p): return hashlib.sha256(open(os.path.join(HS, p), "rb").read()).hexdigest()


def imported():
    code = ("import sys, os; sys.path.insert(0, 'validation'); sys.path.insert(0, '.'); sys.argv = ['x'];\n"
            "import importlib.util as u\n"
            f"for n, p in (('reduce_testR_koa', '{R}/reduce_testR_koa.py'),):\n"
            "    s = u.spec_from_file_location(n, p); m = u.module_from_spec(s); s.loader.exec_module(m)\n"
            "import tests_20260913, paper1_populate_cs_err_20261002\n"
            "hs = os.path.abspath('.')\n"
            "print('\\n'.join(sorted({os.path.relpath(m.__file__, hs) for m in list(sys.modules.values())\n"
            "      if getattr(m, '__file__', None) and os.path.abspath(m.__file__).startswith(hs + os.sep)})))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=HS, capture_output=True, text=True)
    if out.returncode: sys.exit(f"STOP: the import probe failed: {out.stderr[-600:]}")
    return [l for l in out.stdout.splitlines() if l.endswith(".py")]


def main():
    a = [f"{P}/tasks_preview_koa_261011.txt", f"{P}/tasks_preview_koa_261011.sha256", f"{P}/preview_chunk_size.txt",
         f"{R}/build_record_a643bc7.txt", "cluster/gen3_koa_261009/testG_worker_koa.sh", "cluster/koa_env.sh"]
    b = [f"{R}/chunk_logs_koa.sh", f"{R}/merge_koa.py", f"{R}/reduce_testR_koa.py", "cluster/confinement_20261013/reduce_B.py"]
    b += [m for m in imported() if m not in b]
    for name, files in (("preview_inputs.sha256", a), ("preview_merge_inputs.sha256", b)):
        txt = "".join(f"{sha(f)}  {f}\n" for f in files)
        open(os.path.join(HERE, name), "w").write(txt)
        print(f"== {P}/{name} ({len(files)} files)\n{txt}")


if __name__ == "__main__":
    main()
