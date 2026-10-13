#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 c): the input manifest of the stage 1 gate jobs, s1_gate_inputs.sha256 (sha256sum -c
format, paths relative to hspist3/), as make_testR_manifests.py: every file the gate job, the cross-node tasks and the compare job run or
read from the scripts clone -- the gate package's own files, Test G's worker and task list, the a643bc7 binaries' hash file, the cross-node
compare, cluster/resched_gate_261005/audit_runs_261007.py and every module under hspist3/ it and the compares import (found by importing
them in a fresh interpreter), and cluster/koa_env.sh.
usage (from hspist3/): python3 cluster/gen3_s1_gate_261012/make_s1_gate_manifest.py
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
P = os.path.relpath(HERE, HS)


def sha(p): return hashlib.sha256(open(os.path.join(HS, p), "rb").read()).hexdigest()


def imported():
    code = ("import sys, os; sys.argv = ['x']; import importlib.util as u\n"
            "for n, p in (('audit_runs_261007', 'cluster/resched_gate_261005/audit_runs_261007.py'), ('g3_xnode_compare', 'cluster/gen3_koa_261009/g3_xnode_compare.py'),\n"
            f"             ('s1_gate_testg_compare', '{P}/s1_gate_testg_compare.py')):\n"
            "    sys.path.insert(0, os.path.dirname(os.path.abspath(p))); s = u.spec_from_file_location(n, p); m = u.module_from_spec(s); s.loader.exec_module(m)\n"
            "hs = os.path.abspath('.')\n"
            "print('\\n'.join(sorted({os.path.relpath(m.__file__, hs) for m in list(sys.modules.values())\n"
            "      if getattr(m, '__file__', None) and os.path.abspath(m.__file__).startswith(hs + os.sep)})))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=HS, capture_output=True, text=True)
    if out.returncode: sys.exit(f"STOP: the import probe failed: {out.stderr[-800:]}")
    return [l for l in out.stdout.splitlines() if l.endswith(".py")]


def main():
    a = [f"{P}/s1_gate_testg_compare.py", f"{P}/koa_harness_ref_a643bc7.sha256", f"{P}/record_s1_build.sh",
         "cluster/resched_gate_261005/audit_runs_261007.py", "cluster/gen3_koa_261009/testG_worker_koa.sh",
         "cluster/gen3_koa_261009/tasks_testG_koa_261009.txt", "cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256",
         "cluster/gen3_koa_261009/g3_xnode_compare.py", "cluster/koa_env.sh"]
    a += [m for m in imported() if m not in a]
    txt = "".join(f"{sha(f)}  {f}\n" for f in a)
    open(os.path.join(HERE, "s1_gate_inputs.sha256"), "w").write(txt)
    print(f"== {P}/s1_gate_inputs.sha256 ({len(a)} files)\n{txt}")


if __name__ == "__main__":
    main()
