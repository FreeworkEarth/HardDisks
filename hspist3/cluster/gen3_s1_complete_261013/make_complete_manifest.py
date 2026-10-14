#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 47 item 3): WRITES s1_complete_inputs.sha256, the manifest of every file the reduction
completion runs: its own two files, and the stage 1 package's files of the reduction (the lines of s1_merge_inputs.sha256 and koa_env.sh's
line of s1_inputs.sha256, copied unchanged after checking that the files here still have those SHA-256). Run from hspist3/.
usage: python3 cluster/gen3_s1_complete_261013/make_complete_manifest.py
"""
import hashlib, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
S1 = os.path.join(HS, "cluster", "gen3_stage1_261012")
OWN = ["cluster/gen3_s1_complete_261013/complete_s1_koa.py", "cluster/gen3_s1_complete_261013/s1_complete.sbatch"]


def sha(p): return hashlib.sha256(open(os.path.join(HS, p), "rb").read()).hexdigest()


lines = [f"{sha(p)}  {p}" for p in OWN]
pkg = [l.rstrip("\n") for l in open(os.path.join(S1, "s1_merge_inputs.sha256")) if l.strip()]
pkg += [l.rstrip("\n") for l in open(os.path.join(S1, "s1_inputs.sha256")) if l.strip().endswith("  cluster/koa_env.sh")]
for l in pkg:
    h, p = l.split("  ", 1)
    if sha(p) != h: sys.exit(f"STOP: {p} differs from its line in the package's manifests")
    lines.append(l)
open(os.path.join(HERE, "s1_complete_inputs.sha256"), "w").write("\n".join(lines) + "\n")
print(f"s1_complete_inputs.sha256: {len(lines)} files ({len(OWN)} of the completion, {len(pkg)} of the stage 1 package)")
