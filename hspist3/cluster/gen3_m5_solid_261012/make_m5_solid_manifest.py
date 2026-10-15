#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (decision 38 item 6): writes m5_solid_inputs.sha256, the manifest m5_solid_koa.sbatch checks with sha256sum -c
(paths relative to hspist3/): the runner, every module of this repository that the runner's code paths import (found by importing them
in a fresh interpreter: the runner, kr_cs's plot_speed_of_sound_edmd, argmax_nu's paper1_populate_cs_err_20261002), koa_env.sh and the
binaries' hash file k2_binaries_a643bc7.sha256.
usage (from hspist3/): python3 cluster/gen3_m5_solid_261012/make_m5_solid_manifest.py"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
PROBE = ("import sys, os; sys.path.insert(0, %r); import m5_solid_koa as R; R.T.kr_cs(0.74); "
         "import paper1_populate_cs_err_20261002; "
         "print('\\n'.join(sorted(m.__file__ for m in list(sys.modules.values()) if getattr(m, '__file__', None) and "
         "os.path.abspath(m.__file__).startswith(%r))))") % (HERE, HS + os.sep)
out = subprocess.run([sys.executable, "-c", PROBE], capture_output=True, text=True, check=True, cwd=HS).stdout.split()
files = sorted({os.path.relpath(os.path.abspath(p), HS) for p in out} | {"cluster/koa_env.sh", "cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256"})
with open(os.path.join(HERE, "m5_solid_inputs.sha256"), "w") as fh:
    for f in files: fh.write(f"{hashlib.sha256(open(os.path.join(HS, f), 'rb').read()).hexdigest()}  {f}\n")
print(f"m5_solid_inputs.sha256: {len(files)} files"); [print(f"  {f}") for f in files]
