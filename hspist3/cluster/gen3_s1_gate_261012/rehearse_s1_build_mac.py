#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 a; decision 20's rule): THE MAC REHEARSAL OF GATE BLOCK B's BUILD LINES, the ones
rehearse_s1_gate_mac.py does not run (it builds its stand-ins with rehearse_k2_mac.build()): build_s1_koa.sh in a FRESH clone of f72d764,
then record_s1_build.sh, run from the scripts clone at --commit, as in the runsheet. Nothing outside --root is written.
  home/harddisks_gen3_s1   f72d764, a shared clone with KOA's step-1 sparse set and the pinned kissfft (rehearse_k2_mac.clone/kissfft)
  home/harddisks_gen3_k2   --commit with the K2 clone's sparse set
The steps: f72d764's own build_gen3_koa.sh refuses the clone (the reason build_s1_koa.sh exists); build_s1_koa.sh builds and prints BUILD OK;
record_s1_build.sh writes the record; a second build and a second record are refused (never rebuilt); the binaries are compared, for
information, with --standin (the f72d764 stand-in of the gate rehearsal, built by rehearse_k2_mac.build() with the same Mac flags).
SHIMS (each a listed residual difference): SLURM_JOB_ID set; module (nothing), gcc (--version prints 14.3.0; otherwise cc without the x86
flags -march=x86-64-v2 and -mtune=generic), make (the real make with rehearse_k2_mac.build()'s Mac CFLAGS_KOA and BUILD_CF), ldd (otool -L);
HOME; build_s1_koa.sh's /tmp/ld_$$ (the long-double report) rewritten to <root>/tmp. On arm64 the long double is the double.
usage (from hspist3/): python3 cluster/gen3_s1_gate_261012/rehearse_s1_build_mac.py --root R --commit C [--standin <home/harddisks_gen3_s1>]
"""
import argparse, os, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009"))
import rehearse_k2_mac as K2R
F72 = "f72d764a8a7e837644c26bd8abc175d9e2cac268"; P = "cluster/gen3_s1_gate_261012"
BINS = ("00ALLINONE", "bin/00ALLINONE_ld", "bin/gen3_m1", "bin/gen3_m2")
STEPS = []


def step(lab, ok, det): STEPS.append((lab, bool(ok), det)); print(f"[{'PASS' if ok else 'FAIL'}] {lab}: {det}", flush=True)


def shims(R):
    d = os.path.join(R, "shims"); os.makedirs(d)
    real_make = shutil.which("make")
    mac = '-O2 -ffp-contract=off'
    for n, body in (
            ("module", "exit 0\n"),
            ("gcc", 'if [ "${1:-}" = --version ]; then echo "gcc (GCC) 14.3.0"; exit 0; fi\n'
                    'a=(); for x in "$@"; do case "$x" in -march=x86-64-v2|-mtune=generic) ;; *) a+=("$x");; esac; done\nexec cc "${a[@]}"\n'),
            ("make", f'cf="{mac} [mac rehearsal]"; for x in "$@"; do [ "$x" = koa-ld ] && cf="{mac} -DEDMD3_LONG_DOUBLE [mac rehearsal]"; done\n'
                     f'exec {real_make} "$@" \'CFLAGS_KOA={mac} $(INCLUDE)\' "BUILD_CF=$cf"\n'),
            ("ldd", 'exec otool -L "$@"\n')):
        p = os.path.join(d, n); open(p, "w").write("#!/bin/bash\n" + body); os.chmod(p, 0o755)
    return d


def run(cmd, env, cwd=None):
    r = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--commit", required=True)
    ap.add_argument("--standin"); a = ap.parse_args()
    R = os.path.abspath(a.root)
    if os.path.exists(R) and os.listdir(R): sys.exit(f"STOP: {R} exists and is not empty")
    home = os.path.join(R, "home"); os.makedirs(home); os.makedirs(os.path.join(R, "tmp"))
    S1 = os.path.join(home, "harddisks_gen3_s1"); head = K2R.clone(S1, F72, K2R.SPARSE_KOA); K2R.kissfft(S1); HS1 = os.path.join(S1, "hspist3")
    C = os.path.join(home, "harddisks_gen3_k2"); ch = K2R.clone(C, a.commit, K2R.SPARSE_K2); H2 = os.path.join(C, "hspist3")
    step("the fresh campaign clone at f72d764 and the scripts clone", head == F72 and ch.startswith(a.commit[:7]), f"{head[:12]}; {ch[:12]}")
    t = open(os.path.join(H2, P, "build_s1_koa.sh")).read()
    if t.count("/tmp/ld_$$") != 4: sys.exit("STOP: build_s1_koa.sh does not hold /tmp/ld_$$ four times")
    open(os.path.join(H2, P, "build_s1_koa.sh"), "w").write(t.replace("/tmp/ld_$$", os.path.join(R, "tmp", "ld_$$")))
    changed = [l.strip() for l in K2R.sh(["git", "-C", C, "status", "--porcelain", "--untracked-files=no"]).stdout.splitlines() if l.strip()]
    step("the scripts clone's shim (build_s1_koa.sh's /tmp)", changed == [f"M hspist3/{P}/build_s1_koa.sh"], "; ".join(changed))
    env = dict(os.environ, HOME=home, PATH=shims(R) + ":" + os.environ["PATH"], SLURM_JOB_ID="930001")
    rc, out = run(["bash", os.path.join(HS1, "cluster/gen3_koa_261009/build_gen3_koa.sh")], env)
    step("f72d764's own build_gen3_koa.sh in the new clone: refused", rc == 2 and "not harddisks_gen3_s1" in out and
         not os.path.exists(os.path.join(HS1, "00ALLINONE")), out.strip().splitlines()[-1])
    rc, out = run(["bash", os.path.join(H2, P, "build_s1_koa.sh")], env)
    open(os.path.join(R, "build_s1_koa_output.txt"), "w").write(out)
    last = out.strip().splitlines()
    step("block B: build_s1_koa.sh (the build)", rc == 0 and last[-1] == "BUILD OK: ./00ALLINONE (koa), bin/00ALLINONE_ld (koa-ld), bin/gen3_m1 and bin/gen3_m2 of f72d764",
         last[-1])
    print("== build_s1_koa.sh 2>&1 | tail -20\n" + "\n".join(last[-20:]))
    rc, out = run(["bash", os.path.join(H2, P, "record_s1_build.sh")], env)
    step("block B: record_s1_build.sh (the build-time record)", rc == 0 and out.strip().splitlines()[-1].startswith("BUILD RECORD OK: "), out.strip().splitlines()[-1])
    print("== record_s1_build.sh\n" + out)
    rc, out = run(["bash", os.path.join(H2, P, "build_s1_koa.sh")], env)
    step("STOP check: the build a second time", rc == 2 and "never rebuilt" in out and "BUILD OK" not in out, out.strip().splitlines()[-1][:160])
    rc, out = run(["bash", os.path.join(H2, P, "record_s1_build.sh")], env)
    step("STOP check: the record a second time", rc == 2 and "never rebuilt" in out, out.strip().splitlines()[-1][:160])
    rec = open(os.path.join(HS1, "logs", "BUILD_RECORD_S1.txt")).read()
    ok = all(f"sha256 {K2R.sha(os.path.join(HS1, b))} {b}" in rec for b in BINS)
    step("the record holds the four binaries' SHA-256", ok and f"commit {F72}" in rec, f"{sum(1 for l in rec.splitlines() if l.startswith('sha256 '))} sha256 lines")
    if a.standin:
        same = [b for b in BINS if K2R.sha(os.path.join(HS1, b)) == K2R.sha(os.path.join(a.standin, "hspist3", b))]
        print(f"INFORMATION: binaries byte-identical to the gate rehearsal's f72d764 stand-in: {len(same)} of 4 ({', '.join(same) or 'none'})")
    print("== the rehearsal\n\n| step | result | detail |\n|---|---|---|")
    for lab, ok, det in STEPS: print(f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |")
    print(f"\nBUILD REHEARSAL: {'PASS' if all(s[1] for s in STEPS) else 'FAIL'} ({sum(s[1] for s in STEPS)} of {len(STEPS)})")
    print(f"data under {R}: {K2R.sh(['du', '-sh', R]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
