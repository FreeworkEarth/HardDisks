#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 38 item 6; 261012 sec. 4.7.80): THE MAC REHEARSAL of the M5 solid-cell job, end to end, before
Chris types it (decision 20's rule), as rehearse_k2_mac.py rehearsed K2-2. Nothing outside --root is written.
  home/harddisks_gen3     a symlink to --gen3-home: a stand-in for KOA's ~/harddisks_gen3, i.e. a sparse clone of a643bc7 with the Mac
                          builds of the Makefile's koa and koa-ld targets and the build-time hash files, as rehearse_k2_mac.py made it
                          (checked first: both version lines and both hash files; the job only reads it).
  home/harddisks_gen3_k2  a fresh shared clone of --commit (this package's commit), checked out with the K2 clone's sparse set
                          (/hspist3/*.py, /hspist3/cluster/, /hspist3/validation/): runs the job, the fetch and the verdict.
SHIMS (the residual differences from KOA, printed at the end): the Slurm variables (SLURM_SUBMIT_DIR, SLURM_JOB_ID, SLURM_JOB_NAME) and
the log file Slurm would write (logs/%x_%j.out); `module` (a no-op), `gcc --version` (prints 14.3.0; nothing is compiled), `lscpu` (the
Mac's CPU); HOME; four files of the clone, rewritten in place: k2_binaries_a643bc7.sha256 (the Mac builds' SHA-256 instead of KOA's),
m5_solid_inputs.sha256 (that file's new SHA-256), the job's Lustre root (<root>/lustre) and the fetch's DTN source (made local). On arm64
the long double is the 8-byte double, so the koa-ld runs exercise the plumbing only, not the 80-bit numerics (as in K2's rehearsal).
ORDER: the checks of --gen3-home; the clone and the shims; the job (8 simulation processes at most); two STOP checks (the same job id
again; a submission from the wrong clone); the fetch; the verdict.
usage (from hspist3/): python3 cluster/gen3_m5_solid_261012/rehearse_m5_solid_mac.py --root R --commit C --gen3-home H
"""
import argparse, hashlib, os, re, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
P = "cluster/gen3_m5_solid_261012"; KFILE = "cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256"; JOB = "9100001"; NAME = "gen3-m5-solid"


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def sh(cmd, **kw): return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def rewrite(path, old, new):
    t = open(path).read()
    if t.count(old) != 1: sys.exit(f"STOP: {path} does not hold exactly one '{old}'")
    open(path, "w").write(t.replace(old, new))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--commit", required=True)
    ap.add_argument("--gen3-home", required=True); a = ap.parse_args()
    R = os.path.abspath(a.root); G = os.path.abspath(a.gen3_home); B = os.path.join(G, "hspist3")
    if os.path.exists(R) and os.listdir(R): sys.exit(f"STOP: {R} exists and is not empty")
    steps = []
    # ------------------------------------------------------------ the stand-in for ~/harddisks_gen3
    v = sh([os.path.join(B, "00ALLINONE"), "--version"]).split("\n")[0]; vl = sh([os.path.join(B, "bin", "00ALLINONE_ld"), "--version"]).split("\n")[0]
    h = sha(os.path.join(B, "00ALLINONE")); hl = sha(os.path.join(B, "bin", "00ALLINONE_ld"))
    ok = (v == "00ALLINONE  git a643bc7  target koa" and vl == "00ALLINONE  git a643bc7  target koa-ld"
          and open(os.path.join(B, "logs", "BUILD_KOA_LAST.hash")).read().split() == [h, "00ALLINONE"]
          and open(os.path.join(B, "logs", "BUILD_KOA_LD.hash")).read().split() == [hl, "00ALLINONE_ld"])
    steps.append(("the stand-in for ~/harddisks_gen3 (version lines, build-time hash files)", ok, f"{v} {h[:12]}...; {vl} {hl[:12]}..."))
    if not ok: sys.exit(f"STOP: {G} is not a643bc7 with its koa / koa-ld Mac builds")
    # ------------------------------------------------------------ the clone and the shims
    home = os.path.join(R, "home"); os.makedirs(home); os.symlink(G, os.path.join(home, "harddisks_gen3"))
    C = os.path.join(home, "harddisks_gen3_k2"); repo = os.path.dirname(sh(["git", "-C", HS, "rev-parse", "--path-format=absolute", "--git-common-dir"]).strip())
    sh(["git", "clone", "-q", "--no-checkout", "--shared", repo, C]); sh(["git", "-C", C, "config", "core.abbrev", "7"])
    sh(["git", "-C", C, "config", "core.sparseCheckout", "true"])
    with open(os.path.join(C, ".git", "info", "sparse-checkout"), "w") as fh: fh.write("/hspist3/*.py\n/hspist3/cluster/\n/hspist3/validation/\n")
    sh(["git", "-C", C, "-c", "advice.detachedHead=false", "checkout", "-q", a.commit])
    H2 = os.path.join(C, "hspist3"); os.makedirs(os.path.join(H2, "logs"), exist_ok=True)
    head = sh(["git", "-C", C, "log", "-1", "--format=%h %s"]).strip()[:100]
    with open(os.path.join(H2, KFILE), "w") as fh:
        for f in ("00ALLINONE", "bin/00ALLINONE_ld", "bin/gen3_m1", "bin/gen3_m2"): fh.write(f"{sha(os.path.join(B, f))}  {f}\n")
    man = os.path.join(H2, P, "m5_solid_inputs.sha256"); t = open(man).read()
    t2 = re.sub(r"^[0-9a-f]{64}  " + re.escape(KFILE) + "$", f"{sha(os.path.join(H2, KFILE))}  {KFILE}", t, flags=re.M)
    if t2 == t: sys.exit("STOP: the manifest has no line for the binaries' hash file")
    open(man, "w").write(t2)
    rewrite(os.path.join(H2, P, "m5_solid_koa.sbatch"), "/mnt/lustre/koa/scratch/charing", os.path.join(R, "lustre"))
    rewrite(os.path.join(H2, P, "fetch_m5_solid_koa.sh"), 'SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/m5_solid_koa"',
            f'SRC="{os.path.join(R, "lustre", "harddisks_gen3", "m5_solid_koa")}"')
    changed = sh(["git", "-C", C, "status", "--porcelain", "--untracked-files=no"]).split("\n")
    steps.append((f"the clone at {head}; the files the shims changed", len([x for x in changed if x.strip()]) == 4, "; ".join(x.strip() for x in changed if x.strip())))
    shims = os.path.join(R, "shims"); os.makedirs(shims)
    cpu = sh(["sysctl", "-n", "machdep.cpu.brand_string"]).strip()
    for name, body in (("module", "exit 0\n"), ("gcc", 'echo "gcc (GCC) 14.3.0"\n'), ("lscpu", f'echo "Model name:            {cpu}"\n')):
        p = os.path.join(shims, name); open(p, "w").write("#!/bin/bash\n" + body); os.chmod(p, 0o755)
    env = dict(os.environ, HOME=home, PATH=shims + ":" + os.environ["PATH"], SLURM_SUBMIT_DIR=H2, SLURM_JOB_ID=JOB, SLURM_JOB_NAME=NAME)
    # ------------------------------------------------------------ the job
    logp = os.path.join(H2, "logs", f"{NAME}_{JOB}.out")
    with open(logp, "w") as fh: rc = subprocess.run(["bash", f"{P}/m5_solid_koa.sbatch"], cwd=H2, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    log = open(logp).read()
    steps.append(("the job (m5_solid_koa.sbatch)", rc == 0 and "M5 solid cell done" in log and "STOP" not in log, f"exit {rc}"))
    print("== the job's log\n" + log)
    # ------------------------------------------------------------ two STOP checks
    os.makedirs(os.path.join(R, "stopchecks"))
    for lab, e2, want in (("the same job id again", env, "exists -- a job never writes into another job's folder"),
                          ("a submission from ~/harddisks_gen3", dict(env, SLURM_SUBMIT_DIR=os.path.join(home, "harddisks_gen3", "hspist3")),
                           "the M5 job runs from ~/harddisks_gen3_k2, not harddisks_gen3")):
        o = os.path.join(R, "stopchecks", re.sub(r"\W+", "_", lab) + ".out")
        with open(o, "w") as fh: rc2 = subprocess.run(["bash", os.path.join(H2, P, "m5_solid_koa.sbatch")], cwd=H2, env=e2, stdout=fh, stderr=subprocess.STDOUT).returncode
        t = open(o).read().strip()
        steps.append((f"STOP check: {lab}", rc2 == 2 and want in t and "M5 solid cell on" not in t, t.split("\n")[-1][:160]))
    # ------------------------------------------------------------ the fetch and the verdict
    f = subprocess.run(["bash", os.path.join(H2, P, "fetch_m5_solid_koa.sh")], cwd=C, env=env, capture_output=True, text=True)
    steps.append(("the fetch (fetch_m5_solid_koa.sh, source made local)", f.returncode == 0 and "job_" + JOB in f.stdout, f.stdout.strip().split("\n")[-1][:160]))
    vd = subprocess.run(["python3", os.path.join(P, "m5_solid_verdict_261012.py")], cwd=H2, env=env, capture_output=True, text=True)
    steps.append(("the verdict script on the fetched folder", vd.returncode == 0 and "M5 AT THE SOLID CELL, BY THE REGISTERED RULE: the double run's" in vd.stdout,
                  (vd.stdout.strip().split("\n") or [""])[-1][:160]))
    print("== the verdict script's output (the Mac's koa-ld is the double: plumbing only)\n" + vd.stdout + vd.stderr)
    # ------------------------------------------------------------ summary
    print("== the rehearsal\n\n| step | result | detail |\n|---|---|---|")
    for lab, ok, det in steps: print(f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |")
    print(f"\nREHEARSAL: {'PASS' if all(s[1] for s in steps) else 'FAIL'} ({sum(s[1] for s in steps)} of {len(steps)})")
    print("\nResidual differences from KOA (the shims): the Slurm variables and log file; module, gcc --version, lscpu; HOME; the four rewritten "
          "files above; arm64's long double is the double (the koa-ld numerics are not exercised).")
    print(f"data under {R}: {sh(['du', '-sh', R]).split()[0]}")


if __name__ == "__main__":
    main()
