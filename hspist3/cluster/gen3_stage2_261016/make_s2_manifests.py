#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 58; notes 261012 sec. 4.7.142): THE TWO MANIFESTS OF THE S2 PACKAGE (sha256sum -c format, paths
relative to hspist3/), written on the Mac before the package commit. It is make_p2_manifests.py with the
build commit 8308fb00 (the package commit's parent) and the S2 file names:
  s2_build_sources.sha256  the build inputs of 8308fb00, from git itself (git show 8308fb00:<path>, not the working tree): hspist3/Makefile,
                           every tracked *.c and *.h of hspist3/ and of hspist3/edmd_core/, cluster/build_koa.sh and cluster/koa_env.sh;
                           and kissfft's top-level *.c and *.h at febd4ca (git -C hspist3/kissfft show febd4ca:<file>; the Mac clone's
                           local change to kiss_fft_log.h is not in it). record_s2_build.sh checks the KOA clone's files against it.
  s2_inputs.sha256         every file an S2 job reads: the task list, its SHA-256 file, the chunk table, the worker, the build-record
                           script and the build-source manifest, the stage 1 chunk check, cluster/koa_env.sh, and the repository's Python
                           modules that the merge's scripts load (merge_s1_koa.py, reduce_B.py, reduce_s1_koa.py, p2_reduce.py: their import
                           closure inside hspist3/, found by importing them in a fresh interpreter). s2_koa.sbatch and s2_merge.sbatch check it (sha256sum -c).
usage (from hspist3/): python3 cluster/gen3_stage2_261016/make_s2_manifests.py
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL); REPO = os.path.dirname(HS)
BUILD, KISS = "8308fb004fdcb7e9c82197451cc04ab49ecd2f13", "febd4caeed32e33ad8b2e0bb5ea77542c40f18ec"
P = "cluster/gen3_stage2_261016"


def git(*a, cwd=REPO): return subprocess.run(["git", *a], cwd=cwd, capture_output=True, check=True).stdout


def build_sources():
    names = git("ls-tree", "-r", "--name-only", BUILD, "hspist3/").decode().split("\n"); out = []
    for n in sorted(x for x in names if x):
        r = n[len("hspist3/"):]
        top = "/" not in r and r.endswith((".c", ".h"))
        core = r.startswith("edmd_core/") and r.count("/") == 1 and r.endswith((".c", ".h"))
        if r == "Makefile" or top or core or r in ("cluster/build_koa.sh", "cluster/koa_env.sh"):
            out.append((hashlib.sha256(git("show", f"{BUILD}:{n}")).hexdigest(), r))
    kd = os.path.join(HS, "kissfft")
    for n in sorted(git("ls-tree", "--name-only", KISS, cwd=kd).decode().split()):
        if n.endswith((".c", ".h")): out.append((hashlib.sha256(git("show", f"{KISS}:{n}", cwd=kd)).hexdigest(), f"kissfft/{n}"))
    return out


def closure(scripts):
    """the repository's Python files that the merge job loads: a fresh interpreter imports p2_reduce, merge_s1_koa and reduce_s1_koa (their
    main() is guarded) and reduce_B.py's two imports (tests_20260913, paper1_populate_cs_err_20261002; reduce_B runs at import), with
    the search path the scripts give themselves, and lists every loaded module file inside hspist3/"""
    code = ("import sys, os; HS = sys.argv[1]; sys.argv = ['x']; "
            "sys.path[:0] = [os.path.join(HS, 'cluster', 'gen3_stage2_261016'), os.path.join(HS, 'cluster', 'gen3_stage1_261012'), "
            "os.path.join(HS, 'cluster', 'gen3_testR_261011'), os.path.join(HS, 'validation'), HS]; "
            "import p2_reduce, merge_s1_koa, reduce_s1_koa, tests_20260913, paper1_populate_cs_err_20261002; "
            "print('\\n'.join(sorted({os.path.relpath(m.__file__, HS) for m in list(sys.modules.values()) "
            "if getattr(m, '__file__', None) and os.path.abspath(m.__file__).startswith(HS + os.sep)})))")
    out = subprocess.run([sys.executable, "-B", "-c", code, HS], capture_output=True, text=True, check=True).stdout.split()
    return sorted(set(out) | set(scripts))


def write(path, rows):
    open(path, "w").write("".join(f"{h}  {r}\n" for h, r in rows)); return len(rows)


def main():
    bs = build_sources(); nb = write(os.path.join(HERE, "s2_build_sources.sha256"), bs)
    py = closure([f"{P}/p2_reduce.py", "cluster/gen3_stage1_261012/merge_s1_koa.py", "cluster/gen3_stage1_261012/reduce_s1_koa.py",
                  "cluster/confinement_20261013/reduce_B.py"])
    other = [f"{P}/tasks_s2_koa_261017.txt", f"{P}/tasks_s2_koa_261017.sha256", f"{P}/chunks_s2_koa_261017.tsv", f"{P}/p2_worker_koa.sh",
             f"{P}/record_s2_build.sh", f"{P}/s2_build_sources.sha256", "cluster/gen3_stage1_261012/s1_chunk_logs_koa.sh", "cluster/koa_env.sh"]
    rows = [(hashlib.sha256(open(os.path.join(HS, f), "rb").read()).hexdigest(), f) for f in sorted(set(other) | set(py))]
    ni = write(os.path.join(HERE, "s2_inputs.sha256"), rows)
    print(f"s2_build_sources.sha256: {nb} files (8308fb00's tracked build inputs and kissfft at febd4ca)")
    print(f"s2_inputs.sha256: {ni} files, of them {len(py)} Python files of the merge's import closure:")
    for f in py: print(f"  {f}")


if __name__ == "__main__":
    main()
