#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE P2 RUNSHEET (KOA lines for Chris, ready to paste, no
placeholders but job ids), written from the package commit's hash given on the command line (a commit cannot name itself, so the
runsheet goes into the commit after the package). The blocks, in order (decision 57 item 3):
  A0 (login node)  squeue; the sandbox line
  A  (sandbox)     the fresh clone ~/harddisks_gen3_p2 at 8aaf5fef (fetched by hash), its sparse set, kissfft at febd4ca, the build with
                   8aaf5fef's cluster/build_koa.sh (make -B koa, cluster/koa_env.sh: stage 1's target, flags and environment); the fetch
                   and checkout of the package commit; record_p2_build.sh (the build-time record, the build inputs against
                   p2_build_sources.sha256); the package's SHA-256; the import check; the fresh scratch root
  B  (login node)  the guarded one-line submit of the array (decision 42's form: a second paste refuses)
  C  (login node)  squeue and the chunk checks (read only)
  D  (login node)  the merge job (it checks everything itself; never repeated)
  D2 (login node)  the merge's log
  E  (Mac)         the fetch
usage (from hspist3/): python3 cluster/gen3_stage2_261016/make_p2_runsheet.py <package commit, 40 hex> [--out <file>]
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); P = "cluster/gen3_stage2_261016"
BUILD = "8aaf5fef27c92391fc4fadd5133324740c29ef06"; KISS = "febd4caeed32e33ad8b2e0bb5ea77542c40f18ec"
ROOT = "/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/P2"
SPARSE = '"/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"'


def sha(rel): return hashlib.sha256(open(os.path.join(HS, rel), "rb").read()).hexdigest()


def main():
    if len(sys.argv) < 2 or len(sys.argv[1]) != 40: sys.exit(__doc__)
    PKG = sys.argv[1]; out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(HERE, "runsheet_p2_261016.txt")
    subj = subprocess.run(["git", "log", "-1", "--format=%s", PKG], cwd=HS, capture_output=True, text=True).stdout.strip()
    bsubj = subprocess.run(["git", "log", "-1", "--format=%s", BUILD], cwd=HS, capture_output=True, text=True).stdout.strip()
    if not subj: sys.exit(f"STOP: {PKG} is not a commit of this repository")
    chk = [f"{P}/p2_koa.sbatch", f"{P}/p2_merge.sbatch", f"{P}/p2_inputs.sha256", f"{P}/p2_build_sources.sha256", f"{P}/record_p2_build.sh"]
    nin = sum(1 for l in open(os.path.join(HS, P, "p2_inputs.sha256")) if l.strip())
    nbs = sum(1 for l in open(os.path.join(HS, P, "p2_build_sources.sha256")) if l.strip())
    tsha = open(os.path.join(HS, P, "tasks_p2_koa_261016.sha256")).read().strip()
    design = [l.rstrip("\n") for l in open(os.path.join(HS, P, "p2_design_261016.txt")) if l.startswith("- predicted") or l.startswith("- wall")]
    cls = sorted({int(l.split("\t")[9]) for l in list(open(os.path.join(HS, P, "chunks_p2_koa_261016.tsv")))[1:]})
    if len(cls) != 1: sys.exit(f"STOP: the chunk table has {len(cls)} time classes")
    tlim = cls[0]; subj80 = subj if len(subj) <= 80 else subj[:subj.rfind(" ", 0, 80)] + " ..."
    T = []
    T += [f"##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135) -- THE P2 RUNSHEET (KOA): the Stage 2 pilot, N 900, 400 runs.",
          "Written by cluster/gen3_stage2_261016/make_p2_runsheet.py. PREPARED: typed only after the plan author has read Report D57-1 and Chris",
          "has pushed main. Chris types every command.", "",
          "WHAT IT DOES. Block A makes ONE fresh clone ~/harddisks_gen3_p2: at 8aaf5fef it builds the campaign binary with 8aaf5fef's own",
          "cluster/build_koa.sh (make -B koa with cluster/koa_env.sh: stage 1's target, flags and environment); then it fetches the package commit",
          f"{PKG[:8]} (\"{subj80}\") and checks it out -- the package adds files only, the build inputs stay 8aaf5fef's, and",
          "record_p2_build.sh checks every one of them against p2_build_sources.sha256 before it writes the build-time record. Block B submits",
          "ONE array of 400 one-core tasks, at most 64 at once; each task runs one task line with p2_worker_koa.sh (the stage 1 worker plus the",
          "initial dump) after its start-up checks (inputs, task list, binary against the build-time record; no git). Block D submits the merge",
          "(the stage 1 merge and reducers, unchanged, then p2_reduce.py). Block E fetches the tables, logs, psi6 and virial files and the",
          "initial states to the Mac (no traces, no snapshots).",
          "KOA RULES, as before: the login node only for light commands; squeue before every sbatch; never submit a job twice; nothing in",
          "~/harddisks, ~/harddisks_resched, ~/harddisks_resched2, ~/harddisks_gen3, ~/harddisks_gen3_s1 or ~/harddisks_gen3_k2 changes; every job",
          "stops with \"STOP: [<node>] ...\" on a failed start-up check before it runs anything: paste the STOP line, change nothing, wait.",
          "NO GIT AT RUN TIME; every sbatch file carries --exclude=cn-05-08-07. THE SANDBOX LINE is the established one.",
          "",
          "CORE-HOURS AND TIME (Report D56-1 (a), the upper Z):"] + [f"  {d[2:]}" for d in design] + [
          f"  requested: 400 tasks x 0:{tlim:02d} h = {400 * tlim / 60:.0f} core-hours at most (Slurm charges use); the merge 1 core, --time=24:00:00,",
          "  expected well under 1 h.",
          "SCRATCH: about 0.9 GB on Lustre at most (per run, measured in the Mac rehearsal: the trace 0.5 MB until the merge gzips it to 0.07 MB,",
          "the snapshots 1.4 MB, psi6(t) 0.2 MB, the run folder with the initial state 0.12 MB); the fetch about 0.16 GB (0.4 MB per run).", "",
          "P2 BLOCK A0 (KOA login node; opens the sandbox shell; then wait for the cn-... prompt):",
          "squeue -u charing",
          "srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash", "",
          "  Expected: squeue lists no gen3 job (its header only). Paste block A only at the cn-... prompt.", "",
          "P2 BLOCK A (in the sandbox shell: the fresh clone, the build at 8aaf5fef, the package checkout, the build-time record, the checks):",
          "git --version",
          "ls -d ~/harddisks_gen3_p2 2>&1 | head -1",
          "cd ~ && git clone --depth 1 --branch main --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_gen3_p2",
          f"cd ~/harddisks_gen3_p2 && git fetch --depth 1 origin {BUILD}",
          f"git -c advice.detachedHead=false checkout {BUILD}",
          f"git sparse-checkout set --no-cone {SPARSE}",
          "git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft",
          f"git -C hspist3/kissfft -c advice.detachedHead=false checkout {KISS}",
          "git log -1 --format='%h %s' | cut -c1-80",
          "bash hspist3/cluster/build_koa.sh 2>&1 | tail -16",
          f"git fetch --depth 1 origin {PKG}",
          f"git -c advice.detachedHead=false checkout {PKG}",
          "git log -1 --format='%h %s' | cut -c1-80",
          "bash hspist3/cluster/gen3_stage2_261016/record_p2_build.sh",
          "cd hspist3 && sha256sum " + " ".join(chk),
          f"sha256sum -c {P}/p2_inputs.sha256 | grep -c ': OK$'",
          f"cat {P}/tasks_p2_koa_261016.sha256",
          "source cluster/koa_env.sh",
          "python3 -c \"import sys; sys.argv = ['x']; sys.path[:0] = ['cluster/gen3_stage2_261016', 'cluster/gen3_stage1_261012', 'cluster/gen3_testR_261011', "
          "'validation', '.']; import merge_s1_koa, reduce_s1_koa, p2_lattice, p2_reduce; print('the merge and reducer imports: OK')\"",
          f"ls {ROOT} 2>&1 | head -1",
          "exit", "",
          "  THE LINES THAT SHOW IT IS GOOD (paste everything):",
          "    \"git version 2.35.0\" or newer",
          "    \"ls: cannot access '/home/charing/harddisks_gen3_p2': No such file or directory\" (the clone does not exist yet)",
          f"    the clone; the fetch of {BUILD[:8]}...; \"HEAD is now at {BUILD[:7]} ...\"; the kissfft clone and \"HEAD is now at febd4ca ...\"",
          f"    \"{BUILD[:7]} {bsubj[:70]}\" (the log line)",
          "    the build's last lines: \"== checks\", the version lines with \"00ALLINONE  git 8aaf5fe  target koa\", \"== provenance\" with",
          "    \"compiler       gcc (GCC) 14.3.0\" and \"git_commit     8aaf5fef27c92391fc4fadd5133324740c29ef06\", the recorded hash line, \"BUILD OK\"",
          f"    the fetch of {PKG[:8]}...; \"HEAD is now at {PKG[:7]} ...\" (git may first print \"Warning: you are leaving 1 commit behind\": harmless)",
          f"    \"{PKG[:7]} {subj[:70]}\" (the log line)",
          "    the record: \"env: gcc (GCC) 14.3.0 | ...\", its lines (\"commit_build 8aaf5fef27c92391fc4fadd5133324740c29ef06\",",
          f"    \"commit_scripts {PKG}\", \"version 00ALLINONE  git 8aaf5fe  target koa\", \"sha256 <SHA-256> 00ALLINONE\",",
          f"    \"build_inputs {nbs} of {nbs} equal to 8aaf5fef's (...)\", the build_koa lines) and \"BUILD RECORD OK: /home/charing/harddisks_gen3_p2/hspist3/logs/BUILD_RECORD_P2.txt\"",
          "    exactly these five lines:"] + [f"      {sha(c)}  {c}" for c in chk] + [
          f"    {nin}         (the input manifest, each file \": OK\")",
          f"    {tsha}   (the task list's SHA-256)",
          "    \"env: gcc (GCC) 14.3.0 | Python ... at /home/charing/envs/hd/bin/python3 | ...\"",
          "    \"the merge and reducer imports: OK\"",
          f"    \"ls: cannot access '{ROOT}': No such file or directory\" (a fresh root)",
          "  STOP RULE: a STOP line or anything else: type exit if the shell is still open, paste it, submit nothing, wait.", "",
          "P2 BLOCK B (KOA login node, after block A's exit; ONE guarded line: it submits the array once, a second paste refuses):",
          "squeue -u charing",
          f"cd ~/harddisks_gen3_p2/hspist3 && {{ Q=$(squeue -h -u charing -n gen3-p2 -o %i 2>&1) && [ -z \"$Q\" ] || {{ echo \"REFUSED: squeue failed or lists a gen3-p2 job: $Q\"; false; }}; }} "
          f"&& {{ [ ! -e {ROOT} ] || {{ echo \"REFUSED: {ROOT} exists: the array ran before\"; false; }}; }} "
          f"&& {{ grep -qx 'commit_scripts {PKG}' logs/BUILD_RECORD_P2.txt || {{ echo \"REFUSED: logs/BUILD_RECORD_P2.txt does not name the package commit {PKG[:8]}\"; false; }}; }} "
          "&& echo \"CHECKS GOOD\" && sbatch --array=1-400%64 cluster/gen3_stage2_261016/p2_koa.sbatch && echo \"P2: THE ARRAY IS SUBMITTED\" || echo \"P2: THIS PASTE DID NOT SUBMIT THE ARRAY\"", "",
          "  Expected: squeue lists no gen3 job; then \"CHECKS GOOD\", \"Submitted batch job <its number>\" and \"P2: THE ARRAY IS SUBMITTED\". 400 tasks, at most 64 at once,",
          f"  0:{tlim:02d} h each at most; predicted about 50 min of wall time. A second paste prints \"REFUSED: ...\" and \"P2: THIS PASTE DID NOT SUBMIT THE ARRAY\".", "",
          "P2 BLOCK C (KOA login node, when squeue -u charing lists no gen3-p2 job; every line only reads; paste everything):",
          "squeue -u charing",
          "cd ~/harddisks_gen3_p2/hspist3",
          "ls logs/gen3-p2_*.out | wc -l",
          "grep -lx 'chunk [0-9]* done' logs/gen3-p2_*.out | wc -l",
          "cat logs/gen3-p2_*.out | grep -c -e FAILED -e REFUSED -e STOP",
          "grep -h '^chunk [0-9]* on ' logs/gen3-p2_*.out | sed 's/: group.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c",
          f"bash cluster/gen3_stage1_261012/s1_chunk_logs_koa.sh logs gen3-p2 {ROOT}/chunks {P}/chunks_p2_koa_261016.tsv P2", "",
          "  THE LINES THAT SHOW IT IS GOOD:",
          "    squeue: its header line only",
          "    400        (one log per chunk)",
          "    400        (each chunk once with its line \"chunk <j> done\")",
          "    0          (no FAILED, REFUSED or STOP line; a run the worker refused prints FAILED: paste it, it does not stop the merge)",
          "    one line per CPU model, the counts adding up to 400",
          "    \"chunk logs: 400; chunks done: 400 of 400 (group P2; ...); earlier attempts: 0 stopped in the start-up checks, 0 started and not finished\"",
          "  STOP RULE: anything else: paste it, change nothing, wait (a resubmission line follows from CC and the plan author).", "",
          "P2 BLOCK D (KOA login node, after block C was good; the merge job):",
          "squeue -u charing",
          "cd ~/harddisks_gen3_p2/hspist3 && sbatch cluster/gen3_stage2_261016/p2_merge.sbatch", "",
          "  Expected: squeue lists no gen3-p2 job; \"Submitted batch job <its number>\". Submit it ONCE: the merge checks everything itself before it",
          "  moves anything and refuses a second run (exit 3). It takes well under an hour.", "",
          "P2 BLOCK D2 (KOA login node, when squeue -u charing lists no gen3-p2-merge job; reads only; paste everything):",
          "squeue -u charing",
          "cd ~/harddisks_gen3_p2/hspist3 && L=$(ls -t logs/gen3-p2-merge_*.out | head -1) && echo \"$L\" && grep -Ev '^  m_|^reduce_' \"$L\" | tail -25", "",
          "  THE LINES THAT SHOW IT IS GOOD: squeue's header only; the log's path; \"merged 400 chunks: 1600 files moved, 800 files appended, 400 kept",
          "  per chunk (0 failed-run folders), 0 incomplete run folders, 0 locks or temp files left in place, 400 build files equal; build '00ALLINONE  git",
          "  8aaf5fe  target koa'\",",
          "  the lines of p2_reduce.py (40 folders, one line each: e start, e end, Q, stalls, psi6), \"p2_reduce: 40 folders (0 without the stage 1",
          "  tables), 400 trajectories, 0 with an error; <s> s\", \"P2 merge and reduction done; next: the fetch from the Mac (fetch_p2_koa.sh)\" and",
          "  \"this job's log is copied to ...\". STOP RULE: a STOP line: paste it, change nothing, wait.", "",
          "P2 BLOCK E (Mac, after block D2's log is good; one password and Duo push; about 0.16 GB):",
          "cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks && bash hspist3/cluster/gen3_stage2_261016/fetch_p2_koa.sh", "",
          "  Expected, last line: \"fetched into .../hspist3/experiments_gen3_stage2_koa_261016/P2: <size>; merge record: <time>; p2_reduce: p2_reduce: 40 folders",
          "  (0 without the stage 1 tables), 400 trajectories, 0 with an error; <s> s\". Then tell CC: CC runs p2_tables_mac.py and writes the P2 report."]
    open(out, "w").write("\n".join(T) + "\n"); print(f"written {out}: {len(T)} lines")


if __name__ == "__main__":
    main()
