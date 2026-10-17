#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 58; notes 261012 sec. 4.7.142): THE S2 RUNSHEET (KOA lines for Chris, ready to paste, NO placeholders),
written from the package commit's hash given on the command line (a commit cannot name itself, so the runsheet goes into the commit after the
package). Decision 58 item 2: THREE pastes, then read-only checks:
  PASTE A0+A (login node)  squeue; then ONE sandbox job, srun ... bash -s with the block as a here-document (no interactive shell to wait
                           for): the fresh clone ~/harddisks_gen3_s2 at the build commit 8308fb00 (the package's parent, fetched by hash), its
                           sparse set, kissfft at febd4ca, the build with its cluster/build_koa.sh (make -B koa, cluster/koa_env.sh: stage 1's
                           and P2's target, flags and environment); the fetch and checkout of the package commit; record_s2_build.sh; the
                           package's SHA-256; the input manifest; the import check; the fresh scratch root. Each critical step stops the job
                           at once on failure ("STOP: <step> failed"), so a failed clone never writes into the home folder.
  PASTE B (login node)     ONE guarded line: no gen3-s2 / gen3-s2-merge job queued, the scratch root absent, the build record naming the
                           package commit; then the array (sbatch --parsable: the job id J captured in the same line) and the merge held by
                           Slurm until every task of J has ended with exit 0 (--dependency=afterok:$J); a second paste refuses.
  CHECK C, CHECK D2        read only: the array's logs and chunk check once squeue lists no gen3-s2 task; the merge's log once it has ended.
  PASTE E (Mac)            the fetch.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/make_s2_runsheet.py <package commit, 40 hex> [--out <file>]
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); P = "cluster/gen3_stage2_261016"
BUILD = "8308fb004fdcb7e9c82197451cc04ab49ecd2f13"; KISS = "febd4caeed32e33ad8b2e0bb5ea77542c40f18ec"
ROOT = "/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/S2"
SPARSE = '"/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"'
NTASK = 2100


def sha(rel): return hashlib.sha256(open(os.path.join(HS, rel), "rb").read()).hexdigest()


def main():
    if len(sys.argv) < 2 or len(sys.argv[1]) != 40: sys.exit(__doc__)
    PKG = sys.argv[1]; out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(HERE, "runsheet_s2_261017.txt")
    git = lambda *a: subprocess.run(["git", *a], cwd=HS, capture_output=True, text=True).stdout.strip()
    subj, bsubj, parent = git("log", "-1", "--format=%s", PKG), git("log", "-1", "--format=%s", BUILD), git("rev-parse", f"{PKG}^")
    if not subj: sys.exit(f"STOP: {PKG} is not a commit of this repository")
    if parent != BUILD: sys.exit(f"STOP: the package's parent is {parent}, not the build commit {BUILD} (decision 58: build at the package's parent)")
    chk = [f"{P}/s2_koa.sbatch", f"{P}/s2_merge.sbatch", f"{P}/s2_inputs.sha256", f"{P}/s2_build_sources.sha256", f"{P}/record_s2_build.sh"]
    nin = sum(1 for l in open(os.path.join(HS, P, "s2_inputs.sha256")) if l.strip())
    nbs = sum(1 for l in open(os.path.join(HS, P, "s2_build_sources.sha256")) if l.strip())
    tsha = open(os.path.join(HS, P, "tasks_s2_koa_261017.sha256")).read().strip()
    design = [l.rstrip("\n") for l in open(os.path.join(HS, P, "s2_design_261017.txt")) if l.startswith("- predicted") or l.startswith("- wall")]
    cls = sorted({int(l.split("\t")[9]) for l in list(open(os.path.join(HS, P, "chunks_s2_koa_261017.tsv")))[1:]})
    if len(cls) != 1: sys.exit(f"STOP: the chunk table has {len(cls)} time classes")
    tlim = cls[0]; s80 = lambda s: s if len(s) <= 80 else s[:s.rfind(" ", 0, 80)] + " ..."
    B7, P7 = BUILD[:7], PKG[:7]
    T = [f"##CHRIS 2026-10-16 (plan-author decision 58; registration S2, notes 261012 sec. 4.9 and sec. 4.7.142) -- THE S2 RUNSHEET (KOA): Stage 2,",
         f"the hard-disk solid, c at zero strain; N 900 and 1600, {NTASK} runs. Written by cluster/gen3_stage2_261016/make_s2_runsheet.py. PREPARED: typed",
         "only after the plan author has read Report D58-1 and Chris has pushed main. Chris types every command.", "",
         f"WHAT IT DOES. Paste A0+A runs ONE sandbox job (srun, no interactive shell): it makes ONE fresh clone ~/harddisks_gen3_s2, builds the campaign",
         f"binary at {BUILD[:8]} (\"{s80(bsubj)}\", the package commit's parent) with that commit's own cluster/build_koa.sh (make -B koa with",
         f"cluster/koa_env.sh: stage 1's and P2's target, flags and environment; its 50 build inputs equal P2's 8aaf5fef's), then fetches the package",
         f"commit {PKG[:8]} (\"{s80(subj)}\") and checks it out -- the package adds files only -- and record_s2_build.sh checks every build input",
         "against s2_build_sources.sha256 before it writes the build-time record. Paste B submits, in ONE guarded line, the array of 2100 one-core",
         "tasks (at most 64 at once; each runs one task line with P2's worker p2_worker_koa.sh after its start-up checks; no git at run time) AND the",
         "merge, which Slurm holds until every task of the array has ended with exit 0 (--dependency=afterok); the merge is P2's (the stage 1 merge and",
         "reducers, unchanged, then p2_reduce.py) and checks everything itself before it moves anything. Paste E fetches the tables, logs, psi6 and",
         "virial files, the gzipped traces and the initial states to the Mac (no snapshots). CHECK C and CHECK D2 only read.",
         "KOA RULES, as before: the login node only for light commands; squeue before every sbatch; never submit a job twice; nothing in",
         "~/harddisks, ~/harddisks_resched, ~/harddisks_resched2, ~/harddisks_gen3, ~/harddisks_gen3_s1, ~/harddisks_gen3_k2 or ~/harddisks_gen3_p2",
         "changes; every job stops with \"STOP: ...\" on a failed check before it runs anything: paste the STOP line, change nothing, wait. NO GIT AT",
         "RUN TIME; every sbatch file carries --exclude=cn-05-08-07. The sandbox job is srun -p sandbox -t 1:00:00 -c 2 --mem=4G (the established line).", "",
         "CORE-HOURS AND TIME (Report D57-3 task 3, the upper Z):"] + [f"  {d[2:]}" for d in design] + [
         f"  requested: {NTASK} tasks x 0:{tlim:02d} h = {NTASK * tlim / 60:.0f} core-hours at most (Slurm charges use); the merge 1 core, --time=24:00:00,",
         "  expected about 4 h (P2's 47 min for 400 runs, scaled).",
         "SCRATCH: about 5 GB on Lustre at most (per run, P2's measured sizes: the trace 0.5 MB until the merge gzips it, the snapshots 1.4 MB at",
         "N 900 and about 2.5 MB at N 1600, psi6(t) 0.2 MB, the run folder 0.12 MB); the fetch about 1 GB (P2: 0.39 MB per run at N 900).", "",
         "S2 PASTE A0+A (KOA login node; ONE paste: the queue, then the sandbox job runs the whole block and ends by itself; paste everything it prints):",
         "squeue -u charing",
         "srun -p sandbox -t 1:00:00 -c 2 --mem=4G bash -s <<'S2A'",
         "exec 2>&1; export GIT_TERMINAL_PROMPT=0; set -v",
         "git --version",
         "ls -d ~/harddisks_gen3_s2 2>&1 | head -1",
         "cd ~ && git clone --depth 1 --branch main --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_gen3_s2 || { echo 'STOP: the clone failed'; exit 1; }",
         f"cd ~/harddisks_gen3_s2 && git fetch --depth 1 origin {BUILD} || {{ echo 'STOP: the fetch of the build commit failed'; exit 1; }}",
         f"git -c advice.detachedHead=false checkout {BUILD} || {{ echo 'STOP: the checkout of the build commit failed'; exit 1; }}",
         f"git sparse-checkout set --no-cone {SPARSE} || {{ echo 'STOP: the sparse checkout failed'; exit 1; }}",
         "git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft || { echo 'STOP: the kissfft clone failed'; exit 1; }",
         f"git -C hspist3/kissfft -c advice.detachedHead=false checkout {KISS} || {{ echo 'STOP: the kissfft checkout failed'; exit 1; }}",
         "git log -1 --format='%h %s' | cut -c1-80",
         "bash hspist3/cluster/build_koa.sh 2>&1 | tail -16",
         f"git fetch --depth 1 origin {PKG} || {{ echo 'STOP: the fetch of the package commit failed'; exit 1; }}",
         f"git -c advice.detachedHead=false checkout {PKG} || {{ echo 'STOP: the checkout of the package commit failed'; exit 1; }}",
         "git log -1 --format='%h %s' | cut -c1-80",
         "bash hspist3/cluster/gen3_stage2_261016/record_s2_build.sh || { echo 'STOP: the build record was not written'; exit 1; }",
         "cd hspist3 && sha256sum " + " ".join(chk),
         f"sha256sum -c {P}/s2_inputs.sha256 | grep -c ': OK$'",
         f"cat {P}/tasks_s2_koa_261017.sha256",
         "source cluster/koa_env.sh",
         "python3 -c \"import sys; sys.argv = ['x']; sys.path[:0] = ['cluster/gen3_stage2_261016', 'cluster/gen3_stage1_261012', 'cluster/gen3_testR_261011', "
         "'validation', '.']; import merge_s1_koa, reduce_s1_koa, p2_lattice, p2_reduce; print('the merge and reducer imports: OK')\"",
         f"ls {ROOT} 2>&1 | head -1",
         "echo 'S2 PASTE A0+A: END OF THE SANDBOX JOB'",
         "S2A", "",
         "  THE LINES THAT SHOW IT IS GOOD (every command appears once, as it runs, before its output):",
         "    squeue: its header only (no gen3 job); srun may print \"srun: job <n> queued and waiting for resources\" and \"... has been allocated resources\"",
         "    \"git version 2.35.0\" or newer",
         "    \"ls: cannot access '/home/charing/harddisks_gen3_s2': No such file or directory\" (the clone does not exist yet)",
         f"    the clone; the fetch of {BUILD[:8]}...; \"HEAD is now at {B7} ...\"; the kissfft clone and \"HEAD is now at febd4ca ...\"",
         f"    \"{B7} {bsubj[:70]}\" (the log line)",
         f"    the build's last lines: \"== checks\", the version lines with \"00ALLINONE  git {B7}  target koa\", \"== provenance\" with",
         f"    \"compiler       gcc (GCC) 14.3.0\" and \"git_commit     {BUILD}\", the recorded hash line, \"BUILD OK\"",
         f"    the fetch of {PKG[:8]}...; \"HEAD is now at {P7} ...\" (git may first print \"Warning: you are leaving 1 commit behind\": harmless)",
         f"    \"{P7} {subj[:70]}\" (the log line)",
         f"    the record: \"env: gcc (GCC) 14.3.0 | ...\", its lines (\"commit_build {BUILD}\",",
         f"    \"commit_scripts {PKG}\", \"version 00ALLINONE  git {B7}  target koa\", \"sha256 <SHA-256> 00ALLINONE\",",
         f"    \"build_inputs {nbs} of {nbs} equal to {BUILD[:8]}'s (...)\", the build_koa lines) and \"BUILD RECORD OK: /home/charing/harddisks_gen3_s2/hspist3/logs/BUILD_RECORD_S2.txt\"",
         "    exactly these five lines:"] + [f"      {sha(c)}  {c}" for c in chk] + [
         f"    {nin}         (the input manifest, each file \": OK\")",
         f"    {tsha}   (the task list's SHA-256)",
         "    \"env: gcc (GCC) 14.3.0 | Python ... at /home/charing/envs/hd/bin/python3 | ...\"",
         "    \"the merge and reducer imports: OK\"",
         f"    \"ls: cannot access '{ROOT}': No such file or directory\" (a fresh root)",
         "    \"S2 PASTE A0+A: END OF THE SANDBOX JOB\"",
         "  STOP RULE: a STOP line or anything else: paste it, submit nothing, wait.", "",
         "S2 PASTE B (KOA login node, after paste A0+A was good; ONE guarded line submits the array AND the dependent merge; a second paste refuses):",
         "squeue -u charing",
         f"cd ~/harddisks_gen3_s2/hspist3 && {{ Q=$(squeue -h -u charing -n gen3-s2,gen3-s2-merge -o %i 2>&1) && [ -z \"$Q\" ] || {{ echo \"REFUSED: squeue failed or lists a gen3-s2 job: $Q\"; false; }}; }} "
         f"&& {{ [ ! -e {ROOT} ] || {{ echo \"REFUSED: {ROOT} exists: the array ran before\"; false; }}; }} "
         f"&& {{ grep -qx 'commit_scripts {PKG}' logs/BUILD_RECORD_S2.txt || {{ echo \"REFUSED: logs/BUILD_RECORD_S2.txt does not name the package commit {PKG[:8]}\"; false; }}; }} "
         f"&& echo \"CHECKS GOOD\" && J=$(sbatch --parsable --array=1-{NTASK}%64 {P}/s2_koa.sbatch) && J=${{J%%;*}} && echo \"S2: THE ARRAY IS SUBMITTED: job $J\" "
         f"&& {{ sbatch --dependency=afterok:$J {P}/s2_merge.sbatch && echo \"S2: THE MERGE IS SUBMITTED, held by Slurm until every task of job $J has ended with exit 0\" "
         "|| echo \"S2: THE ARRAY $J IS SUBMITTED, THE MERGE IS NOT: paste this, submit nothing, wait\"; } || echo \"S2: THIS PASTE DID NOT SUBMIT THE ARRAY\"", "",
         f"  Expected: squeue lists no gen3 job; then \"CHECKS GOOD\", \"S2: THE ARRAY IS SUBMITTED: job <J>\", \"Submitted batch job <its number>\" and \"S2: THE MERGE",
         f"  IS SUBMITTED, held by Slurm until every task of job <J> has ended with exit 0\". {NTASK} tasks, at most 64 at once, 0:{tlim:02d} h each at most;",
         "  predicted about 5 h of wall time, then the merge (about 4 h). A second paste prints \"REFUSED: ...\" and \"S2: THIS PASTE DID NOT SUBMIT THE ARRAY\".",
         "  While the array runs, squeue lists its tasks and the merge as PENDING with reason (Dependency). If the merge shows",
         "  (DependencyNeverSatisfied), a task ended with a non-zero exit: run CHECK C, paste it, change nothing, wait.", "",
         "S2 CHECK C (KOA login node, read only, when squeue -u charing lists no gen3-s2 task; the merge may be running or done; paste everything):",
         "squeue -u charing",
         "cd ~/harddisks_gen3_s2/hspist3",
         "ls logs/gen3-s2_*.out | wc -l",
         "grep -lx 'chunk [0-9]* done' logs/gen3-s2_*.out | wc -l",
         "cat logs/gen3-s2_*.out | grep -c -e FAILED -e REFUSED -e STOP",
         "grep -h '^chunk [0-9]* on ' logs/gen3-s2_*.out | sed 's/: group.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c",
         f"bash cluster/gen3_stage1_261012/s1_chunk_logs_koa.sh logs gen3-s2 {ROOT}/chunks {P}/chunks_s2_koa_261017.tsv S2", "",
         "  THE LINES THAT SHOW IT IS GOOD:",
         "    squeue: the merge gen3-s2-merge RUNNING (or nothing, when it has ended); no gen3-s2 task",
         f"    {NTASK}       (one log per chunk)",
         f"    {NTASK}       (each chunk once with its line \"chunk <j> done\")",
         "    0          (no FAILED, REFUSED or STOP line; a run the worker refused prints FAILED: paste it, it does not stop the merge)",
         f"    one line per CPU model, the counts adding up to {NTASK}",
         f"    \"chunk logs: {NTASK}; chunks done: {NTASK} of {NTASK} (group S2; ...); earlier attempts: 0 stopped in the start-up checks, 0 started and not finished\"",
         "  STOP RULE: anything else: paste it, change nothing, wait (a resubmission line follows from CC and the plan author).", "",
         "S2 CHECK D2 (KOA login node, read only, when squeue -u charing lists no gen3-s2-merge job; paste everything):",
         "squeue -u charing",
         "cd ~/harddisks_gen3_s2/hspist3 && L=$(ls -t logs/gen3-s2-merge_*.out | head -1) && echo \"$L\" && grep -c '^stage2/gen3/.* errors 0;' \"$L\" && grep -Ev '^  m_|^reduce_|^stage2/gen3/' \"$L\" | tail -25", "",
         f"  THE LINES THAT SHOW IT IS GOOD: squeue's header only; the log's path; 210 (the (cell, H, mass) folders reduced without an error); \"merged {NTASK}",
         f"  chunks: {4 * NTASK} files moved, {2 * NTASK} files appended, {NTASK} kept per chunk (0 failed-run folders), 0 incomplete run folders, 0 locks or temp",
         f"  files left in place, {NTASK} build files equal; build '00ALLINONE  git {B7}  target koa'\"; \"p2_reduce: 210 folders (0 without the stage 1",
         f"  tables), {NTASK} trajectories, 0 with an error; <s> s\"; \"S2 merge and reduction done; next: the fetch from the Mac (fetch_s2_koa.sh)\" and \"this",
         "  job's log is copied to ...\". STOP RULE: a STOP line: paste it, change nothing, wait.", "",
         "S2 PASTE E (Mac, after CHECK D2 was good; one password and Duo push; about 1 GB):",
         "cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks && bash hspist3/cluster/gen3_stage2_261016/fetch_s2_koa.sh", "",
         "  Expected, last line: \"fetched into .../hspist3/experiments_gen3_stage2_koa_261016/S2: <size>; merge record: <time>; p2_reduce: p2_reduce: 210 folders",
         f"  (0 without the stage 1 tables), {NTASK} trajectories, 0 with an error; <s> s\". Then tell CC: CC runs s2_tables_mac.py and writes the S2 report."]
    open(out, "w").write("\n".join(T) + "\n"); print(f"written {out}: {len(T)} lines")


if __name__ == "__main__":
    main()
