#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.8 and 4.8.10; gate G1): WRITES THE STAGE 1 RUNSHEET
(Part B) from the chunk table and the task list, so that every array line, count and expected line is computed, not typed. The blocks:
  A, B   the sandbox: the update in place of ~/harddisks_gen3_k2 to the stage 1 commit (the commit after the gate that adds
         build_record_f72d764.txt, a copy of the gate's fetched build-time record, and the two input manifests), its checks, the unit test;
  C, D   chunk 1 alone and its checks (gate G6 (e));
  per group, in the registered order (sec. 4.8.8: B100, B400, D, B900, S, B1600):
  E-<g>  the group's arrays, one per time class (one array has one time limit), throttled so that at most 64 run at once over all of
         them (the slots given greedily by predicted core-hours per slot, at most the class's chunk count); the next group only after
         the previous group's merge;
  F-<g>  the checks of the chunk logs and the resubmission line(s) of chunks that stopped in their start-up checks (Test R's rule);
  G-<g>  the merge;   H-<g>  its record;   I-<g>  the fetch (Mac).
  FG-<g> ##CHRIS 2026-10-13 (plan-author decision 45 addendum 2): for the remaining groups (D, B900, S, B1600) blocks F and G are ONE
         paste: F's lines, then "if <F's five good conditions>; then <block G's guarded line>; else echo 'F-<g> NOT CLEAN -- the merge
         was not submitted'; fi". The five: logs = chunks, done logs = chunks, 0 FAILED or REFUSED, the chunk check exit 0, "nothing to
         resubmit". B100's and B400's blocks F and G stay as they were typed.
The rehearsal (rehearse_s1_mac.py, gate G4) runs this script on its mini chunk table and runs the blocks it writes.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/make_s1_runsheet.py [--commit <stage 1 commit>] [--chunks T] [--tasks L] [--out F]
  without --commit the update lines name AFTER_THE_GATE (the commit exists only after the gate's fetch) and the SHA-256 lines of the build
  record and the manifests are written after the gate
"""
import argparse, collections, hashlib, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); P = os.path.relpath(HERE, HS)
GROUPS = ["B100", "B400", "D", "B900", "S", "B1600"]
LUS = "/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1"
CAP = 64
BUILD = "00ALLINONE  git f72d764  target koa"
PKG = ["s1_koa.sbatch", "s1_worker_koa.sh", "s1_merge.sbatch", "merge_s1_koa.py", "reduce_s1_koa.py", "s1_chunk_logs_koa.sh",
       "s1_stopped_chunks.sh", "test_merge_s1_plan.py", "tasks_stage1_koa_261012.sha256", "chunks_stage1_koa_261012.tsv"]
AFTER = ["build_record_f72d764.txt", "s1_inputs.sha256", "s1_merge_inputs.sha256"]
FG_ONE = ("D", "B900", "S", "B1600")                  # ##CHRIS 2026-10-13 (decision 45 addendum 2): the groups whose F and G are one paste


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def hm(minutes): return f"{minutes // 60}:{minutes % 60:02d}:00"


def ranges(js):
    js = sorted(js); out = []; a = b = js[0]
    for j in js[1:]:
        if j == b + 1: b = j
        else: out.append(f"{a}-{b}" if b > a else f"{a}"); a = b = j
    out.append(f"{a}-{b}" if b > a else f"{a}")
    return ",".join(out)


def check(test, reason):
    """##CHRIS 2026-10-12 (plan-author decision 42 item 5): one guard of a submit line -- the test, or "REFUSED: <reason>" and false"""
    return "{ " + test + " || { echo \"REFUSED: " + reason + "\"; false; }; }"


def guarded(checks, cmds, ok, fail):
    """decision 42 item 5: one command line -- every check, then the commands unchanged and in their order; prints ok, or fail (a check
    refused, or a command failed: then the lines "Submitted batch job" show how many were submitted)"""
    return ("cd ~/harddisks_gen3_k2/hspist3 && " + " && ".join(checks) + " && echo \"CHECKS GOOD\" && " + " && ".join(cmds)
            + f" && echo \"{ok}\" || echo \"{fail}\"")


def sq(names):
    return check(f"Q=$(squeue -h -u charing -n {names} -o %i 2>&1) && [ -z \"$Q\" ]", f"squeue failed or lists a {names.replace(',', ' or ')} job: $Q")


def head_ok(commit):
    return check(f"[ \"$(cat ../.git/HEAD)\" = {commit} ]", f"the scripts clone is not at {commit[:7]}")


def awk_str(s):
    if "'" in s or "\\" in s: sys.exit("STOP: a guard text holds a single quote or a backslash")
    return '"' + s.replace('"', '\\"') + '"'


def f_awk(g, commit, P):
    """block F's resubmission line, guarded (decision 42 item 5): s1_stopped_chunks.sh (unchanged, 7bdfb67e) prints its plain lines; this
    awk prints them as they are and then ONE guarded line that holds all of them (a second line of its own would find the first one's
    job in squeue and refuse): squeue empty for the group, HEAD, and no log of a listed chunk with its done line"""
    PFX = f"gen3-s1-{g}"
    pre = ("cd ~/harddisks_gen3_k2/hspist3 && " + sq(PFX) + " && " + head_ok(commit) + " && { ! grep -Elx \"chunk (")
    mid = (f") done\" logs/{PFX}_*.out >/dev/null 2>&1 || {{ echo \"REFUSED: a listed chunk has a log with its done line\"; false; }}; }} && echo \"CHECKS GOOD\" && ")
    prog = ("{ print } /^sbatch --job-name=" + PFX + " / { n++; L[n] = $0; a = $0; sub(/.*--array=/, \"\", a); sub(/ .*/, \"\", a); "
            "k = split(a, p, \",\"); for (i = 1; i <= k; i++) { m = split(p[i], q, \"-\"); lo = q[1] + 0; hi = (m == 2) ? q[2] + 0 : lo; "
            "for (j = lo; j <= hi; j++) J = J (J == \"\" ? \"\" : \"|\") j } } "
            "END { if (n) { s = \"\"; for (i = 1; i <= n; i++) s = s L[i] \" && \"; "
            "print \"THE GUARDED RESUBMISSION LINE (one line for all of them; type it once, after the go of the plan author):\"; "
            "print " + awk_str(pre) + " J " + awk_str(mid) + " s " + awk_str(f"echo \"{g}: ALL ") + " n "
            + awk_str(f" RESUBMISSION LINE(S) SUBMITTED\" || echo \"{g}: THIS PASTE DID NOT SUBMIT ALL ") + " n " + awk_str(" RESUBMISSION LINE(S)\"") + " } }")
    if "'" in prog: sys.exit("STOP: the awk program holds a single quote")
    return f"bash {P}/s1_stopped_chunks.sh logs {PFX} {P}/chunks_stage1_koa_261012.tsv {g} | awk '{prog}'"


def g_line(g, R, PFX, commit):
    """block G's guarded line (decision 42 item 5), as before; ##CHRIS 2026-10-13: a function, so block FG holds the same line"""
    return guarded([sq(f"{PFX},gen3-s1-merge-{g}"), head_ok(commit), check(f"[ ! -e {R}/merged ]", f"{R}/merged exists")],
                   [f"sbatch --job-name=gen3-s1-merge-{g} {P}/s1_merge.sbatch {g}"], f"{g}: THE MERGE SUBMITTED", f"{g}: THIS PASTE DID NOT SUBMIT THE MERGE")


def fg_line(g, nch, R, PFX, commit):
    """##CHRIS 2026-10-13 (plan-author decision 45 addendum 2): the last line of block FG -- block G's guarded line, unchanged, only if F's
    five good conditions hold: logs = chunks, done logs = chunks, 0 FAILED or REFUSED, the chunk check exit 0, 'nothing to resubmit'; else
    "F-<g> NOT CLEAN -- the merge was not submitted". The five are computed again here (F's lines above only print them)."""
    ct = f"{P}/chunks_stage1_koa_261012.tsv"
    five = [f"[ \"$(ls logs/{PFX}_*.out 2>/dev/null | wc -l)\" -eq {nch} ]",
            f"[ \"$(grep -lx 'chunk [0-9]* done' logs/{PFX}_*.out 2>/dev/null | wc -l)\" -eq {nch} ]",
            f"[ \"$(cat logs/{PFX}_*.out 2>/dev/null | grep -c -e FAILED -e REFUSED)\" -eq 0 ]",
            f"bash {P}/s1_chunk_logs_koa.sh logs {PFX} {R}/chunks {ct} {g} >/dev/null 2>&1",
            f"bash {P}/s1_stopped_chunks.sh logs {PFX} {ct} {g} 2>/dev/null | grep -qx 'nothing to resubmit'"]
    return ("if cd ~/harddisks_gen3_k2/hspist3 && " + " && ".join(five) + "; then " + g_line(g, R, PFX, commit)
            + f"; else echo 'F-{g} NOT CLEAN -- the merge was not submitted'; fi")


def throttles(cls):
    """cls: {class: (number of chunks, predicted core-h)} -> {class: throttle}; each >= 1 and <= its chunk count, the sum <= CAP;
    the slots left are given one by one to the class with the largest predicted core-h per slot (greedy: this minimises the predicted
    wall time of the slowest class)"""
    k = {c: 1 for c in cls}; left = CAP - len(cls)
    if left < 0: sys.exit(f"STOP: {len(cls)} classes in one group exceed {CAP} cores at one task each")
    while left > 0:
        cand = [c for c in cls if k[c] < cls[c][0]]
        if not cand: break
        c = max(cand, key=lambda c: (cls[c][1] / k[c], -c)); k[c] += 1; left -= 1
    return k


def load(chunks, tasks):
    CT = [l.rstrip("\n").split("\t") for l in open(chunks)][1:]
    TL = [l.split() for l in open(tasks) if l.strip()]
    G = collections.OrderedDict()
    for g in GROUPS:
        rows = [c for c in CT if c[4] == g]
        if not rows: continue
        lines = [TL[i - 1] for c in rows for i in range(int(c[1]), int(c[2]) + 1)]
        G[g] = dict(rows=rows, nch=len(rows), traj=len(lines), snap=sum(1 for l in lines if float(l[12]) != 0.0),
                    pairs=len({(l[1]) for l in lines}), cells=len({l[1].split("/")[2] for l in lines}),
                    core_h=sum(float(c[7]) for c in rows) / 3600.0, req_h=sum(int(c[9]) for c in rows) / 60.0,
                    first=int(rows[0][0]), last=int(rows[-1][0]))
    return CT, TL, G


def arrays(g, info, skip=()):
    by = collections.defaultdict(list)
    for c in info["rows"]:
        if int(c[0]) not in skip: by[int(c[9])].append(c)
    k = throttles({cl: (len(rs), sum(float(c[7]) for c in rs) / 3600.0) for cl, rs in by.items()})
    out = []
    for cl in sorted(by):
        js = [int(c[0]) for c in by[cl]]
        out.append((cl, len(js), k[cl], sum(float(c[7]) for c in by[cl]) / 3600.0,
                    f"sbatch --job-name=gen3-s1-{g} --time={hm(cl)} --array={ranges(js)}%{k[cl]} {P}/s1_koa.sbatch"))
    return out


def block_b(commit, here):
    files = [f"{P}/{f}" for f in PKG]; after = [f"{P}/{f}" for f in AFTER]
    have_after = all(os.path.exists(os.path.join(HS, f)) for f in after)
    shown = files + (after if have_after else [])
    nman = (sum(1 for l in open(os.path.join(HS, P, "s1_inputs.sha256")) if l.strip())
            + sum(1 for l in open(os.path.join(HS, P, "s1_merge_inputs.sha256")) if l.strip())) if have_after else None
    tsha = open(os.path.join(HS, P, "tasks_stage1_koa_261012.sha256")).read().strip()
    L = ["git --version", "cd ~/harddisks_gen3_k2", "git status --porcelain --untracked-files=no | wc -l",
         f"git fetch --depth 1 origin {commit}", f"git -c advice.detachedHead=false checkout {commit}",
         "git log -1 --format='%h %s' | cut -c1-80",
         "cd hspist3 && sha256sum " + " ".join(shown),
         f"sha256sum -c {P}/s1_inputs.sha256 {P}/s1_merge_inputs.sha256 | grep -c ': OK$'",
         f"cat {P}/tasks_stage1_koa_261012.sha256",
         f"cmp {P}/build_record_f72d764.txt ~/harddisks_gen3_s1/hspist3/logs/BUILD_RECORD_S1.txt && echo \"the committed build record equals the build-time record\"",
         "(cd ~/harddisks_gen3_s1/hspist3 && sha256sum 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2)",
         "~/harddisks_gen3_s1/hspist3/00ALLINONE --version | head -1",
         "source cluster/koa_env.sh",
         f"mkdir -p logs/s1_unit_tmp && TMPDIR=$PWD/logs/s1_unit_tmp python3 {P}/test_merge_s1_plan.py > logs/s1_unit_merge.txt 2>&1; tail -1 logs/s1_unit_merge.txt",
         f"python3 -c \"import sys; sys.argv = ['x']; sys.path[:0] = ['{P}', 'cluster/gen3_testR_261011', 'validation', '.']; import merge_s1_koa, reduce_s1_koa; print('the merge and reducer imports: OK')\"",
         f"ls {LUS} 2>&1 | head -2", "exit"]
    good = ["\"git version 2.35.0\" or newer", "0          (no tracked file of the scripts clone has changed)",
            (f"the fetch of {commit} and \"HEAD is now at {commit[:7]} ...\" and the log line starting with {commit[:7]}" if commit != "AFTER_THE_GATE" else
             "the fetch of the stage 1 commit and \"HEAD is now at <its first 7 characters> ...\" and its log line")
            + " (git may first print \"Warning: you are leaving <n> commit(s) behind ...\": harmless)"]
    good.append(f"exactly these {len(shown)} lines:")
    good += [f"  {sha(os.path.join(HS, f))}  {f}" for f in shown]
    if not have_after: good.append("  and the SHA-256 lines of build_record_f72d764.txt, s1_inputs.sha256, s1_merge_inputs.sha256 (WRITTEN AFTER THE GATE: these three files are the stage 1 commit's)")
    good += [f"{nman if nman else '<after the gate>'}         (the two manifests, each file \": OK\")", f"{tsha}   (the registered task list's SHA-256)",
             "\"the committed build record equals the build-time record\"",
             "four lines \"<SHA-256>  <binary>\" equal to the record's four sha256 lines (block B of the gate printed them)",
             f"\"{BUILD}\"", "\"env: gcc (GCC) 14.3.0 | Python ... at /home/charing/envs/hd/bin/python3 | ...\"",
             "\"ALL PASS (8 of 8 cases)\"", "\"the merge and reducer imports: OK\"",
             f"\"ls: cannot access '{LUS}': No such file or directory\" (a fresh root)"]
    return L, good


def write(args):
    CT, TL, G = load(args.chunks, args.tasks)
    c1 = CT[0]
    if c1[0] != "1": sys.exit("STOP: the chunk table's first row is not chunk 1")
    g1 = c1[4]; n1 = int(c1[3]); cid1 = c1[5]; m1 = c1[6]
    out = []
    w = out.append
    tot_core = sum(v["core_h"] for v in G.values()); tot_req = sum(v["req_h"] for v in G.values())
    w(f"##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.8 and 4.8.10) -- THE STAGE 1 RUNSHEET (KOA), Part B")
    w("(blocks B, D and S). Written by cluster/gen3_stage1_261012/make_s1_runsheet.py from the chunk table and the task list. PREPARED: typed")
    w("only after the gate (G6 (c), (d)), its verdict, the stage 1 commit and the plan author's go (G7). Chris types every command.")
    w("")
    w("WHAT IT DOES. Each chunk is one array task of one core: it runs its task lines (consecutive lines of one (cell, mass)) one after the")
    w(f"other with s1_worker_koa.sh and the campaign binary ~/harddisks_gen3_s1/hspist3/00ALLINONE ({BUILD}), checked at start")
    w("against the committed build record and the build-time record, without git. Chunk 1 runs alone first (G6 (e)); then the groups in")
    w("the registered order (sec. 4.8.8): " + ", ".join(G) + ". A group's arrays are one per time class (--time is the class), typed together;")
    w(f"their throttles add up to {CAP} at most, so at most {CAP} chunks run at once over all stage 1 arrays (sec. 4.8.8; 128 only with the plan")
    w("author's go on a pasted sinfo line with >= 10 % of the nodes idle). The next group's arrays are typed only after the previous group's")
    w("merge has ended. Every submit line of the groups (block E, block F's resubmission line, block G) is ONE guarded command line")
    w("(decision 42): it submits its sbatch lines, unchanged and in their order, only if its checks hold, and a second paste refuses.")
    w("For " + ", ".join(g for g in G if g in FG_ONE) + " blocks F and G are ONE paste, block FG (decision 45 addendum 2): F's lines, then block G's guarded")
    w("line only if F's five good conditions hold, else \"F-<g> NOT CLEAN -- the merge was not submitted\".")
    w("The merge (one per group) checks every chunk's log and folder, plans before it moves anything, reduces every")
    w("(cell, mass) and compresses the traces; the fetch brings the summaries, run records, logs, every snapshot and the whole (N 100, pi/8) cell.")
    w("KOA RULES, as before: the login node only for light commands; squeue before every sbatch; never submit a job twice; nothing in")
    w("~/harddisks, ~/harddisks_resched, ~/harddisks_resched2, ~/harddisks_gen3 or ~/harddisks_gen3_s1 changes; every job stops with")
    w("\"STOP: [<node>] ...\" on a failed start-up check before it runs anything: paste the STOP line, change nothing, wait. NO GIT AT RUN")
    w("TIME; every sbatch file carries --exclude=cn-05-08-07. THE SANDBOX LINE is the established one: srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash")
    if args.commit == "AFTER_THE_GATE":
        w("THE STAGE 1 COMMIT does not exist yet: after the gate's fetch CC copies the fetched build-time record s1_gate/BUILD_RECORD_S1.txt to")
        w("cluster/gen3_stage1_261012/build_record_f72d764.txt, writes the two manifests (make_s1_manifests.py), commits, and writes this runsheet")
        w("again with that commit (--commit): block B then names it and holds the SHA-256 lines of those three files. Nothing else changes.")
    w("")
    w("CORE-HOURS AND TIME (predicted: the design's model on the old E5 nodes; requested: the chunks' time limits; Slurm charges use):")
    w("| group | chunks | trajectories | (cell, mass) folders | predicted E5 core-h | x 1.5 | requested core-h | arrays (class: chunks %throttle) | wall at 64 cores, predicted |")
    w("|---|---|---|---|---|---|---|---|---|")
    for g, v in G.items():
        ar = arrays(g, v, skip=(1,) if g == g1 else ())
        wall = max(a[3] / a[2] for a in ar) if ar else 0.0
        w(f"| {g} | {v['nch']} | {v['traj']} | {v['pairs']} | {v['core_h']:,.1f} | {1.5 * v['core_h']:,.1f} | {v['req_h']:,.0f} | "
          + "; ".join(f"{hm(a[0])[:-3]}: {a[1]} %{a[2]}" for a in ar) + f" | {wall:.1f} h |")
    w(f"| all | {sum(v['nch'] for v in G.values())} | {sum(v['traj'] for v in G.values())} | {sum(v['pairs'] for v in G.values())} | "
      f"{tot_core:,.1f} | {1.5 * tot_core:,.1f} | {tot_req:,.0f} | | |")
    w(f"(group {g1}'s arrays leave out chunk 1, run alone in block C)")
    w("")
    w("S1 BLOCK A (KOA login node; opens the sandbox shell; then wait for the cn-... prompt):")
    w("squeue -u charing"); w("srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash"); w("")
    w("  Expected: squeue lists no gen3 job. Paste block B only at the cn-... prompt.")
    w("")
    L, good = block_b(args.commit, HS)
    w("S1 BLOCK B (in the sandbox shell: the update in place of the scripts clone to the stage 1 commit, its checks, the unit test):")
    for l in L: w(l)
    w("")
    w("  THE LINES THAT SHOW IT IS GOOD (paste everything):")
    for l in good: w("    " + l)
    w("  STOP RULE: anything else: type exit if the shell is still open, paste it, submit nothing, wait.")
    w("")
    w(f"S1 BLOCK C (KOA login node, after block B's exit; chunk 1 alone, G6 (e): {n1} task lines of {cid1}, M {m1}; 1 core; time limit {hm(int(c1[9]))[:-3]}):")
    w("squeue -u charing")
    w(f"cd ~/harddisks_gen3_k2/hspist3 && sbatch --job-name=gen3-s1-{g1} --time={hm(int(c1[9]))} --array=1 {P}/s1_koa.sbatch")
    w("")
    w("  Expected: squeue lists no gen3 job; \"Submitted batch job <its number>\". Submit it ONCE. Predicted "
      f"{float(c1[7]) / 60:.1f} min on an E5 node.")
    w("")
    c = f"{LUS}/{g1}/chunks/chunk_1/stage1/gen3/{cid1}/m_{m1}"
    w(f"S1 BLOCK D (KOA login node, when squeue -u charing lists no gen3-s1-{g1} job; every line only reads; paste everything):")
    w("squeue -u charing")
    w("cd ~/harddisks_gen3_k2/hspist3")
    w(f"cat logs/gen3-s1-{g1}_*_1.out")
    w(f"c={c}; echo \"traces $(ls $c/wall_x_positions_* | wc -l) psi6 $(ls $c/psi6_t_*.csv.gz | wc -l) virial $(ls $c/virial_blocks_*.csv.gz | wc -l) "
      f"snapshots $(ls $c/snapshots_* 2>/dev/null | wc -l) runs $(grep -c '^##RUN' $c/run.log) clean $(grep -c '^\\[EDMD3-HEALTH\\] .*: clean=1 ' $c/run.log) build '$(cat $c/.build_git)'\"")
    w("grep -h '^\\[EDMD3-RUN\\]' $c/run.log | grep -o 'box=[^ ]* x [^ ]* px exact_box=[01] eta=[^ ]* N=[0-9]* N_L=[0-9]* N_R=[0-9]*' | sort | uniq -c")
    nsnap1 = sum(1 for i in range(1, n1 + 1) if float(TL[i - 1][12]) != 0.0)
    w("")
    w("  THE LINES THAT SHOW IT IS GOOD (all of them, and nothing else):")
    w("    squeue: its header line only")
    w(f"    the log: \"env: gcc (GCC) 14.3.0 | ...\"; \"chunk 1 on cn-... (...): group {g1}, lines 1-{n1}; {BUILD}\"; \"chunk 1 done\"; no line")
    w("    containing FAILED, REFUSED or STOP")
    w(f"    traces {n1} psi6 {n1} virial {n1} snapshots {nsnap1} runs {n1} clean {n1} build '{BUILD}'")
    w(f"    one line \"{n1} box={2 * float(TL[0][5]) * 24:.17g} x {float(TL[0][6]) * 24:.17g} px exact_box=1 eta=<the recorded eta> N={2 * int(TL[0][7])} N_L={TL[0][7]} N_R={TL[0][7]}\"")
    w(f"    (the box 2 L0 x H in pixels of 1/24 sigma from the requested L0 {TL[0][5]} and H {TL[0][6]}, the exact box, the recorded eta -- within")
    w("    1e-12 of the requested one, as gate G2 checked for every cell -- and the equal halves)")
    w("  STOP RULE: anything else: paste it, change nothing, wait. The groups' blocks are typed only after the plan author's go (G7).")
    for g, v in G.items():
        ar = arrays(g, v, skip=(1,) if g == g1 else ())
        R = f"{LUS}/{g}"; PFX = f"gen3-s1-{g}"
        w("")
        w(f"S1 BLOCK E-{g} (KOA login node; group {g}: {v['nch'] - (1 if g == g1 else 0)} chunks in {len(ar)} array(s), throttles adding up to "
          f"{sum(a[2] for a in ar)}; {v['traj']} trajectories; predicted {v['core_h']:,.1f} E5 core-h):")
        w("squeue -u charing"); w("sinfo -p shared -s")
        if g == g1:
            ec = [sq(PFX), head_ok(args.commit),
                  check(f"[ \"$(ls -d {R}/chunks/chunk_* 2>/dev/null)\" = \"{R}/chunks/chunk_1\" ]", f"the chunk folders of {g} are not exactly chunk_1"),
                  check(f"[ \"$(ls logs/{PFX}_*.out 2>/dev/null | wc -l)\" -eq 1 ]", f"logs/{PFX}_*.out are not exactly one log"),
                  check(f"grep -qx 'chunk 1 done' logs/{PFX}_*_1.out", "the log of chunk 1 has no line 'chunk 1 done'")]
        else:
            ec = [sq(PFX), head_ok(args.commit), check(f"[ ! -e {R} ]", f"{R} exists"),
                  check(f"! ls logs/{PFX}_*.out >/dev/null 2>&1", f"logs/{PFX}_*.out exist")]
        w(guarded(ec, [a[4] for a in ar], f"{g}: ALL {len(ar)} ARRAY LINE(S) SUBMITTED", f"{g}: THIS PASTE DID NOT SUBMIT ALL {len(ar)} ARRAY LINE(S)"))
        w("")
        w(f"  Expected: squeue lists no gen3-s1 job ({'chunk 1 has ended and the plan author has read block D' if g == g1 else 'the previous group has ended: its block H was clean'}); the sinfo line, for the record;")
        w(f"  then \"CHECKS GOOD\", {len(ar)} line(s) \"Submitted batch job <number>\" and \"{g}: ALL {len(ar)} ARRAY LINE(S) SUBMITTED\". The arrays: "
          + "; ".join(f"{hm(a[0])[:-3]} h: {a[1]} chunks, at most {a[2]} at once" for a in ar) + ".")
        w(f"  THE GUARD (decision 42): the last line is ONE command line; it submits the {len(ar)} sbatch line(s) inside it, unchanged and in their order,")
        w(f"  only if squeue lists no {PFX} job (a failed squeue refuses), the scripts clone is at {args.commit[:7]}, and "
          + ("exactly chunk 1's folder and its one log with \"chunk 1 done\" exist (block C ran chunk 1)." if g == g1 else f"no {R} and no logs/{PFX}_*.out exist."))
        w(f"  A second paste prints \"REFUSED: ...\" and \"{g}: THIS PASTE DID NOT SUBMIT ALL {len(ar)} ARRAY LINE(S)\" and submits nothing. If it prints that line")
        w("  after some \"Submitted batch job\" lines, a submission failed half way: paste everything, type nothing more, wait.")
        w("")
        one = g in FG_ONE                                     # ##CHRIS 2026-10-13 (plan-author decision 45 addendum 2): F and G as one paste
        if one:
            w(f"S1 BLOCK FG-{g} (KOA login node, when squeue -u charing lists no {PFX} job; blocks F and G as ONE paste, decision 45 addendum 2: F's "
              "lines only read, the last line submits the merge only if F is clean; the merge: one core, at most 6 h; paste everything, once):")
        else:
            w(f"S1 BLOCK F-{g} (KOA login node, when squeue -u charing lists no {PFX} job; every line only reads; paste everything):")
        w("squeue -u charing"); w("cd ~/harddisks_gen3_k2/hspist3")
        w(f"ls logs/{PFX}_*.out | wc -l")
        w(f"grep -lx 'chunk [0-9]* done' logs/{PFX}_*.out | wc -l")
        w(f"cat logs/{PFX}_*.out | grep -c -e FAILED -e REFUSED")
        w(f"grep -h '^chunk [0-9]* on ' logs/{PFX}_*.out | sed 's/: group.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c")
        w(f"bash {P}/s1_chunk_logs_koa.sh logs {PFX} {R}/chunks {P}/chunks_stage1_koa_261012.tsv {g}")
        w(f_awk(g, args.commit, P))
        if one: w(fg_line(g, v["nch"], R, PFX, args.commit))
        w("")
        w("  THE LINES THAT SHOW IT IS GOOD:")
        w("    squeue: its header line only")
        w(f"    {v['nch']}        (the chunk logs: one per chunk{', chunk 1 included' if g == g1 else ''})")
        w(f"    {v['nch']}        (each chunk once with its line \"chunk <j> done\")")
        w("    0          (no line with FAILED or REFUSED; a run the worker refused prints FAILED: it is excluded and counted by the verdict, never")
        w("               rerun -- paste it, it does not stop the merge)")
        w(f"    one line per CPU model, the counts adding up to {v['nch']}")
        w(f"    \"chunk logs: {v['nch']}; chunks done: {v['nch']} of {v['nch']} (group {g}; each with its folder and exactly one log with 'chunk <j> done'); earlier attempts: 0 stopped in the start-up checks, 0 started and not finished\"")
        w(f"    \"group {g} without a done log: stopped once 0, stopped twice or more 0, started and not finished 0, no log 0\" and \"nothing to resubmit\"")
        if one:
            w(f"    then \"CHECKS GOOD\", \"Submitted batch job <its number>\" and \"{g}: THE MERGE SUBMITTED\" (block G's guarded line, inside the last line)")
            w(f"  THE LAST LINE (decision 45 addendum 2) types block G's guarded line, unchanged, only if all five hold: {v['nch']} chunk logs, {v['nch']} logs")
            w("  with their done line, 0 lines with FAILED or REFUSED, the chunk check exits 0, and \"nothing to resubmit\". Otherwise it prints")
            w(f"  \"F-{g} NOT CLEAN -- the merge was not submitted\" and submits nothing: paste everything and wait. G's own guard (decision 42) still holds:")
            w(f"  it submits only if squeue lists no {PFX} and no gen3-s1-merge-{g} job (a failed squeue refuses), the scripts clone is at {args.commit[:7]}")
            w(f"  and {R}/merged does not exist; so a second paste prints \"REFUSED: ...\" and \"{g}: THIS PASTE DID NOT SUBMIT THE MERGE\". The merge job")
            w("  checks for itself, before it moves anything, and stops with exit 2 and NOTHING moved on any failure (its inputs, the task list, no")
            w("  job of the group pending or running, the chunk check, merged/ absent: exit 3); then the plan, and only a clean plan is executed.")
        w(f"  IF the {'command before the last' if one else 'last command'} prints \"THE RESUBMISSION LINE(S)\" (sbatch --job-name=... --time=... --array=... --exclude=cn-05-08-07,...) and then")
        w("  \"THE GUARDED RESUBMISSION LINE\" and one long line: those chunks stopped in their start-up checks and ran nothing (decision 28, as")
        w("  decision 25: submitted once more, and listed). Paste everything and wait. After the plan author's go, type ONLY the one guarded line,")
        w("  once (not the plain sbatch lines above it): it submits all of them, unchanged and in their order, only if squeue lists no")
        w(f"  {PFX} job, the scripts clone is at {args.commit[:7]} and no listed chunk has a log with its done line; it prints \"CHECKS GOOD\", the")
        w(f"  \"Submitted batch job\" lines and \"{g}: ALL <n> RESUBMISSION LINE(S) SUBMITTED\"; a second paste prints \"REFUSED: ...\" and submits")
        if one:
            w(f"  nothing. The last line has then printed \"F-{g} NOT CLEAN -- the merge was not submitted\". When squeue lists no {PFX} job after the")
            w(f"  resubmission, the plan author gives the next step: the log count then includes the earlier attempts, so block FG-{g} prints NOT CLEAN by design.")
        else:
            w(f"  nothing. When squeue lists no {PFX} job, block F-{g} again (the log counts then include the earlier attempts).")
        w("  STOP RULE: a chunk \"started and did not finish\", \"stopped twice or more\", \"no log\", \"NOT DONE\", or other counts: paste, change")
        w("  nothing, wait.")
        w("")
        if not one:
            w(f"S1 BLOCK G-{g} (KOA login node; the merge of group {g}; one core, at most 6 h):")
            w("squeue -u charing")
            w(g_line(g, R, PFX, args.commit))
            w("")
            w(f"  Expected: no {PFX} job in squeue; \"CHECKS GOOD\", \"Submitted batch job <its number>\" and \"{g}: THE MERGE SUBMITTED\". THE GUARD (decision 42):")
            w(f"  it submits only if squeue lists no {PFX} and no gen3-s1-merge-{g} job (a failed squeue refuses), the scripts clone is at {args.commit[:7]} and")
            w(f"  {R}/merged does not exist; a second paste prints \"REFUSED: ...\" and \"{g}: THIS PASTE DID NOT SUBMIT THE MERGE\". The merge job checks for itself, before it moves anything,")
            w("  and stops with exit 2 and NOTHING moved on any failure (its inputs, the task list, no job of the group pending or running, the")
            w("  chunk check, merged/ absent: exit 3); then the plan, and only a clean plan is executed.")
            w("")
        mv = 3 * v["traj"] + v["snap"]
        w(f"S1 BLOCK H-{g} (KOA login node, when squeue -u charing lists no gen3-s1-merge-{g} job; paste everything):")
        w("squeue -u charing")
        w(f"L=$(ls -t ~/harddisks_gen3_k2/hspist3/logs/gen3-s1-merge-{g}_*.out | head -1) && echo \"$L\"")
        w("grep -v -e '^| stage1/gen3/' -e '^reduce_B stage1/gen3/.*: exit 0 ' -e 'trajectories reduced *$' -e '^reduce_s1 stage1/gen3/.*: exit 0 ' -e ': [0-9]* traces compressed$' \"$L\"")
        w("echo \"destinations $(grep -c '^| stage1/gen3/' \"$L\"); reduce_B exit 0: $(grep -c '^reduce_B stage1/gen3/.*: exit 0 ' \"$L\"); reduce_s1 exit 0: $(grep -c '^reduce_s1 stage1/gen3/.*: exit 0 ' \"$L\"); compressed: $(grep -c ': [0-9]* traces compressed$' \"$L\")\"")
        w("grep -h '^reduce_s1 ' \"$L\" | grep -o '([0-9]* fit failures)' | sort | uniq -c")
        w(f"cat {R}/merged/.merge_complete")
        w(f"ls {R}/job_logs | wc -l")
        w("")
        w("  THE LINES THAT SHOW IT IS GOOD:")
        w("    squeue: its header line only; the log's name")
        w("    the log without its repeated per-folder lines:")
        w("      - \"env: gcc (GCC) 14.3.0 | ...\"")
        w(f"      - \"merge of {R} by the scripts of ~/harddisks_gen3_k2 (no git at run time; both input manifests: OK); the scripts it runs:\"")
        w("        and four SHA-256 lines (merge_s1_koa.py, s1_chunk_logs_koa.sh, reduce_s1_koa.py, reduce_B.py)")
        w(f"      - \"group {g}: {v['nch']} chunks of the chunk table; no {PFX} job pending or running\"")
        w(f"      - \"chunk logs in {R}/job_logs/ (copied, no overwrite): {v['nch']}\" (more if chunks were resubmitted)")
        w(f"      - \"chunk logs: {v['nch']}; chunks done: {v['nch']} of {v['nch']} (group {g}; ...); earlier attempts: 0 stopped in the start-up checks, 0 started and not finished\"")
        w("      - \"FAILED lines in the chunk logs (logs/): 0\"")
        w(f"      - \"the task list: {len(TL)} lines; group {g}: {v['nch']} chunks of the chunk table ({v['first']}-{v['last']}); top-level folder stage1\"")
        w(f"      - \"scan: {v['nch']} chunks, {mv + 2 * v['nch'] + v['traj'] + v['nch']} plan entries, 0 problem(s); plan written to ...\"; the plan table's header")
        w(f"      - \"plan entries by treatment: APPEND {2 * v['nch']}, EQUAL {v['nch']}, MOVE {mv}, PER CHUNK {v['traj']}\" (if no run failed)")
        w("      - one line \"names inside the moved run folders .done_runs/run<r> (...): ...\"")
        w(f"      - \"merged {v['nch']} chunks: {mv} files moved, {2 * v['nch']} files appended, {v['traj']} kept per chunk (0 failed-run folders), 0 incomplete")
        w(f"        run folders, 0 locks or temp files left in place, {v['nch']} build files equal; build '{BUILD}'\"")
        w(f"      - \"merge done; next: the fetch (fetch_s1_koa.sh {g}) and the tables on the Mac\" and \"this job's log is copied to {R}/job_logs/ as the last command\"")
        w(f"    \"destinations {v['pairs']}; reduce_B exit 0: {v['cells']}; reduce_s1 exit 0: {v['pairs']}; compressed: {v['pairs']}\"")
        w(f"    the fit-failure counts of the {v['pairs']} (cell, mass) folders, \"<n> (<k> fit failures)\" (reported; the verdict counts them)")
        w(f"    .merge_complete: the end time and the same \"merged {v['nch']} chunks: ...\" line")
        w(f"    {v['nch'] + 1}        (the job logs copied: {v['nch']} chunk logs and the merge log; more if chunks were resubmitted)")
        w("  STOP RULE: an exit other than 0, a STOP line, a problem list, other counts: paste, change nothing, wait.")
        w("")
        w(f"S1 BLOCK I-{g} (Mac, after block H-{g}; one password and Duo push):")
        w(f"cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3 && bash hspist3/{P}/fetch_s1_koa.sh {g}")
        w("")
        w(f"  Expected, last line: \"fetched into .../experiments_gen3_stage1_koa_261012/{g}: <size>; merge record: <the end time>\". Then tell CC")
        w("  (CC prints the group's INTERIM tables with s1_verdict_261012.py --data --interim). The next group's block E follows once block")
        w(f"  H-{g} is clean; the fetch and the tables run beside it.")
    open(args.out, "w").write("\n".join(out) + "\n")
    print(f"wrote {args.out}: {len(out)} lines; groups {', '.join(G)}; commit {args.commit}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", default="AFTER_THE_GATE"); ap.add_argument("--chunks", default=os.path.join(HERE, "chunks_stage1_koa_261012.tsv"))
    ap.add_argument("--tasks", default=os.path.join(HERE, "tasks_stage1_koa_261012.txt")); ap.add_argument("--out", default=os.path.join(HERE, "runsheet_stage1_261012.txt"))
    write(ap.parse_args())


if __name__ == "__main__":
    main()
