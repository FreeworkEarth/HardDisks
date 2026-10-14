#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 47 item 3): THE MAC REHEARSAL OF THE REDUCTION COMPLETION, "on the mini data with a merge
stopped half way (a pretended kill)". It reuses the stage 1 rehearsal (cluster/gen3_stage1_261012/rehearse_s1_mac.py: the build, the
substitutions, the shims, the mini list and table): group B100's mini chunks run (blocks C and E-B100: 4 chunks, 31 trajectories, 2 cells,
3 (cell, mass) folders), then the group's merge job (s1_merge.sbatch through the sbatch shim, as Chris submits it) runs several times on a
copy of the same chunks:
  REF   to its end (the reference tree);
  K1-K5 KILLED with SIGKILL (the whole process group, as Slurm's time limit kills it) the moment a chosen file appears: K1 the first
        reduce_s1 table of the first folder; K2 the first .gz of the first folder (its compression); K3 reduce_B's red_nu.csv of the
        second cell's first folder; K4 the first reduce_s1 table of the last folder; K5 the first .gz of the last folder, then a cut gzip
        written for its next trace (a .csv.gz of 64 bytes beside its .csv, as a gzip killed at its start leaves it);
  then the completion job (s1_complete.sbatch through the sbatch shim, from the clone home/harddisks_gen3_rc, as Chris would submit it).
CHECKS per killed case: the completion exits 0; the completed tree equals REF file by file (CSV and text byte for byte; .npz array by
array; .gz by its content and its bytes), apart from the job logs, the plan file's timestamped name, .merge_complete's time line, the
new .reduction_completed and the cut gzips kept aside; every file that already equalled REF before the completion keeps its mtime and
inode (not rewritten); no file disappears (a cut gzip is kept under <name>.partial_<job id>).
REFUSALS: REF's tree (nothing to complete: exit 0, nothing changed); K1 a second time (exit 3: .reduction_completed); K2 with a pretended
running merge job in squeue (exit 2, nothing changed); a tree whose merged/ has no .merge_complete (exit 2, nothing changed).
usage (from hspist3/): python3 cluster/gen3_s1_complete_261013/rehearse_complete_mac.py --root <empty dir> --pkg <engine-gen3 commit>
"""
import argparse, gzip, hashlib, os, re, shutil, signal, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage1_261012"))
import rehearse_s1_mac as G
RES = []
GRP = "B100"; C1, C2 = "B_N0100_eta0.050000", "B_N0100_eta0.392699"


def res(step, crit, ok, ev):
    RES.append((step, crit, "PASS" if ok else "**FAIL**", ev)); print(f"    -> {step}: {'PASS' if ok else '**FAIL**'} ({ev})", flush=True)


def digest(p):
    if p.endswith(".npz"):
        z = np.load(p); h = hashlib.sha256()
        for k in sorted(z.files): h.update(k.encode()); a = np.asarray(z[k]); h.update(str(a.dtype).encode()); h.update(a.tobytes())
        return "npz:" + h.hexdigest()
    if p.endswith(".gz"):
        try: return "gz:" + hashlib.sha256(gzip.open(p, "rb").read()).hexdigest()
        except Exception: return "gz-unreadable:" + hashlib.sha256(open(p, "rb").read()).hexdigest()
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def tree(root):
    out = {}
    for dp, dn, fn in os.walk(root):
        for f in fn:
            p = os.path.join(dp, f); st = os.stat(p)
            out[os.path.relpath(p, root)] = dict(d=digest(p), raw=hashlib.sha256(open(p, "rb").read()).hexdigest(), mt=st.st_mtime_ns, ino=st.st_ino)
    return out


def excluded(rel):
    return (rel.startswith("job_logs/") or re.match(r"^merge_plan_\d{8}_\d{6}\.tsv$", rel) or rel == "merged/.merge_complete"
            or rel == "merged/.reduction_completed" or ".partial_" in rel)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--pkg", required=True); a = ap.parse_args()
    R = os.path.abspath(a.root)
    if not os.path.isdir(R) or os.listdir(R): sys.exit(f"STOP: {R} must be an existing empty folder")
    PKG = G.sh(["git", "-C", HS, "rev-parse", a.pkg]).stdout.strip()
    for d in ("lustre", "logs", "cases"): os.makedirs(os.path.join(R, d))
    print(f"# the reduction completion rehearsal (decision 47 item 3): package {PKG} (engine-gen3), root {R}\n", flush=True)
    b = G.sh([sys.executable, os.path.join(HS, "cluster/gen3_s1_gate_261012/rehearse_s1_build_mac.py"), "--root", os.path.join(R, "build"),
              "--commit", PKG], cwd=HS, log=os.path.join(R, "logs", "build_rehearsal.out"), check=False)
    bt = open(os.path.join(R, "logs", "build_rehearsal.out")).read()
    res("the build (rehearse_s1_build_mac.py)", "BUILD REHEARSAL: PASS", b.returncode == 0 and "BUILD REHEARSAL: PASS" in bt, [l for l in bt.splitlines() if l.startswith("BUILD REHEARSAL")][:1])
    home = os.path.join(R, "build", "home"); k2 = os.path.join(home, "harddisks_gen3_k2"); hs = os.path.join(k2, "hspist3")
    shim = G.shims(R); E = G.env_for(R, home, shim); os.makedirs(os.path.join(hs, "logs"), exist_ok=True)
    mini, rows, rew, man = G.substitute(R, hs, home)
    sc = os.path.join(hs, "cluster", "gen3_s1_complete_261013", "s1_complete.sbatch"); t = open(sc).read()
    open(sc, "w").write(t.replace(G.LUSTRE, os.path.join(R, "lustre")))
    mm = G.sh([sys.executable, os.path.join("cluster", "gen3_s1_complete_261013", "make_complete_manifest.py")], cwd=hs).stdout.strip()
    G.sh(["git", "-C", k2, "add", "-u"]); G.sh(["git", "-C", k2, "add", f"hspist3/{G.P}/build_record_f72d764.txt", f"hspist3/{G.P}/s1_inputs.sha256", f"hspist3/{G.P}/s1_merge_inputs.sha256"])
    G.sh(["git", "-C", k2, "-c", "user.name=rehearsal", "-c", "user.email=rehearsal@localhost", "commit", "-q", "-m", "REHEARSAL stand-in: mini list and table, Lustre path, manifests"])
    Q = G.sh(["git", "-C", k2, "rev-parse", "HEAD"]).stdout.strip()
    rc_ = os.path.join(home, "harddisks_gen3_rc"); G.sh(["git", "clone", "-q", "--shared", k2, rc_]); G.sh(["git", "-C", rc_, "-c", "advice.detachedHead=false", "checkout", "-q", Q])
    rcs = os.path.join(rc_, "hspist3"); os.makedirs(os.path.join(rcs, "logs"))
    chk = G.sh(["bash", "-c", "shasum -a 256 -c cluster/gen3_stage1_261012/s1_inputs.sha256 cluster/gen3_stage1_261012/s1_merge_inputs.sha256 cluster/gen3_s1_complete_261013/s1_complete_inputs.sha256 | grep -c ': OK$'"], cwd=rcs, check=False).stdout.strip()
    res("the substitutions; the completion's Lustre path and manifest; the clone home/harddisks_gen3_rc", "the three manifests pass in the rc clone",
        chk.isdigit() and int(chk) > 21, f"{mm}; OK lines {chk}; stand-in {Q[:12]}")
    rs = os.path.join(R, "logs", "runsheet_stage1_mini.txt"); G.sh([sys.executable, os.path.join(G.P, "make_s1_runsheet.py"), "--commit", Q, "--out", rs], cwd=hs)
    B = G.blocks_of(rs)
    rc, out = G.run_block(R, "C", B["C"], E); rc2, out2 = G.run_block(R, f"E-{GRP}", B[f"E-{GRP}"], E)
    groot = os.path.join(R, "lustre", "harddisks_gen3", "stage1", GRP)
    ok = rc == 0 and rc2 == 0 and f"{GRP}: ALL" in out2 and sorted(os.listdir(os.path.join(groot, "chunks"))) == ["chunk_1", "chunk_2", "chunk_3", "chunk_4"]
    res(f"group {GRP}'s mini chunks (blocks C and E-{GRP})", "4 chunk folders, every chunk exit 0", ok, G.last_exit(R))
    pre = os.path.join(R, "cases", "pre_merge"); shutil.copytree(groot, pre, symlinks=True, copy_function=shutil.copy2)
    EX = os.path.join(R, "logs", "sbatch_shim_exits.txt")
    case_no = [0]

    def restore(name_prev):
        if os.path.exists(groot): os.rename(groot, os.path.join(R, "cases", name_prev))
        shutil.copytree(pre, groot, symlinks=True, copy_function=shutil.copy2)

    def merge(kill_when=None, tag=""):
        """the merge job through the sbatch shim; killed (its whole process group, SIGKILL) the moment kill_when() is true"""
        p = subprocess.Popen(["bash", "-c", f"sbatch --job-name=gen3-s1-merge-{GRP} cluster/gen3_stage1_261012/s1_merge.sbatch {GRP}"], cwd=hs, env=E,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
        killed = False
        while p.poll() is None:
            if kill_when and kill_when():
                os.killpg(p.pid, signal.SIGKILL); killed = True; break
            time.sleep(0.005)
        p.wait()
        ml = sorted((f for f in os.listdir(os.path.join(hs, "logs")) if f.startswith(f"gen3-s1-merge-{GRP}_")), key=lambda f: os.stat(os.path.join(hs, "logs", f)).st_mtime)[-1]
        return killed, os.path.join(hs, "logs", ml)

    def complete(mlog, env=None, tag=""):
        n0 = len(open(EX).read().splitlines())
        q = subprocess.run(["bash", "-c", f"sbatch --time=24:00:00 --job-name=gen3-s1-complete-{GRP} cluster/gen3_s1_complete_261013/s1_complete.sbatch {GRP} {mlog}"],
                           cwd=rcs, env=env or E, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        line = [l for l in open(EX).read().splitlines()[n0:] if " submit " not in l]
        jl = sorted((f for f in os.listdir(os.path.join(rcs, "logs")) if f.startswith(f"gen3-s1-complete-{GRP}_")), key=lambda f: os.stat(os.path.join(rcs, "logs", f)).st_mtime)[-1]
        shutil.copy2(os.path.join(rcs, "logs", jl), os.path.join(R, "logs", f"completion_{tag}.out"))
        return (int(line[-1].split()[-1]) if line else None), open(os.path.join(rcs, "logs", jl)).read()

    # ---- REF: the merge to its end (on the chunks as block E left them; pre_merge is their copy for the other cases)
    k, mlog = merge(None, "REF"); ref = tree(groot)
    res("REF: the merge to its end", ".merge_complete; every folder reduced and compressed", not k and os.path.exists(os.path.join(groot, "merged", ".merge_complete"))
        and not [r for r in ref if re.search(r"wall_x_positions_.*\.csv$", r)], f"{len(ref)} files; log {os.path.basename(mlog)}")
    before = tree(groot); x, jt = complete(mlog, tag="REF")
    after = tree(groot)
    res("REF: the completion on a finished merge", "exit 0, 'nothing to complete', nothing changed but the logs copied", x == 0 and "nothing to complete" in jt
        and {r: v for r, v in after.items() if not r.startswith("job_logs/")} == {r: v for r, v in before.items() if not r.startswith("job_logs/")}, f"exit {x}")
    refd = {r: v["d"] for r, v in ref.items() if not excluded(r)}
    refraw = {r: v["raw"] for r, v in ref.items() if not excluded(r)}
    plan_ref = [v["d"] for r, v in ref.items() if re.match(r"^merge_plan_\d{8}_\d{6}\.tsv$", r)]
    mc_ref = open(os.path.join(groot, "merged", ".merge_complete")).read().splitlines()[1:]
    prev = "REF"
    S = lambda *p: os.path.join(groot, "merged", "stage1", "gen3", *p)
    first_gz = lambda d: lambda: os.path.isdir(d) and any(f.endswith(".csv.gz") and f.startswith("wall_x_positions_") for f in os.listdir(d))
    cases = [("K1", "the first reduce_s1 table of the first folder", lambda: os.path.exists(S(C1, "m_50", "nud_runs.csv"))),
             ("K2", "the first .gz of the first folder", first_gz(S(C1, "m_50"))),
             ("K3", "reduce_B's red_nu.csv of the second cell's first folder", lambda: os.path.exists(S(C2, "m_50", "red_nu.csv"))),
             ("K4", "the first reduce_s1 table of the last folder", lambda: os.path.exists(S(C2, "m_500", "nud_runs.csv"))),
             ("K5", "the first .gz of the last folder, then a cut gzip for its next trace", first_gz(S(C2, "m_500")))]
    for name, what, cond in cases:
        restore(prev); prev = name
        k, mlog = merge(cond, name)
        if name == "K5":
            md = S(C2, "m_500"); un = sorted(f for f in os.listdir(md) if f.startswith("wall_x_positions_") and f.endswith(".csv") and not os.path.exists(os.path.join(md, f + ".gz")))
            if un:
                full = subprocess.run(["gzip", "-9", "-c", os.path.join(md, un[0])], stdout=subprocess.PIPE).stdout
                open(os.path.join(md, un[0] + ".gz"), "wb").write(full[:64])
        st = subprocess.run([sys.executable, os.path.join(rcs, "cluster", "gen3_s1_complete_261013", "complete_s1_koa.py"), groot, "--plan-only"],
                            cwd=rcs, env=E, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
        todo = [l for l in st.splitlines() if l.startswith("| stage1/") and "nothing (finished)" not in l]
        print(f"  {name}: killed {k} at '{what}'; to complete: " + "; ".join(re.sub(r"\s+", " ", l) for l in todo), flush=True)
        if name == "K2":
            q9 = os.path.join(R, "logs", "pretend_squeue_merge.txt"); open(q9, "w").write(f"940990 RUNNING gen3-s1-merge-{GRP}\n")
            b9 = tree(groot); x9, j9 = complete(mlog, env=G.env_for(R, home, shim, REHEARSAL_SQUEUE=q9), tag="K2_pending_merge")
            res("K2: the completion with a pretended running merge job", "exit 2, nothing changed", x9 == 2 and tree(groot) == b9 and "a merge or another completion" in j9, f"exit {x9}")
        bt_ = tree(groot); x, jt = complete(mlog, tag=name); at = tree(groot)
        got = {r: v["d"] for r, v in at.items() if not excluded(r)}
        same = got == refd; rawsame = {r: v["raw"] for r, v in at.items() if not excluded(r)} == refraw
        untouched = [r for r, v in bt_.items() if refd.get(r) == v["d"]]
        rewritten = [r for r in untouched if r in at and (at[r]["mt"], at[r]["ino"]) != (bt_[r]["mt"], bt_[r]["ino"])]
        vanished = [r for r in bt_ if r not in at and not any(q.startswith(r + ".partial_") and at[q]["raw"] == bt_[r]["raw"] for q in at)
                    and not (r + ".gz" in at and at[r + ".gz"]["d"] == "gz:" + bt_[r]["raw"])]      # a trace gzipped: its .csv is in the .gz
        extra = [r for r in at if r not in bt_ and r not in refd and not excluded(r)]
        plan_ok = [v["d"] for r, v in at.items() if re.match(r"^merge_plan_\d{8}_\d{6}\.tsv$", r)] == plan_ref
        mc_ok = open(os.path.join(groot, "merged", ".merge_complete")).read().splitlines()[1:] == mc_ref
        kept = [r for r in at if ".partial_" in r]
        okc = [x == 0, "completion done" in jt, same, not rewritten, not vanished, not extra, plan_ok, mc_ok, os.path.exists(os.path.join(groot, "merged", ".reduction_completed"))]
        if name == "K5": okc.append(len(kept) >= 1)
        res(f"{name}: killed at {what}; then the completion", "exit 0; the tree equals REF file by file; the files already equal to REF not rewritten; "
            "nothing vanished; no extra file" + ("; the cut gzip kept aside" if name == "K5" else ""), all(okc),
            f"killed {k}; checks {okc}; {len(untouched)} files already final and untouched; raw bytes equal to REF: {rawsame}; kept aside: {kept}")
        if name == "K1":
            b8 = tree(groot); x8, j8 = complete(mlog, tag="K1_again")
            res("K1: the completion a second time", "exit 3 (.reduction_completed), nothing changed but its log", x8 == 3
                and {r: v for r, v in tree(groot).items() if not r.startswith("job_logs/")} == {r: v for r, v in b8.items() if not r.startswith("job_logs/")}, f"exit {x8}")
    # ---- a merge stopped during its moves: merged/ without .merge_complete
    restore(prev); prev = "K6"
    os.makedirs(os.path.join(groot, "merged", "stage1")); b6 = tree(groot)
    x6, j6 = complete(os.path.join(R, "logs", "no_such_merge_log.out"), tag="K6")
    res("K6: merged/ without .merge_complete (a merge stopped during its moves)", "exit 2, nothing changed", x6 == 2 and tree(groot) == b6
        and ".merge_complete is missing" in j6, f"exit {x6}")
    os.rename(groot, os.path.join(R, "cases", prev))
    print("\n| step | criterion | result | evidence |\n|---|---|---|---|")
    for r in RES: print("| " + " | ".join(str(x) for x in r) + " |")
    nb = sum(r[2] != "PASS" for r in RES)
    print(f"\nRESIDUAL DIFFERENCES: the shims (Slurm commands, module, gcc --version, lscpu, hostname, sha256sum); the Mac binaries; the mini list "
          f"and table; the Lustre path; the kills by SIGKILL to the shim's process group at a file's appearance (on KOA: Slurm's time limit)")
    print(f"\nCOMPLETION REHEARSAL: {'PASS' if nb == 0 else 'FAIL'} ({len(RES) - nb} of {len(RES)})")
    print(f"data under {R}: {G.sh(['du', '-sh', R]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
