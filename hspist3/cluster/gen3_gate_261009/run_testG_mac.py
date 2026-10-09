#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.15, Test G): the Mac runner. Runs a registered task list with testG_worker.sh, at most
--jobs processes at once (programme rule 1: <= 12), IN TASK ORDER (gen2 and gen3 of one seed are consecutive lines, so the engines
run interleaved), each trajectory in its own folder; then reduces every cell (reduce_B.py, unchanged) and compresses the cell's
traces losslessly (gzip -9; the SHA-256 of each uncompressed trace in the cell's .sha256_uncompressed first). Refuses to start if
the task list's SHA-256 or the binary's version line or SHA-256 differ from the expected ones. Never reruns a failed task.
KEEP_EV: the first 10 AF seeds per engine (task order) keep their event log on disk (testG_worker.sh); the others stream it.
usage (from hspist3/): python3 cluster/gen3_gate_261009/run_testG_mac.py --tasks <list> --expect-sha <sha256 of the list>
          --bin <frozen 00ALLINONE> --expect-build "<first line of --version>" --expect-bin-sha <sha256> [--data <root>] [--jobs 12]
"""
import argparse, collections, hashlib, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
LOC = os.path.join(HS, "experiments_gen3_gate_261009")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True); ap.add_argument("--expect-sha", required=True)
    ap.add_argument("--bin", required=True); ap.add_argument("--expect-build", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--data", default=LOC); ap.add_argument("--jobs", type=int, default=12); ap.add_argument("--keep-ev", type=int, default=10)
    a = ap.parse_args()
    if a.jobs > 12: sys.exit("STOP: at most 12 processes (programme rule 1)")
    tasks, binp, data = os.path.abspath(a.tasks), os.path.abspath(a.bin), os.path.abspath(a.data)
    if sha(tasks) != a.expect_sha: sys.exit(f"STOP: task list SHA-256 {sha(tasks)} != registered {a.expect_sha}")
    if sha(binp) != a.expect_bin_sha: sys.exit(f"STOP: binary SHA-256 {sha(binp)} != {a.expect_bin_sha}")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    if build != a.expect_build: sys.exit(f"STOP: binary version line '{build}' != '{a.expect_build}'")
    os.makedirs(data, exist_ok=True)
    log = open(os.path.join(data, f"runner_{os.path.basename(tasks)}.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} run_testG_mac.py: tasks {os.path.relpath(tasks, HS)} (SHA-256 {a.expect_sha}); binary {binp} "
        f"(SHA-256 {a.expect_bin_sha}; {build}); data {data}; jobs {a.jobs}")
    say(subprocess.run(["df", "-h", data], capture_output=True, text=True).stdout.strip())
    lines = [l.split() for l in open(tasks) if l.strip()]
    af_seen = collections.Counter()
    env0 = dict(os.environ, HD_BIN=binp, HD_DATA=data, HD_BUILD=build)
    env0.pop("KEEP_EV", None)
    run, done, t0 = [], [], time.time()
    def reap():
        for j in list(run):
            p, i, ts = j
            if p.poll() is not None:
                run.remove(j); out = p.stdout.read().decode(errors="ignore").strip()
                done.append((i, p.returncode))
                if p.returncode != 0 or out: say(f"task {i + 1}: exit {p.returncode} {out}")
    for i, f in enumerate(lines):
        while len(run) >= a.jobs:
            reap(); time.sleep(0.2)
        env = dict(env0)
        if f[0] == "AF":
            eng = f[10]; af_seen[eng] += 1
            if af_seen[eng] <= a.keep_ev: env["KEEP_EV"] = "1"
            args = f[:10] + [eng]
        else:
            args = f[:11]
        p = subprocess.Popen(["bash", os.path.join(HERE, "testG_worker.sh")] + args, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        run.append((p, i, time.time()))
        if (i + 1) % 200 == 0: say(f"  started {i + 1} of {len(lines)} tasks, {time.time() - t0:.0f} s")
    while run: reap(); time.sleep(0.2)
    nf = sum(rc != 0 for _, rc in done)
    say(f"# all {len(lines)} tasks finished in {time.time() - t0:.0f} s; non-zero exits: {nf}")
    cells = sorted({os.path.dirname(os.path.join(data, f[1])) for f in lines if f[0] == "B"})
    for c in cells:
        r = subprocess.run([sys.executable, os.path.join(CL, "confinement_20261013", "reduce_B.py"), c], capture_output=True, text=True)
        say(f"reduce_B {os.path.relpath(c, data)}: exit {r.returncode} {r.stdout.strip()} {r.stderr.strip()[-300:]}")
        if r.returncode != 0: continue
        for m in sorted(d for d in os.listdir(c) if d.startswith("m_")):
            md = os.path.join(c, m); tr = sorted(x for x in os.listdir(md) if x.startswith("wall_x_positions_") and x.endswith(".csv"))
            with open(os.path.join(md, ".sha256_uncompressed"), "a") as fh:
                for x in tr: fh.write(f"{sha(os.path.join(md, x))}  {x}\n")
            for x in tr: subprocess.run(["gzip", "-9", os.path.join(md, x)], check=True)
            say(f"  {os.path.relpath(md, data)}: {len(tr)} traces compressed (SHA-256 of each uncompressed trace in .sha256_uncompressed)")
    say(subprocess.run(["df", "-h", data], capture_output=True, text=True).stdout.strip())
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")


if __name__ == "__main__":
    main()
