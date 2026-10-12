#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 34, the rest; decision 38 item 4): the records of the evidence, as decision 19's
(D19_SHA256_FILES.txt, D19_RUN_RECORDS.txt): every file of the bin and evidence folders with its SHA-256 and size; the key lines of every
run log (the logs are .gitignored) and the runner's log in full; and the list of the files under 1 MB that git does not ignore (the ones
committed; the larger ones are recorded here only).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d34s_261012/s_records.py <bin dir> <evidence dir>   (writes D34S_SHA256_FILES.txt, D34S_RUN_RECORDS.txt, D34S_COMMIT_LIST.txt)
"""
import hashlib, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = re.compile(r"^(\[EDMD3[-\w]*\]|\[EDMD-CONTACT\]|\[GUI-STATS\]|\[GEN3-GUI\]|\[GUI-LIMIT\]|\[GUI-STOP\]|\[GUI\] |\[GUI-SHOT\]|\[GUI-KEYS\]|\[ENERGY\]|\[DEMO\] screenshot|"
                 r"Scene reset queued|.*spring: k = |.*Energy measurement initialized|.*Invalid wall_x|STOP|.*Failed|\[warn\])")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    B, O = (os.path.abspath(x) for x in sys.argv[1:3])
    files = []
    for top in (B, O):
        for r, _, fs in os.walk(top):
            for f in fs:
                p = os.path.join(r, f)
                if os.path.islink(p): continue
                files.append(p)
    files.sort(key=lambda p: os.path.relpath(p, HERE))
    with open(os.path.join(HERE, "D34S_SHA256_FILES.txt"), "w") as fh:
        for p in files: fh.write(f"{sha(p)}  {os.path.getsize(p):10d}  {os.path.relpath(p, HERE)}\n")
    with open(os.path.join(HERE, "D34S_RUN_RECORDS.txt"), "w") as fh:
        fh.write("# Every key line ([EDMD3...], [EDMD-CONTACT], [GUI-STATS], [GEN3-GUI], [GUI-LIMIT], [GUI-SHOT], [GUI-KEYS], [ENERGY], the scene "
                 "reset and spring lines, 'Invalid wall_x', STOP, Failed, [warn]) of the decision 34 (the rest) runs (the logs are .gitignored), "
                 f"and the runner's log in full; extracted {time.strftime('%Y-%m-%d %H:%M %Z')}\n\n")
        for p in files:
            if os.path.basename(p) != "run.log": continue
            d = os.path.dirname(p)
            ec = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "-"
            lines = [l.rstrip("\n") for l in open(p, errors="ignore") if KEY.match(l)]
            fh.write(f"\n== {os.path.relpath(p, HERE)} (exit {ec})\n" + "\n".join(lines) + "\n")
        fh.write(f"\n== {os.path.relpath(os.path.join(O, 'runner.log'), HERE)} (in full)\n" + open(os.path.join(O, "runner.log")).read())
    small = [p for p in files if os.path.getsize(p) < 1_000_000]
    ign = set()
    if small:
        r = subprocess.run(["git", "check-ignore", "--stdin"], input="\n".join(small), capture_output=True, text=True, cwd=HERE)
        ign = set(r.stdout.split("\n"))
    keep = [p for p in small if p not in ign]
    with open(os.path.join(HERE, "D34S_COMMIT_LIST.txt"), "w") as fh:
        for p in keep: fh.write(os.path.relpath(p, os.path.dirname(os.path.dirname(HERE))) + "\n")
    big = [p for p in files if os.path.getsize(p) >= 1_000_000]
    print(f"files {len(files)} ({sum(os.path.getsize(p) for p in files) / 1e6:.1f} MB); under 1 MB and not ignored: {len(keep)} "
          f"({sum(os.path.getsize(p) for p in keep) / 1e6:.2f} MB); ignored: {len(small) - len(keep)}; 1 MB or more: {len(big)} "
          f"({sum(os.path.getsize(p) for p in big) / 1e6:.1f} MB): " + ", ".join(os.path.relpath(p, HERE) for p in big))


if __name__ == "__main__":
    main()
