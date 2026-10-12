#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 38 item 1; 261012 sec. 4.7.75): REPORT S-K, a READ-ONLY audit, no code change.
The defect (b) of sec. 4.7.74 A.5: with --spring-wall but without --spring-k / --spring-k-sigma, main() warns "no spring force applied",
yet the struct default k = 1000 /px^2 (= 576000 kT/sigma^2) stays and acts (00ALLINONE.c: the --spring-wall branch of main sets no
constant; initialize_energy_measurement sets it only when both flags are given; the spring acts only in the spring output mode, i.e.
not with --eff-output=wall-ke, where spring_enabled is 0 and the engine gets no divider_k).
For every energy-transfer data folder that is kept (both worktrees, every file on disk, tracked or not), from the recorded commands:
  - every summary CSV with a "command" column: one row per run, with the run's full command line and the constant it ran with
    (column spring_k = energy_measurement.spring_constant at the end of the run, per px^2; the summary writer, 00ALLINONE.c);
  - every 00_COMMAND.md (the command of the folder's invocations) and run_params.json;
  - every task list or runner script under hspist3/ that contains "--spring-wall" (commands launched elsewhere, e.g. on KOA).
Counts per source and per folder: energy-transfer runs, runs with --spring-wall, of those without an explicit spring-constant flag;
every such run listed with the constant it ran with (from its summary row where one exists, else from the code's rule) and its output
mode. Paper 2 = hspist3/experiments_energy_transfer/level*; everything else is listed separately.
usage: python3 hspist3/validation/audit_spring_k_261012.py   (from the main worktree's root)
"""
import csv, glob, json, os, re, sys
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)
ROOTS = ["/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks",
         "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3"]
PPS2 = 24.0 * 24.0
K_DEFAULT_PX = 1000.0
ET = re.compile(r"--experiment[= ]energy_transfer|--energy-transfer-summary|--energy-transfer-trace|--spring-wall|--eff-output|--energy-measurement")


def flags(cmd):
    sw = re.search(r"--spring-wall(?:=|\s+)(\S+)", cmd)
    k = re.search(r"--spring-k(?:-sigma)?(?:=|\s+)(\S+)", cmd)
    ks = re.search(r"--spring-k-sigma(?:=|\s+)(\S+)", cmd)
    eff = re.search(r"--eff-output(?:=|\s+)(\S+)", cmd)
    return dict(spring_wall=sw.group(1) if sw else None, k_flag=(k.group(0) if k else None), k_sigma=ks.group(1) if ks else None,
                eff=(eff.group(1) if eff else "spring"))


def area(path):
    rel = path
    for r in ROOTS:
        if path.startswith(r + "/"): rel = os.path.basename(r) + "/" + path[len(r) + 1:]
    m = re.search(r"experiments_energy_transfer/([^/]+)", rel)
    paper2 = bool(m and m.group(1).startswith("level"))
    return rel, (m.group(1) if m else rel.split("/")[2] if rel.count("/") >= 2 else rel), paper2


def walk(names):
    for r in ROOTS:
        for dp, dn, fn in os.walk(r):
            if "/.git" in dp or dp.endswith(".git"): dn[:] = []; continue
            dn[:] = [d for d in dn if d != ".git"]
            for f in fn:
                if names(f): yield os.path.join(dp, f)


def main():
    runs = []            # (source, path, row id, command, recorded spring_k or None)
    for p in walk(lambda f: f.startswith("summary") and f.endswith(".csv") or f.startswith("energy_transfer_runs") and ".csv" in f):
        try:
            with open(p, newline="", errors="ignore") as fh:
                rd = csv.DictReader(fh)
                if not rd.fieldnames or "command" not in rd.fieldnames: continue
                for i, row in enumerate(rd):
                    cmd = row.get("command") or ""
                    if not ET.search(cmd): continue
                    runs.append(("summary", p, f"row {i + 1} seed {row.get('seed', '?')}", cmd, row.get("spring_k")))
        except Exception as e:
            print(f"(unreadable {p}: {e})")
    for p in walk(lambda f: f == "00_COMMAND.md"):
        t = open(p, errors="ignore").read()
        for i, cmd in enumerate(re.findall(r"```sh\n(.*?)\n```", t, re.S)):
            if ET.search(cmd): runs.append(("00_COMMAND.md", p, f"command {i + 1}", cmd.strip(), None))
    for p in walk(lambda f: f == "run_params.json"):
        try: d = json.load(open(p))
        except Exception: continue
        s = json.dumps(d)
        if "spring" in s or "energy_transfer" in s:
            cmd = d.get("command") or d.get("argv") or ""
            runs.append(("run_params.json", p, "-", cmd if isinstance(cmd, str) else " ".join(map(str, cmd)), None))
    tl = []
    for p in walk(lambda f: f.endswith((".txt", ".tsv", ".sh", ".sbatch", ".py", ".md")) and "COMMAND" not in f):
        if "/hspist3/" not in p or "/validation/audit_spring_k" in p: continue
        try: t = open(p, errors="ignore").read()
        except Exception: continue
        if "--spring-wall" not in t: continue
        for ln in t.splitlines():   # a run's command line: the flag with the run's arguments (not a --help check or an echo)
            if "--spring-wall" in ln and ("--particles" in ln or "--experiment" in ln): tl.append((p, ln.strip()))

    print("# Report S-K: the spring-wall audit (READ-ONLY; printed by hspist3/validation/audit_spring_k_261012.py)\n")
    print(f"sources: {len(runs)} recorded energy-transfer commands (summary rows, 00_COMMAND.md, run_params.json) in the two worktrees; "
          f"{len(tl)} command lines with --spring-wall in task lists and runner scripts\n")
    by = defaultdict(Counter); flagged = []; kcheck = Counter(); kbad = []
    for src, p, rid, cmd, k in runs:
        rel, folder, paper2 = area(p); f = flags(cmd)
        key = (("Paper 2 " if paper2 else "other ") + folder, src)
        by[key]["runs"] += 1
        if f["spring_wall"] is not None:
            by[key]["spring_wall"] += 1
            if not f["k_flag"]:
                by[key]["no_k"] += 1
                flagged.append((src, rel, rid, f, k, cmd))
            elif src == "summary" and k not in (None, ""):   # the recorded constant against the flag's value
                kf = re.search(r"--spring-k-sigma(?:=|\s+)(\S+)", cmd); kp = re.search(r"--spring-k(?:=|\s+)(\S+)", cmd)
                want = float(kf.group(1)) / PPS2 if kf else float(kp.group(1))
                ok_ = abs(float(k) - want) <= 1e-6 * max(1.0, abs(want))
                kcheck["equal" if ok_ else "different"] += 1
                if not ok_: kbad.append((rel, rid, k, want))
    print("## Per folder and source: energy-transfer runs, with --spring-wall, and of those without --spring-k / --spring-k-sigma\n")
    print("| folder | source | runs | with --spring-wall | of those without a spring-constant flag |\n|---|---|---|---|---|")
    for (folder, src), c in sorted(by.items()):
        print(f"| {folder} | {src} | {c['runs']} | {c['spring_wall']} | {c['no_k']} |")
    for paper2 in (True, False):
        for src in ("summary", "00_COMMAND.md", "run_params.json"):
            cs = [c for (f_, s_), c in by.items() if s_ == src and f_.startswith("Paper 2") == paper2]
            if cs:
                print(f"| **{'Paper 2' if paper2 else 'other'}, all** | {src} | {sum(c['runs'] for c in cs)} | {sum(c['spring_wall'] for c in cs)} | {sum(c['no_k'] for c in cs)} |")
    print(f"\n## The runs with --spring-wall and a spring-constant flag: the constant in their summary row (spring_k, per px^2, the value the run "
          f"used) against the flag's value (--spring-k-sigma / 24^2, or --spring-k): equal {kcheck['equal']}, different {kcheck['different']}"
          + "".join(f"\n- {r} {i}: recorded {k}, flag {w}" for r, i, k, w in kbad[:20]))
    print("\n## Every recorded run with --spring-wall and without a spring-constant flag\n")
    if not flagged:
        print("none")
    else:
        print("| source | file | run | --spring-wall | output mode | the constant it ran with |\n|---|---|---|---|---|---|")
        for src, rel, rid, f, k, cmd in flagged:
            if k not in (None, ""):
                kk = float(k); eff_k = f"{kk:g} /px^2 = {kk * PPS2:g} kT/sigma^2 (its summary row)"
            else:
                eff_k = f"{K_DEFAULT_PX:g} /px^2 = {K_DEFAULT_PX * PPS2:g} kT/sigma^2 (the code's rule; no summary row)"
            if f["eff"] != "spring": eff_k += "; NOT applied: wall-ke output mode, no spring acts"
            print(f"| {src} | {rel} | {rid} | {f['spring_wall']} | {f['eff']} | {eff_k} |")
    print("\n## Task lists and runner scripts: command lines with --spring-wall\n")
    tc = Counter(); tno = []
    for p, ln in tl:
        f = flags(ln); tc["lines"] += 1
        if not f["k_flag"]: tno.append((p, ln))
    print(f"lines with --spring-wall: {tc['lines']}; of those without a spring-constant flag: {len(tno)}")
    for p, ln in tno: print(f"- {area(p)[0]}: `{ln[:300]}`")
    n2 = sum(1 for x in flagged if area(x[1])[2])
    wt = Counter((os.path.basename(next(r for r in ROOTS if p.startswith(r + "/"))), src, "--spring-wall" in cmd)
                 for src, p, rid, cmd, k in runs if area(p)[2])
    uniq = {(cmd, rid.split(" seed ")[-1]) for src, p, rid, cmd, k in runs if area(p)[2] and src == "summary"}
    print("\nPaper 2 records by worktree (the engine-gen3 worktree holds tracked copies of some level folders, so the counts above include "
          "those copies): " + "; ".join(f"{w} {s_}: {c} ({'with' if sw else 'without'} --spring-wall)" for (w, s_, sw), c in sorted(wt.items()))
          + f"; distinct Paper 2 summary runs (command and seed): {len(uniq)}, of them with --spring-wall: "
          + f"{sum(1 for c, sd in uniq if '--spring-wall' in c)}")
    print(f"\nRECORDED RUNS WITH --spring-wall AND NO SPRING-CONSTANT FLAG: Paper 2 {n2}, other {len(flagged) - n2}; "
          f"task-list and script lines without a spring-constant flag: {len(tno)}")


if __name__ == "__main__":
    main()
