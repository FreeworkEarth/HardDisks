#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 41 items 3 and 17; registration 261012 sec. 4.8.8, amendment A3): THE TIMING TABLE of the watch
rule, on the Mac, from fetched files only (nothing runs, nothing on KOA is read):
  per chunk     the sum of its trajectories' ##RUN seconds (the worker's "##RUN <date> <time> run <r> seed <s> node <host> (<sec> s)" in the
                merged run.log of the chunk's (cell, mass)) / the chunk table's predicted_s_E5, with the node and CPU model of the chunk's log
                that holds "chunk <j> done" (job_logs/<job name>_<array id>_<j>.out of the fetch);
  per CPU model the minimum, median and maximum of that ratio;
  the watch     "If the E5 maximum of a finished group is 2.0 or more, I decide before the next group's block E" (decision 41 item 3): the
                maximum over the chunks on E5 models is printed against 2.0.
The ##RUN seconds count the trajectories that finished (the worker appends a ##RUN section only for a clean run); a chunk whose run
records are fewer than its task lines is marked, never filled in. A chunk run twice (A3) is timed by its successful log.
usage (from hspist3/):
  python3 cluster/gen3_stage1_261012/s1_timing_261012.py <group> [<group> ...] [--root R] [--chunks T] [--tasks L] [--out O] [--all]
      R: the fetched folder (default experiments_gen3_stage1_koa_261012); T, L: the chunk table and task list (default the package's);
      the full per-chunk table goes to <O>/timing_<group>.tsv (--out O, default experiments_gen3_s1_261012/timing; never into the fetched
      folder); --all prints every chunk (default: every chunk up to 60, else the 20 largest ratios)
  python3 cluster/gen3_stage1_261012/s1_timing_261012.py --paste <file> --chunk <j> [--chunks T] [--tasks L]
      chunk j's ratio from any text holding its ##RUN lines (e.g. a block D paste that prints them, or the chunk's run.log)
"""
import argparse, collections, os, re, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
RUN = re.compile(r"^##RUN \S+ \S+ run (\d+) seed (\d+) node (\S+) \((\d+) s\)\s*$")
ON = re.compile(r"^chunk (\d+) on (\S+) \((.*)\): group (\S+), lines (\d+)-(\d+); ")


def rp(p):
    """a path for printing: relative to hspist3/ inside it, absolute outside"""
    r = os.path.relpath(os.path.abspath(p), HS)
    return os.path.abspath(p) if r.startswith("..") else r


def load(chunks, tasks):
    CT = {}
    for l in open(chunks).readlines()[1:]:
        c = l.rstrip("\n").split("\t")
        CT[int(c[0])] = dict(j=int(c[0]), first=int(c[1]), last=int(c[2]), g=c[4], cid=c[5], M=c[6], pred=float(c[7]), cls=int(c[9]))
    TL = [l.split() for l in open(tasks) if l.strip()]
    for c in CT.values():
        rows = TL[c["first"] - 1:c["last"]]
        rels = {r[1] for r in rows}
        if len(rels) != 1: sys.exit(f"STOP: chunk {c['j']}'s lines name {len(rels)} cell folders")
        c["rel"] = rels.pop(); c["runs"] = [(int(r[3]), int(r[4])) for r in rows]          # (s, exact seed) per line
    return CT


def runs_of(text):
    """{run index: (seed, node, seconds)} from ##RUN lines; a run that appears twice keeps its last record (a rerun after A3)"""
    out = {}
    for l in text.splitlines():
        m = RUN.match(l)
        if m: out[int(m.group(1))] = (int(m.group(2)), m.group(3), int(m.group(4)))
    return out


def chunk_logs(d, g):
    """{chunk: (node, CPU model)} from the fetched chunk logs that hold 'chunk <j> done'"""
    out, jl = {}, os.path.join(d, "job_logs")
    for f in sorted(os.listdir(jl)) if os.path.isdir(jl) else []:
        if not f.startswith(f"gen3-s1-{g}_") or not f.endswith(".out"): continue
        t = open(os.path.join(jl, f), errors="replace").read()
        for l in t.splitlines():
            m = ON.match(l)
            if m and re.search(rf"^chunk {m.group(1)} done$", t, re.M): out[int(m.group(1))] = (m.group(2), m.group(3).replace(" (Mac rehearsal)", ""))
    return out


def ratio_row(c, recs):
    have = [(s, recs[s]) for s, seed in c["runs"] if s in recs and recs[s][0] == seed]
    sec = sum(r[2] for _, r in have)
    return dict(j=c["j"], cid=c["cid"], M=c["M"], cls=c["cls"], pred=c["pred"], n=len(c["runs"]), found=len(have), sec=sec,
                ratio=sec / c["pred"] if c["pred"] > 0 else float("nan"))


def group_table(root, g, CT, show_all, outdir):
    d = os.path.join(root, g); base = os.path.join(d, "merged", "stage1", "gen3")
    if not os.path.isdir(base): print(f"## group {g}: no fetched folder {rp(base)}\n"); return None
    cache, rows, logs = {}, [], chunk_logs(d, g)
    for c in sorted((c for c in CT.values() if c["g"] == g), key=lambda c: c["j"]):
        p = os.path.join(d, "merged", c["rel"], "run.log")
        if p not in cache: cache[p] = runs_of(open(p, errors="replace").read()) if os.path.exists(p) else {}
        r = ratio_row(c, cache[p]); r["node"], r["cpu"] = logs.get(c["j"], ("-", "no log with its done line")); rows.append(r)
    os.makedirs(outdir, exist_ok=True); tsv = os.path.join(outdir, f"timing_{g}.tsv")
    with open(tsv, "w") as fh:
        fh.write("chunk\tcell\tM\tclass_min\tpredicted_s_E5\tlines\truns_found\tsum_run_s\tratio\tnode\tcpu\n")
        for r in rows: fh.write(f"{r['j']}\t{r['cid']}\t{r['M']}\t{r['cls']}\t{r['pred']:.1f}\t{r['n']}\t{r['found']}\t{r['sec']}\t{r['ratio']:.4f}\t{r['node']}\t{r['cpu']}\n")
    print(f"## group {g}: {len(rows)} chunks of the chunk table; run records found for {sum(r['found'] for r in rows)} of {sum(r['n'] for r in rows)} "
          f"task lines; the full table: {rp(tsv)}\n")
    shown = rows if (show_all or len(rows) <= 60) else sorted(rows, key=lambda r: -r["ratio"])[:20]
    print(("every chunk:" if shown is rows else "the 20 largest ratios:") + "\n")
    print("| chunk | cell | M | class [min] | predicted s (E5) | runs found / lines | sum of ##RUN s | ratio | node | CPU model |\n|---|---|---|---|---|---|---|---|---|---|")
    for r in shown:
        print(f"| {r['j']} | {r['cid']} | {r['M']} | {r['cls']} | {r['pred']:.1f} | {r['found']} / {r['n']}{' INCOMPLETE' if r['found'] < r['n'] else ''} | "
              f"{r['sec']} | {r['ratio']:.3f} | {r['node']} | {r['cpu']} |")
    by = collections.defaultdict(list)
    for r in rows:
        if r["found"]: by[r["cpu"]].append(r["ratio"])
    print("\n| CPU model | chunks | minimum | median | maximum |\n|---|---|---|---|---|")
    for cpu, v in sorted(by.items()): print(f"| {cpu} | {len(v)} | {min(v):.3f} | {statistics.median(v):.3f} | {max(v):.3f} |")
    e5 = [x for cpu, v in by.items() if "E5" in cpu for x in v]
    if e5:
        mx = max(e5)
        print(f"\nWATCH RULE (decision 41 item 3): the E5 maximum of group {g} is {mx:.3f}: "
              + ("2.0 OR MORE -- the plan author decides before the next group's block E" if mx >= 2.0 else "below 2.0"))
    else:
        print(f"\nWATCH RULE (decision 41 item 3): no chunk of group {g} ran on an E5 model (no E5 maximum)")
    print()
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("groups", nargs="*"); ap.add_argument("--root", default=os.path.join(HS, "experiments_gen3_stage1_koa_261012"))
    ap.add_argument("--chunks", default=os.path.join(HERE, "chunks_stage1_koa_261012.tsv")); ap.add_argument("--tasks", default=os.path.join(HERE, "tasks_stage1_koa_261012.txt"))
    ap.add_argument("--out", default=os.path.join(HS, "experiments_gen3_s1_261012", "timing"))
    ap.add_argument("--all", action="store_true"); ap.add_argument("--paste"); ap.add_argument("--chunk", type=int); a = ap.parse_args()
    CT = load(a.chunks, a.tasks)
    print(f"# The timing table (decision 41, the watch rule), from fetched files: chunk table {rp(a.chunks)}, task list {rp(a.tasks)}\n")
    if a.paste:
        if a.chunk not in CT: sys.exit(f"STOP: chunk {a.chunk} is not in the chunk table")
        c = CT[a.chunk]; r = ratio_row(c, runs_of(open(a.paste, errors="replace").read()))
        print(f"chunk {c['j']} ({c['cid']}, M {c['M']}): run records {r['found']} of {r['n']}; sum of ##RUN seconds {r['sec']}; predicted (E5) "
              f"{c['pred']:.1f} s; ratio {r['ratio']:.3f}; time limit {c['cls']} min = {60 * c['cls'] / c['pred']:.2f} x predicted"
              + ("" if r["found"] == r["n"] else " -- INCOMPLETE: not every task line has its run record"))
        return
    if not a.groups: sys.exit("usage: s1_timing_261012.py <group> [...] | --paste <file> --chunk <j>")
    for g in a.groups: group_table(a.root, g, CT, a.all, a.out)


if __name__ == "__main__":
    main()
