#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 53 item 2; 261012 sec. 4.7.118): THE ORIENTATION PIN of Part A (A5, option b) and its
feasibility check. The plan author's decision: every stencil position of a cell starts with the crystal orientation of the cell's j = 0
arrangement (--gen3-lattice-orientation), so that the stencil measures one state.
  1. Every block B cell, the default lattice start at j = 0, the start only (60 held steps): each compartment's arrangement from the run
     header ([EDMD3-RUN] left=[...] right=[...]) and its crystal orientation as g3_lat_cryst in 00ALLINONE.c computes it (0: close-packed
     lines along x; turned: along y). At j = 0 the two compartments are equal by symmetry; checked.
  2. Every cell at the five positions with --gen3-lattice-orientation = that orientation, the start only: the start must exist (no
     "[EDMD3-M4] INFEASIBLE"), the run must exit 0 with a clean run record, and the [EDMD3-M4C] line must report the requested
     orientation in both compartments. The pinned arrangement is listed against j = 0's (the orientation is pinned; the rows and
     n_s x n_r may still follow the width).
Writes pa_orient_261014.tsv (cell, N, eta, the pinned orientation, j = 0's arrangement), which pa_tasks.py reads, and REPORT.md.
usage (from hspist3/): python3 cluster/gen3_partA_261014/pa_orient.py <binary> <new run folder> [--table PATH]
"""
import concurrent.futures as cf, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pa_tasks import cells, P, POS                         # block B's cells as stage 1 writes them, the registered functions

ENV = dict(os.environ, SDL_VIDEODRIVER="dummy", HD_CONTACT_AUDIT="1")


def cmd(binary, d, N, H, L0, xw, seed, pin):
    a = [binary, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0", "--seed-drift-order=drift-first",
         f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv", "--trace-every=6", f"--particles={N}",
         f"--particles-boxes={N // 2},{N // 2}", "--particle-radius=0.5", f"--l0={L0}", f"--height={H}", "--num-walls=1",
         f"--wall-positions={xw}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
         "--eff-output=wall-ke", "--wall-hold-steps=60", "--steps=6", "--fixed-dt=0.4", "--kbt1", f"--seed={seed}", "--engine=gen3",
         "--gen3-exact-box", "--gen3-seeding=lattice"]
    return a + ([f"--gen3-lattice-orientation={pin}"] if pin else [])


def run(binary, d, N, H, L0, xw, seed, pin):
    os.makedirs(d)                                         # refuses an existing folder
    with open(os.path.join(d, "stdout.log"), "w") as o:
        rc = subprocess.run(cmd(binary, d, N, H, L0, xw, seed, pin), stdout=o, stderr=subprocess.STDOUT, env=ENV, cwd=d).returncode
    log = open(os.path.join(d, "stdout.log"), errors="replace").read()
    arr = []
    m = re.search(r"^\[EDMD3-RUN\].*$", log, re.M)
    if m:
        for side in ("left", "right"):
            a = re.search(side + r"=\[orient (\w+) rows (\w+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) ", m.group(0))
            o, rows, ns, nr, a_s, a_r = a.group(1), a.group(2), int(a.group(3)), int(a.group(4)), float(a.group(5)), float(a.group(6))
            along = 1 if ns < 2 else (0 if nr < 2 else int(a_r < 2.0 * a_s))
            arr.append(dict(orient=o, rows=rows, ns=ns, nr=nr, cryst="0" if along == int(o == "perpendicular") else "turned"))
    m4c = re.search(r"^\[EDMD3-M4C\] lattice orientation requested: (\S+) .*?; left (.+?), right (.+)$", log, re.M)
    clean = bool(re.search(r"^\[EDMD3-HEALTH\] .*: clean=1 ", log, re.M))
    return dict(rc=rc, arr=arr, m4c=m4c.groups() if m4c else None, infeasible="[EDMD3-M4] INFEASIBLE" in log, clean=clean)


def fmt(a):
    return f"{a['orient'][:4]} {a['rows'][:3]} {a['ns']}x{a['nr']} {'turned' if a['cryst'] == 'turned' else 'x'}"


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    binary, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    table = sys.argv[sys.argv.index("--table") + 1] if "--table" in sys.argv else os.path.join(HERE, "pa_orient_261014.tsv")
    if os.path.exists(out):
        sys.exit(f"REFUSED: {out} exists")
    os.makedirs(out)
    C = []
    for stream, (cell, N, elab, L0s, Hs, NS) in enumerate(cells()):
        L0, H = float(L0s), float(Hs)
        dL = 0.5 * P.leff(L0) / math.sqrt(2 * NS * P.cs2(P.eta_of(NS, H, L0)))
        C.append((cell, N, elab, L0s, Hs, {j: (L0s if j == 0 else f"{L0 + j * dL:.17g}") for _, j in POS}, 7790000 + 10 * stream))
    # 1. j = 0, the default start
    with cf.ThreadPoolExecutor(6) as ex:
        base = list(ex.map(lambda c: run(binary, os.path.join(out, f"{c[0]}_j+0_default"), c[1], c[4], c[3], c[5][0], c[6] + 2, None), C))
    pins, problems = {}, []
    for c, b in zip(C, base):
        if b["rc"] != 0 or not b["clean"] or len(b["arr"]) != 2:
            problems.append(f"{c[0]}: the default start at j = 0 did not run (exit {b['rc']})"); continue
        if b["arr"][0]["cryst"] != b["arr"][1]["cryst"]:
            problems.append(f"{c[0]}: at j = 0 the compartments' orientations differ ({fmt(b['arr'][0])} / {fmt(b['arr'][1])})"); continue
        pins[c[0]] = b["arr"][0]["cryst"]
    # 2. the five positions, pinned
    jobs = [(c, j) for c in C if c[0] in pins for _, j in POS]
    with cf.ThreadPoolExecutor(6) as ex:
        res = list(ex.map(lambda cj: run(binary, os.path.join(out, f"{cj[0][0]}_j{cj[1]:+d}_pin{pins[cj[0][0]]}"), cj[0][1], cj[0][4],
                                         cj[0][3], cj[0][5][cj[1]], cj[0][6] + cj[1] + 2, pins[cj[0][0]]), jobs))
    by = {}
    for (c, j), r in zip(jobs, res):
        by.setdefault(c[0], {})[j] = r
    want = {"0": "lines along x", "turned": "lines along y"}
    rows, nfeas, nsame = [], 0, 0
    for c, b in zip(C, base):
        if c[0] not in pins:
            rows.append(f"| {c[0]} | - | - | - | - | - | - | NOT RUN |"); continue
        p = pins[c[0]]; d = by[c[0]]
        ok = all(d[j]["rc"] == 0 and d[j]["clean"] and not d[j]["infeasible"] and d[j]["m4c"] and d[j]["m4c"][1] == want[p]
                 and d[j]["m4c"][2] == want[p] and len(d[j]["arr"]) == 2 and d[j]["arr"][0]["cryst"] == p and d[j]["arr"][1]["cryst"] == p
                 for _, j in POS)
        same = ok and all(d[j]["arr"][0] == b["arr"][0] and d[j]["arr"][1] == b["arr"][1] for _, j in POS)
        nfeas += ok; nsame += same
        if not ok:
            for _, j in POS:
                if d[j]["infeasible"] or d[j]["rc"] != 0:
                    problems.append(f"{c[0]} j = {j:+d}: {'INFEASIBLE' if d[j]['infeasible'] else 'exit ' + str(d[j]['rc'])} with the pin {p}")
        cellsj = " | ".join((f"{fmt(d[j]['arr'][0])} / {fmt(d[j]['arr'][1])}" if len(d[j]["arr"]) == 2 else "-") for _, j in POS)
        rows.append(f"| {c[0]} | {p} | {cellsj} | {'FEASIBLE' if ok else 'NOT FEASIBLE'}, {'as j = 0' if same else 'arrangement differs'} |")
    with open(table, "w") as f:
        f.write("cell\tN\teta\torientation\tarrangement_j0\n")
        for c, b in zip(C, base):
            if c[0] in pins:
                f.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{pins[c[0]]}\t{fmt(b['arr'][0])}\n")
    rep = [f"# pa_orient.py: the orientation pin (decision 53 item 2) with {os.path.relpath(binary)}", "",
           "| cell | pinned orientation (j = 0's) | j = -2: left / right | j = -1 | j = 0 | j = +1 | j = +2 | the pinned start |",
           "|---|---|---|---|---|---|---|---|"] + rows + ["",
           f"{len(pins)} of {len(C)} cells have an orientation at j = 0; {nfeas} of {len(pins)} are feasible at all five positions with it; "
           f"in {nsame} the pinned arrangement equals j = 0's at every position.",
           "Entries: rows' orientation (para/perp), rows (equal/alt), n_s x n_r, crystal orientation (x: close-packed lines along x; "
           "turned: along y).", "", "Problems: " + ("none" if not problems else ""), *[f"- {x}" for x in problems],
           "", f"The table: {os.path.relpath(table)}"]
    open(os.path.join(out, "REPORT.md"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
