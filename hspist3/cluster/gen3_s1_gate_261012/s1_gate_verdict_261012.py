#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; registration 261012 sec. 4.8.10, G6 (c) and (d)): THE STAGE 1 GATE VERDICT, read on the
Mac from the fetch (fetch_s1_gate_koa.sh), nothing runs:
  G6 (c)  with the campaign binary: E0 -- "audit vs plain" and "plain vs ref" IDENTICAL for ctrl_min, ctrl_leg, afix, afix_leg, with the
          default build and with --engine=gen2 (8 of 8; e0_*/report.txt); the M1/M2 harness hashes -- the SHA-256 of the three fetched
          outputs equal to the a643bc7 gate's KOA outputs (koa_harness_ref_a643bc7.sha256; 3 of 3); cross-node determinism -- the
          compare job's lines (every compared file IDENTICAL across nodes, and at least two nodes: else NOT TESTED);
  G6 (d)  the 12 gen3 lines of Test G chunk 1, the gate binary against the campaign binary: event hashes and counts equal in 12 of 12
          (s1_gate_testg_compare.py on the fetched run.logs);
  and the build-time record, from which build_record_f72d764.txt is written (its 00ALLINONE line and version line).
usage (from hspist3/): python3 cluster/gen3_s1_gate_261012/s1_gate_verdict_261012.py [experiments_gen3_s1_gate_koa_261012] [--e0-cases a,b]
  --e0-cases  the E0 cases the verdict requires (default: all four, as the gate job runs them; the Mac rehearsal runs two, a listed shim)
"""
import glob, hashlib, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))


def main():
    av = sys.argv[1:]; cases = ("ctrl_min", "ctrl_leg", "afix", "afix_leg")
    if "--e0-cases" in av: cases = tuple(av[av.index("--e0-cases") + 1].split(",")); i = av.index("--e0-cases"); av = av[:i] + av[i + 2:]
    D = os.path.abspath(av[0] if av else os.path.join(HS, "experiments_gen3_s1_gate_koa_261012"))
    G, X = os.path.join(D, "s1_gate"), os.path.join(D, "s1_xnode")
    print(f"# The stage 1 gate (G6 (c) and (d)), read from {os.path.relpath(D, HS)}\n")
    rec = open(os.path.join(G, "BUILD_RECORD_S1.txt")).read() if os.path.exists(os.path.join(G, "BUILD_RECORD_S1.txt")) else ""
    print("## The build-time record (BUILD_RECORD_S1.txt)\n\n```\n" + (rec.strip() or "MISSING") + "\n```\n")
    log = open(os.path.join(G, "job_log.out"), errors="replace").read() if os.path.exists(os.path.join(G, "job_log.out")) else ""
    node = re.search(r"^s1 gate on (.*?): out ", log, re.M)
    print(f"gate job node: {node.group(1) if node else '(no job log)'}\n")
    e0 = []
    for tag in ("default", "gen2flag"):
        p = os.path.join(G, f"e0_{tag}", "report.txt"); t = open(p, errors="replace").read() if os.path.exists(p) else ""
        for case in cases:
            m = re.search(rf"^\| {case} \| git.*$", t, re.M)
            cells = [c.strip() for c in m.group(0).split("|")] if m else []
            ok = bool(cells) and cells[-3] == "IDENTICAL" and cells[-2] == "IDENTICAL"
            e0.append((tag, case, ok, f"audit vs plain {cells[-3]}, plain vs ref {cells[-2]}" if cells else "no row"))
    print("## G6 (c): E0 (rule 4) with the campaign binary\n\n| build | case | result | detail |\n|---|---|---|---|")
    for tag, case, ok, det in e0: print(f"| {tag} | {case} | {'IDENTICAL' if ok else '**NOT IDENTICAL**'} | {det} |")
    ref = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, "koa_harness_ref_a643bc7.sha256")) if l.strip() and not l.startswith("#")}
    print("\n## G6 (c): the M1/M2 harness outputs against the a643bc7 gate's KOA outputs\n\n| output | campaign SHA-256 | a643bc7 gate | equal |\n|---|---|---|---|")
    hs = []
    for f, want in ref.items():
        p = os.path.join(G, "harness", f); k = hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else "missing"
        hs.append(k == want); print(f"| {f} | {k} | {want} | {'yes' if k == want else '**no**'} |")
    print("\n## G6 (c): cross-node determinism with the campaign binary (the compare job's log)\n")
    cl = sorted(glob.glob(os.path.join(X, "job_logs", "gen3-s1-xnode-compare_*.out")))
    ct = open(cl[-1], errors="replace").read() if cl else ""
    print("```\n" + (ct.strip() or "NO COMPARE LOG") + "\n```")
    xn_ok = bool(ct) and "NOT TESTED" not in ct and "DIFFERENT" not in ct and "IDENTICAL" in ct
    print("\n## G6 (d): Test G chunk 1, the 12 gen3 lines, the gate binary against the campaign binary\n")
    r = subprocess.run([sys.executable, os.path.join(HERE, "s1_gate_testg_compare.py"), os.path.join(G, "testg_gate"), os.path.join(G, "testg_campaign"),
                        os.path.join(G, "testg_chunk1_gen3_lines.txt")], capture_output=True, text=True)
    print(r.stdout.strip())
    m = re.search(r"event hashes equal in (\d+) of (\d+)", r.stdout); d_ok = bool(m) and m.group(1) == m.group(2) == "12"
    c_ok = all(x[2] for x in e0) and len(e0) == 2 * len(cases) and all(hs) and len(hs) == 3 and xn_ok
    print(f"\nG6 (c): {'PASS' if c_ok else 'FAIL'} (E0 {sum(x[2] for x in e0)} of {2 * len(cases)} IDENTICAL{'' if len(cases) == 4 else ' -- ONLY ' + ','.join(cases) + ' (rehearsal)'}; harness {sum(hs)} of 3 equal; cross-node "
          f"{'IDENTICAL across nodes' if xn_ok else 'NOT shown (see the compare log)'})")
    print(f"G6 (d): {'PASS' if d_ok else 'FAIL'} ({m.group(0) if m else 'no comparison'})")


if __name__ == "__main__":
    main()
