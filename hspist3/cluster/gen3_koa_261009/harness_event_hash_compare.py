#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 20, task 6; 261012 sec. 4.7.47): the EVENT HASHES inside KOA's M1/M2 harness outputs (the
gate job's part (b), job 15040861: m1_audit_output.txt, m2_audit_quick_output.txt, m2_audit_output.txt) against the Mac's committed
outputs (engine-gen3: experiments_gen3_m1_261008/m1_audit_output.txt, experiments_gen3_m2_261009/m2_audit_{quick_,}output.txt). The
files' SHA-256 cannot agree (each file names its compiler in its "build:" line); the event hashes are the comparison that means
something. Per file and cell: the rows A (audits on), B (audits off), C (audits off, again) -- event hash and event counts -- KOA
against the Mac, and A = B = C within each machine (the audits do not steer; a repeat is identical). Each cell is marked by how its
initial velocities are made: "libm" if the harness draws them with log/cos/sin (gauss_velocities: Box-Muller), which macOS's libm and
glibc may round differently in the last bit, else "exact" (equal speeds on a lattice, the cradles). Read-only.
usage (from hspist3/ of the engine-gen3 worktree, after the fetch):
  python3 cluster/gen3_koa_261009/harness_event_hash_compare.py [--koa experiments_gen3_gate_koa_261009/gate/harness]
  python3 cluster/gen3_koa_261009/harness_event_hash_compare.py --rows <a text file of KOA's "audits on" rows, in m2 quick, m2, m1 order>
          (the preliminary check: the 32 rows of Chris's read-only grep of 2026-10-10, before the fetch)
"""
import argparse, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
MAC = {"m1_audit_output.txt": os.path.join(HS, "experiments_gen3_m1_261008", "m1_audit_output.txt"),
       "m2_audit_quick_output.txt": os.path.join(HS, "experiments_gen3_m2_261009", "m2_audit_quick_output.txt"),
       "m2_audit_output.txt": os.path.join(HS, "experiments_gen3_m2_261009", "m2_audit_output.txt")}
EXACT = ("cradle", "tie")          # cells whose initial velocities are not drawn with libm (M1 make_tie; M2's three cradles)
ROW = re.compile(r"^\| ([ABC]) gen3, audits (on|off|off, again) \| ([0-9a-f]{16}) \| (.*)$")


def cells(path):
    out, name = [], None
    for l in open(path, errors="ignore"):
        m = re.match(r"^### (?:Cell|Run|Case) (\S+?):?\s", l)
        if m: name = m.group(1)
        r = ROW.match(l.rstrip("\n"))
        if r: out.append((name, r.group(1), r.group(3), r.group(4).rstrip(" |")))
    return out


def kind(name): return "exact" if name and any(name.startswith(e) for e in EXACT) else "libm"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--koa", default=os.path.join(HS, "experiments_gen3_gate_koa_261009", "gate", "harness"))
    ap.add_argument("--rows"); a = ap.parse_args()
    print("# The event hashes in KOA's M1/M2 harness outputs against the Mac's (261012 sec. 4.7.47; decision 20, task 6)\n")
    tot, same_n = 0, 0
    if a.rows:
        koa = [l.split("|")[2].strip() for l in open(a.rows) if "audits on" in l]
        mac = []
        for f in ("m2_audit_quick_output.txt", "m2_audit_output.txt", "m1_audit_output.txt"):
            mac += [(f, n, h) for n, v, h, _ in cells(MAC[f]) if v == "A"]
        print("PRELIMINARY (Chris's grep of the A rows, lines cut at 90 characters; the hash column complete):\n")
        print("| # | file | cell | initial velocities | KOA | Mac | equal |\n|---|---|---|---|---|---|---|")
        for i, (k, (f, n, h)) in enumerate(zip(koa, mac)):
            tot += 1; same_n += k == h
            print(f"| {i + 1} | {f} | {n} | {kind(n)} | {k} | {h} | {'yes' if k == h else 'no'} |")
        if len(koa) != len(mac): print(f"\n**the counts differ: KOA {len(koa)} rows, Mac {len(mac)}**")
    else:
        for f, mp in MAC.items():
            kp = os.path.join(a.koa, f)
            if not os.path.exists(kp): print(f"## {f}: **not fetched** ({kp})\n"); continue
            K, M = cells(kp), cells(mp)
            print(f"## {f}\n\n| cell | row | initial velocities | KOA hash | Mac hash | equal | KOA counts = Mac counts |\n|---|---|---|---|---|---|---|")
            for (kn, kv, kh, kc), (mn, mv, mh, mc) in zip(K, M):
                tot += 1; same_n += kh == mh
                print(f"| {kn} | {kv} | {kind(kn)} | {kh} | {mh} | {'yes' if kh == mh else 'no'} | {'yes' if kc == mc else 'no'} |")
            for lab, X in (("KOA", K), ("Mac", M)):
                by = {}
                for n, v, h, c in X: by.setdefault(n, set()).add((h, c))
                print(f"\n{lab}: A = B = C (hash and counts) in {sum(len(s) == 1 for s in by.values())} of {len(by)} cells")
            if len(K) != len(M): print(f"**the row counts differ: KOA {len(K)}, Mac {len(M)}**")
            print()
    print(f"\nequal event hashes: {same_n} of {tot}")


if __name__ == "__main__":
    main()
