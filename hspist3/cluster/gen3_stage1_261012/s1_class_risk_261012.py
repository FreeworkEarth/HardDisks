#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 45 item 6, addendum 1 item 3): THE TIME-CLASS RISK TABLE OF B900 AND B1600, on the Mac, from the
chunk table and the measured ratios of the fetched groups (nothing on KOA is read, nothing runs). A measured ratio is a chunk's sum of
##RUN seconds over its predicted_s_E5, with the CPU model of its chunk log: s1_timing_261012.py's timing_<group>.tsv.
  the measured maximum ratio per eta and CPU model, over the measured groups' chunks of that eta (block B: every group has the same 16 etas);
  per target group and time class: its chunks; the smallest class / predicted (with its chunk); the measured maximum ratio at that chunk's
  eta on each CPU model; and the largest share of its class that any chunk of the class reaches with its predicted time times the
  measured maximum at its eta (over the CPU models);
  AT RISK (decision 45 item 6): a chunk whose measured maximum (at its eta, on any CPU model) x its predicted time exceeds 0.8 of its
  class; listed, or "none". A chunk at an eta without a measured ratio is listed as not assessed.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_class_risk_261012.py [--measured B100,B400] [--target B900,B1600]
                           [--timing O] [--chunks T]
"""
import argparse, collections, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
LIMIT = 0.8


def rp(p):
    """a path for printing: relative to hspist3/ inside it, absolute outside (as s1_timing_261012.rp)"""
    r = os.path.relpath(os.path.abspath(p), HS)
    return os.path.abspath(p) if r.startswith("..") else r


def eta_of(cid):
    m = re.search(r"_eta([0-9.]+)", cid); return round(float(m.group(1)), 6) if m else None


def short(cpu):
    m = re.search(r"(E5-\d+ v\d+|Gold \d+\w*|Platinum \d+\w*|Silver \d+\w*|Apple M\d+(?: \w+)?)", cpu)
    return m.group(1) if m else cpu


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--measured", default="B100,B400"); ap.add_argument("--target", default="B900,B1600")
    ap.add_argument("--timing", default=os.path.join(HS, "experiments_gen3_s1_261012", "timing"))
    ap.add_argument("--chunks", default=os.path.join(HERE, "chunks_stage1_koa_261012.tsv")); a = ap.parse_args()
    CT = [l.rstrip("\n").split("\t") for l in open(a.chunks)][1:]          # chunk first last lines group cell M predicted_s_E5 limit_min class_min
    print(f"# INTERIM -- the time-class risk table of {a.target.replace(',', ' and ')} (s1_class_risk_261012.py; at risk: measured maximum x "
          f"predicted > {LIMIT} x class)\n")
    meas, used, per = {}, [], {}
    for g in a.measured.split(","):
        p = os.path.join(a.timing, f"timing_{g}.tsv")
        if not os.path.exists(p): print(f"- {g}: no timing table {rp(p)} (run s1_timing_261012.py {g} after its fetch)"); continue
        rows = [l.rstrip("\n").split("\t") for l in open(p)][1:]
        n = 0
        for r in rows:                                                     # chunk cell M class_min predicted lines runs_found sum ratio node cpu
            if int(r[6]) == 0 or r[10] == "no log with its done line": continue
            k = (eta_of(r[1]), short(r[10])); v = (float(r[8]), g, int(r[0]))
            if k not in meas or v[0] > meas[k][0]: meas[k] = v
            if (g,) + k not in per or v[0] > per[(g,) + k][0]: per[(g,) + k] = v
            n += 1
        used.append(g); print(f"- {g}: {n} chunks timed ({rp(p)})")
    cpus = sorted({k[1] for k in meas})
    print(f"\nthe measured maximum ratio per eta and CPU model, per measured group (ratio (chunk)), and over {', '.join(used) or 'no group'}:\n")
    if cpus:
        cols = [(g, c) for c in cpus for g in used]
        print("| eta | " + " | ".join(f"{c}, {g}" for g, c in cols) + " | " + " | ".join(f"{c}, all" for c in cpus) + " |\n|---|"
              + "---|" * (len(cols) + len(cpus)))
        for e in sorted({k[0] for k in meas}):
            print(f"| {e:.6f} | " + " | ".join(f"{per[(g, e, c)][0]:.3f} ({per[(g, e, c)][2]})" if (g, e, c) in per else "-" for g, c in cols) + " | "
                  + " | ".join(f"**{meas[(e, c)][0]:.3f}**" if (e, c) in meas else "-" for c in cpus) + " |")
    else:
        print("(none yet)")
    risk, unassessed = [], []
    for g in a.target.split(","):
        rows = [dict(j=int(c[0]), cid=c[5], M=c[6], pred=float(c[7]), cls=int(c[9]), eta=eta_of(c[5])) for c in CT if c[4] == g]
        by = collections.defaultdict(list)
        for r in rows: by[r["cls"]].append(r)
        print(f"\n## {g}: {len(rows)} chunks\n")
        print("| class [min] | chunks | smallest class / predicted (chunk: cell, M, predicted s) | 0.8 x that | measured maximum at its eta: "
              + " | ".join(cpus or ["(no CPU model measured)"]) + " | largest share of the class (chunk, CPU model) |\n|---|---|---|---|---|" + "---|" * max(1, len(cpus)))
        for cls in sorted(by):
            rr = by[cls]; mn = min(rr, key=lambda r: 60 * cls / r["pred"]); q = 60 * cls / mn["pred"]
            best = None
            for r in rr:
                for c in cpus:
                    if (r["eta"], c) in meas:
                        sh = meas[(r["eta"], c)][0] * r["pred"] / (60 * cls)
                        if best is None or sh > best[0]: best = (sh, r["j"], c)
            print(f"| {cls} | {len(rr)} | {q:.2f} ({mn['j']}: {mn['cid']}, {mn['M']}, {mn['pred']:.0f}) | {LIMIT * q:.2f} | "
                  + " | ".join(f"{meas[(mn['eta'], c)][0]:.3f}" if (mn["eta"], c) in meas else "-" for c in (cpus or [None])) + " | "
                  + (f"{best[0]:.3f} ({best[1]}, {best[2]})" if best else "-") + " |")
        for r in rows:
            got = [(c, meas[(r["eta"], c)][0]) for c in cpus if (r["eta"], c) in meas]
            if not got: unassessed.append((g, r)); continue
            for c, x in got:
                if x * r["pred"] > LIMIT * 60 * r["cls"]: risk.append((g, r, c, x))
    print(f"\n## AT RISK (measured maximum at its eta x predicted > {LIMIT} x class, on any CPU model)\n")
    if risk:
        print("| group | chunk | cell | M | predicted s | class [min] | CPU model | measured maximum | projected s | share of the class |\n|---|---|---|---|---|---|---|---|---|---|")
        for g, r, c, x in risk: print(f"| {g} | {r['j']} | {r['cid']} | {r['M']} | {r['pred']:.0f} | {r['cls']} | {c} | {x:.3f} | {x * r['pred']:.0f} | {x * r['pred'] / (60 * r['cls']):.3f} |")
    elif cpus and not unassessed:
        print("none")
    else:
        print("none among the chunks assessed")
    if unassessed:
        print(f"\nnot assessed (no measured ratio at their eta): {len(unassessed)} chunks, etas "
              + ", ".join(f"{e:.6f}" for e in sorted({r['eta'] for _, r in unassessed})))


if __name__ == "__main__":
    main()
