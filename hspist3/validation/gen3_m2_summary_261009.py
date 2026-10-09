#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.4 M2 acceptance, sec. 4.7.6): one row per cell of the M2 harness output
(experiments_gen3_m2_261009/m2_audit_output.txt, printed by edmd_core/tests/gen3_m2_harness.c audit), so the notes quote
script-printed tables: reproducibility and the schedule audit; health and the contact audit beside gen2; the ledgers;
the observables (information only).
usage (from hspist3/):  python3 validation/gen3_m2_summary_261009.py [output file]
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HS, "experiments_gen3_m2_261009", "m2_audit_output.txt")


def num(s):
    try:
        return float(s)
    except ValueError:
        return float("nan")


def parse(path):
    cells, c = [], None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"### Cell (\w+): (.*)", line)
        if m:
            c = {"name": m.group(1), "what": m.group(2), "runs": [], "audit": {}, "ledger": {}, "obs": []}
            cells.append(c); continue
        if c is None:
            continue
        m = re.match(r"N = (\d+), box .*T = (\d+) sigma-time", line)
        if m:
            c["N"], c["T"] = int(m.group(1)), int(m.group(2))
        m = re.match(r"\| ([ABC]) gen3[^|]*\| (\w+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\S+) \|", line)
        if m:
            c["runs"].append({"run": m.group(1), "hash": m.group(2), "pair": int(m.group(3)), "wall": int(m.group(4)), "div": int(m.group(5)),
                              "pis": int(m.group(6)), "band": int(m.group(7)), "cross": int(m.group(8)), "stale": int(m.group(9)), "eq": m.group(10)})
        m = re.match(r"audits do not steer \(A = B\): (\w+); same-seed bit identity \(B = C\): (\w+)", line)
        if m:
            c["steer"], c["ident"] = m.group(1), m.group(2)
        m = re.match(r"schedule audit \(gen3, run A\): (\d+) audited states", line)
        if m:
            c["audits"] = int(m.group(1))
        m = re.match(r"\| (pairs|outer walls|crossings|divider faces|pistons) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\S+) \| (\S+) \| (\S+) \(at horizon (\S+)\) \| (\S+) \(at horizon (\S+)\) \|", line)
        if m:
            c["audit"][m.group(1)] = {"cmp": int(m.group(2)), "miss": int(m.group(3)), "extra": int(m.group(4)), "dt": int(m.group(5)), "dtrel": int(m.group(6)),
                                      "early": m.group(8), "maxdt": num(m.group(9)), "dthz": num(m.group(10)), "maxrel": num(m.group(11)), "hz": num(m.group(12))}
        m = re.match(r"bands: missing (\d+), extra (\d+), short (\d+); second live crossing of one disk (\d+); duplicate disagreements (\d+); disks outside their cell (\d+)", line)
        if m:
            c["bands"] = tuple(int(x) for x in m.groups())
        m = re.match(r"\| (gen3 \(A\)|gen2) \| (\d+) \| (\S+) \| (\S+) \| (\S+) \| (\S+) \|$", line)
        if m:
            c["contact_" + ("g3" if m.group(1).startswith("gen3") else "g2")] = [num(x) for x in m.groups()[2:]]
        if line.startswith("[EDMD3-HEALTH]"):
            kv = dict(re.findall(r"(\w+)=(\S+)", line)); c["health"] = kv
        m = re.match(r"run flag \(edmd3_health_clean\): A (\d), B (\d), C (\d)", line)
        if m:
            c["clean"] = m.group(1) + m.group(2) + m.group(3)
        if line.startswith("[EDMD-HEALTH gen2]"):
            kv = dict(re.findall(r"(\w+)=(\S+)", line)); c["gen2"] = kv
            m2 = re.search(r"worst surface gap (\S+) px, states with a gap below -2\.4e-5 px: (\d+)", line)
            m3 = re.search(r"overlap check of (\d+) sampled states", line)
            if m2: c["gen2_worst"], c["gen2_bad"] = num(m2.group(1)), int(m2.group(2))
            if m3: c["gen2_nsamp"] = int(m3.group(1))
        m = re.search(r"v_ref = (\S+) px/unit, K = \S+ -> c_tol = (\S+) px\^2 .*tol_face = (\S+) px", line)
        if m:
            c["tol"] = (num(m.group(1)), num(m.group(2)), num(m.group(3)))
        m = re.match(r"\| (x|y|energy) \| (\S+)(?: \(E - E0\))? \| (\S+)(?: \(W\))? \| (\S+) \| (\S+) \| (\S+) \|", line)
        if m:
            c["ledger"][m.group(1)] = (num(m.group(2)), num(m.group(3)), num(m.group(4)), num(m.group(5)), m.group(6))
        if line.startswith("| engine | Z (pair virial)"):
            c["obs_cols"] = [x.strip() for x in line.strip("|").split("|")]
        m = re.match(r"\| (gen3 \(B\)|gen2) \| (\S+) \| (\S+) \|(.*)", line)
        if m and "obs_cols" in c:
            vals = [x.strip() for x in line.strip().strip("|").split("|")]
            c["obs"].append(dict(zip(c["obs_cols"], vals)))
        m = re.match(r"static method .*F_left = (\S+) kT/px \(Z = F L / \(N kT\) = (\S+),.*F_right = (\S+) kT/px \(Z = (\S+),", line)
        if m:
            c["static"] = tuple(num(x) for x in m.groups())
    return cells


def main():
    cells = parse(OUT)
    print(f"# M2 harness summary, printed by validation/gen3_m2_summary_261009.py from {os.path.relpath(OUT, HS)}\n")
    print("## 1. Reproducibility and the schedule audit (gen3)\n")
    print("| cell | N | T | events A: pair / wall / divider / piston / band / crossings | A = B (audits do not steer) | B = C (same seed) | audited states | "
          "missing (all classes) | extra | deferred earlier than eligible | bands missing / extra / short | second live crossing, duplicate disagreements, disks outside their cell |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in cells:
        a = c["runs"][0]; au = c["audit"]
        miss = sum(v["miss"] for v in au.values()); ext = sum(v["extra"] for v in au.values())
        early = sum(int(v["early"]) for v in au.values() if v["early"] not in ("-",))
        b = c.get("bands", (0,) * 6)
        print(f"| {c['name']} | {c['N']} | {c['T']} | {a['pair']} / {a['wall']} / {a['div']} / {a['pis']} / {a['band']} / {a['cross']} | {c['steer']} | {c['ident']} | "
              f"{c['audits']} | {miss} | {ext} | {early} | {b[0]} / {b[1]} / {b[2]} | {b[3]}, {b[4]}, {b[5]} |")
    print("\n## 2. Amendment e: per class, over ALL matched events: max |dt| [units] (at the horizon of that event); max |dt| / horizon (at that horizon); horizon = max(t_bruteforce, t_heap) - now [units]\n")
    print("| cell | pairs | outer walls | crossings | divider faces | pistons |\n|---|---|---|---|---|---|")
    for c in cells:
        au = c["audit"]; row = []
        for k in ("pairs", "outer walls", "crossings", "divider faces", "pistons"):
            v = au.get(k)
            row.append("-" if not v or v["cmp"] == 0 else f"{v['maxdt']:.3g} (horizon {v['dthz']:.3g}); {v['maxrel']:.3g} (horizon {v['hz']:.3g})")
        print(f"| {c['name']} | " + " | ".join(row) + " |")
    print("\n## 3. Health (gen3 run A; safety nets and validator must be 0), contact audit beside gen2 (max |gap| at executed events, px)\n")
    print("| cell | safety nets and validator (sum) | run flag A, B, C | contact_now / wall / body (c_min px^2) | gen3: pairs / walls / divider / pistons | "
          "gen2: pairs / walls / divider / pistons | gen2 safety nets: overlap_repair, clamp_repair, wall_overdue, forced, past | gen2 sampled states with an overlap (worst gap px) |\n"
          "|---|---|---|---|---|---|---|---|")
    nets = ("overlap_repair", "wall_overdue", "obj_overlap_repair", "past_event", "clamp_repair", "cell_repair", "grid_escape", "stagnation",
            "local_findings", "full_findings", "body_findings")
    for c in cells:
        h = c["health"]; g2 = c["gen2"]
        s = sum(int(h[k]) for k in nets)
        g3c = c["contact_g3"]; g2c = c["contact_g2"]
        print(f"| {c['name']} | {s} | {', '.join(c['clean'])} | {h['contact_now']} / {h['wall_contact_now']} / {h['obj_contact_now']} ({num(h['contact_c_min']):.3g}) | "
              f"{g3c[0]:.3g} / {g3c[1]:.3g} / {g3c[2]:.3g} / {g3c[3]:.3g} | {g2c[0]:.3g} / {g2c[1]:.3g} / {g2c[2]:.3g} / {g2c[3]:.3g} | "
              f"{g2['overlap_repair']}, {g2['clamp_repair']}, {g2['wall_overdue']}, {g2['forced_advance']}, {g2['past_event']} | "
              f"{c.get('gen2_bad', 0)} of {c.get('gen2_nsamp', 0)} ({c.get('gen2_worst', 0.0):.3g}) |")
    print("\n## 4. Ledgers (gen3 run B, end of run): residual and its rounding scale\n")
    print("| cell | x: P - P0 | x: residual / scale | y: P - P0 | y: residual / scale | E - E0 [kT] | W [kT] | energy: residual / scale |\n|---|---|---|---|---|---|---|---|")
    for c in cells:
        L = c["ledger"]
        print(f"| {c['name']} | {L['x'][0]:.6g} | {L['x'][2]:.3g} / {L['x'][3]:.3g} | {L['y'][0]:.6g} | {L['y'][2]:.3g} / {L['y'][3]:.3g} | "
              f"{L['energy'][0]:.6g} | {L['energy'][1]:.6g} | {L['energy'][2]:.3g} / {L['energy'][3]:.3g} |")
    worst = max(((abs(c["ledger"][k][2]) / c["ledger"][k][3], c["name"], k) for c in cells for k in ("x", "y", "energy") if c["ledger"][k][3] > 0))
    print(f"\nlargest |residual| / scale: {worst[0]:.3g} ({worst[1]}, {worst[2]})")
    print("\n## 5. Tolerances in force (end of run A) against the measured contact errors\n")
    print("| cell | v_ref [px/unit] | c_tol [px^2] | 2 d x max pair gap / c_tol | tol_face [px] | max face gap / tol_face |\n|---|---|---|---|---|---|")
    for c in cells:
        v, ct, tf = c["tol"]; g = c["contact_g3"]
        print(f"| {c['name']} | {v:.4g} | {ct:.4g} | {48.0 * g[0] / ct:.3g} | {tf:.4g} | {max(g[1], g[2], g[3]) / tf:.3g} |")
    print("\n## 6. Observables, information only (single trajectories; block SEs ignore slow correlations)\n")
    print("| cell | Z gen3 (SE) | Z gen2 (SE) | divider mean x gen3 / gen2 [px] | divider SD gen3 / gen2 [px] | divider period gen3 / gen2 [sigma-time] "
          "(periods in the record) | static method Z left / right (gen3) | right piston work gen3 / gen2 [kT] |\n|---|---|---|---|---|---|---|---|")
    for c in cells:
        rows = {r.get("engine"): r for r in c["obs"]}
        g3, g2 = rows.get("gen3 (B)"), rows.get("gen2")
        if not g3 or not g2:
            continue
        def col(r, prefix):
            for k, v in r.items():
                if k.startswith(prefix):
                    return v
            return None
        zk = [k for k in g3 if k.startswith("Z (pair virial)")][0]
        dm = f"{col(g3, 'divider mean x')} / {col(g2, 'divider mean x')}" if col(g3, "divider mean x") else "-"
        ds = f"{col(g3, 'SD')} / {col(g2, 'SD')}" if col(g3, "SD") else "-"
        def per(v):   # "P (n)": a peak at the search bound n = 3 resolves no oscillation in the record
            if v is None:
                return "-"
            m = re.match(r"(\S+) \((\d+)\)", v)
            return v if not m or int(m.group(2)) > 3 else "not resolved (peak at the 3-period bound)"
        dp = f"{per(col(g3, 'period'))} / {per(col(g2, 'period'))}" if col(g3, "period") else "-"
        pw = f"{col(g3, 'work of the right piston')} / {col(g2, 'work of the right piston')}" if col(g3, "work of the right piston") else "-"
        st = c.get("static"); sm = f"{st[1]:.4g} / {st[3]:.4g}" if st else "-"
        print(f"| {c['name']} | {g3[zk]} ({g3['SE']}) | {g2[zk]} ({g2['SE']}) | {dm} | {ds} | {dp} | {sm} | {pw} |")

if __name__ == "__main__":
    main()
