#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 48, item 5): the test runs' results, read from their stdout/stderr (tests/<name>.out/.err):
per gen3 build its cell width; per gen3 state (before the APPLY rebuild, and the scene's last) the clean flag, the executed events,
the worst surface gaps, the event hash and t_end; the scene line of the APPLY (N, r, eta); the STOP / STOPPED lines.
usage: results_d48.py <tests folder> <name> ...   (prints a markdown table and one line per run)"""
import re, sys, os

EV = sys.argv[1]
KINDS = ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")


def health(line):
    d = dict(re.findall(r"(\w+)=([^ ]+)", line))
    d["events"] = sum(int(d[k]) for k in KINDS)
    d["worst"] = min(float(d["local_worst"]), float(d["full_worst"]))
    return d


rows = []
for name in sys.argv[2:]:
    out = open(os.path.join(EV, name + ".out"), encoding="utf-8", errors="replace").read()
    err = open(os.path.join(EV, name + ".err"), encoding="utf-8", errors="replace").read()
    cells = re.findall(r"\[EDMD3\] built #\d+ .*?cell ([0-9.]+) px", out)
    init = re.search(r"\[EDMD\] initialized: N=(\d+) .*?R=([0-9.]+)σ -> ([0-9.]+)px", out + err)
    scene = re.search(r"Scene reset queued: .*? N=(\d+) r=([0-9.]+)σ eta=([0-9.]+)", out)
    h = {m.group(1): health(m.group(0)) for m in re.finditer(r"\[EDMD3-HEALTH\] (scene(?:, before the rebuild)?): clean=\d overlap_repair.*", out)}
    nostate = re.search(r"\[EDMD3-HEALTH\] scene: clean=0 no gen3 state", out)
    stop = re.search(r"(STOPPED \(--engine=gen3, in the window[^\n]*|STOP \(--engine=gen3\)[^\n]*)", out + err)
    rc = re.search(r"# exit code (\d+)", open(os.path.join(EV, name + ".cmd")).read())
    rows.append(dict(name=name, cells=cells, init=init, scene=scene, h=h, nostate=nostate, stop=stop, rc=rc.group(1) if rc else "?"))

print("| test | N | radius [σ] (R [px]) | η (APPLY) | cell width of the builds [px] | events before APPLY / after | clean before / after | worst gap [px] | event hash before / after | exit |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    N = r["init"].group(1) if r["init"] else "-"
    rad = f"{float(r['init'].group(2)):.4f} ({float(r['init'].group(3)):.3f})" if r["init"] else "-"
    eta = r["scene"].group(3) if r["scene"] else "- (no APPLY)"
    cells = ", ".join(r["cells"]) if r["cells"] else ("refused: " + r["stop"].group(0).split("edmd_gen3: ")[-1].split(" (sec.")[0] if r["stop"] else "-")
    b, a = r["h"].get("scene, before the rebuild"), r["h"].get("scene")
    ev = f"{b['events']:,} / {a['events']:,}" if b and a else (f"- / {a['events']:,}" if a else ("0 (no gen3 state)" if r["nostate"] else "-"))
    cl = f"{b['clean']} / {a['clean']}" if b and a else (a["clean"] if a else ("0 (no gen3 state)" if r["nostate"] else "-"))
    worst = f"{min(b['worst'], a['worst']):.2e}" if b and a else (f"{a['worst']:.2e}" if a else "-")
    hs = f"{b['hash']} / {a['hash']}" if b and a else "-"
    print(f"| {r['name']} | {N} | {rad} | {eta} | {cells} | {ev} | {cl} | {worst} | {hs} | {r['rc']} |")
