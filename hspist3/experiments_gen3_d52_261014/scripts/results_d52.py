#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 52, addendum 1 item 3; 261012 sec. 4.7.120): the return-map panel's test results.
  H: per scene, the [GUI-HASH] lines (gen3's event hash, time and executed events at steps 600, 1200 and 1800) of the four runs
     (closed, open, toggled, keys) side by side; EQUAL when all four agree at every step. Also each run's [GEN3-GUI] notes (that the
     panel did open and its keys did act) and its run record's clean flag.
  F: the frames (BMP from the scratch run folder) as PNG into tests/: the whole window and the panel cropped.
usage (from hspist3/): python3 experiments_gen3_d52_261014/scripts/results_d52.py <prefix>   -> tests/<prefix>_RESULTS.md
"""
import glob, os, re, sys
from PIL import Image

E = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = os.path.join(E, "tests")
SP = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad"


def run_facts(name):
    out = open(os.path.join(T, name + ".out"), errors="replace").read()
    hashes = {int(m.group(1)): (m.group(2), m.group(3), m.group(4))
              for m in re.finditer(r"^\[GUI-HASH\] step (\d+): hash=(\w+) t=(\S+) events=(\d+)", out, re.M)}
    notes = re.findall(r"^\[GEN3-GUI\] (return map[^\n]*|I: return map closed)", out, re.M)
    clean = re.findall(r"^\[EDMD3-HEALTH\] scene: clean=(\d)", out, re.M)
    return hashes, notes, clean


def main():
    p = sys.argv[1]
    L = [f"# The return-map panel's tests ({p}), the binary in tests/{p}_H_e030_closed.cmd", ""]
    for sc in ("e030", "e069"):
        L += [f"## H, eta {'0.30' if sc == 'e030' else '0.69'} (N 400): gen3's event hash at exactly steps 600, 1200 and 1800", "",
              "| run | step 600: hash, t, events | step 1200 | step 1800 | panel notes | run record clean |", "|---|---|---|---|---|---|"]
        runs = {}
        for v in ("closed", "open", "toggled", "keys"):
            h, notes, clean = run_facts(f"{p}_H_{sc}_{v}")
            runs[v] = h
            cells = " | ".join((f"{h[s][0]}, {h[s][1]}, {h[s][2]}" if s in h else "-") for s in (600, 1200, 1800))
            L.append(f"| {v} | {cells} | {len(notes)}: {'; '.join(notes)[:150]} | {','.join(clean) or '-'} |")
        eq = all(runs[v].get(s) == runs["closed"].get(s) and s in runs["closed"] for v in runs for s in (600, 1200, 1800))
        L += ["", f"verdict: {'EQUAL' if eq else 'DIFFERENT'} (the four runs at the three steps)", ""]
    L += ["## F, the frames (PNG in tests/)", ""]
    for d in sorted(glob.glob(os.path.join(SP, "d52", "runs", f"{p}_F_*"))):
        b = os.path.join(d, "frame.bmp"); n = os.path.basename(d)
        if not os.path.exists(b):
            L.append(f"- {n}: no frame"); continue
        im = Image.open(b).convert("RGB"); w, h = im.size
        im.save(os.path.join(T, f"{n}.png"))
        im.crop((int(w * 0.6), int(h * 0.66), int(w * 0.8), h)).save(os.path.join(T, f"{n}_panel.png"))
        L.append(f"- {n}: {w} x {h} px window -> tests/{n}.png and tests/{n}_panel.png (the panel)")
    open(os.path.join(T, f"{p}_RESULTS.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
