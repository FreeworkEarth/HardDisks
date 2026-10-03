#!/usr/bin/env python3
"""##CHRIS 2026-10-02 (Task L2): every LOCAL module imported by a TRACKED script that is itself NOT tracked -- the scan that
found plot_speed_of_sound_edmd.py on 2026-10-02 (tracked since 4f5fd31), now as a script. Scope: the tracked .py files in
hspist3/, hspist3/validation/, hspist3/cluster/ and hspist3/cluster/confinement_20261013/; a "local module" is a .py file
in one of those directories. Imports anywhere in the file count (top level or inside functions).
usage: python3 hspist3/validation/untracked_imports_scan_261002.py
"""
import ast, os, subprocess, sys
REPO = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip()
DIRS = ["hspist3", "hspist3/validation", "hspist3/cluster", "hspist3/cluster/confinement_20261013"]
tracked = set(subprocess.run(["git", "ls-files", "hspist3"], capture_output=True, text=True, cwd=REPO).stdout.split())
local = {}
for d in DIRS:
    for f in os.listdir(os.path.join(REPO, d)):
        if f.endswith(".py"): local.setdefault(f[:-3], []).append(f"{d}/{f}")
scripts = sorted(f for f in tracked if f.endswith(".py") and os.path.dirname(f) in DIRS)
miss, nparsed = {}, 0
for f in scripts:
    try: tree = ast.parse(open(os.path.join(REPO, f), errors="replace").read())
    except SyntaxError: print(f"(not parsed: {f})"); continue
    nparsed += 1
    for n in ast.walk(tree):
        names = [a.name.split(".")[0] for a in n.names] if isinstance(n, ast.Import) else \
                ([n.module.split(".")[0]] if isinstance(n, ast.ImportFrom) and n.module else [])
        for m in names:
            if m in local and not any(p in tracked for p in local[m]):
                miss.setdefault(m, set()).add(f)
print(f"tracked scripts scanned: {nparsed} (of {len(scripts)}); local modules known: {len(local)}")
print("| untracked module | file(s) | imported by (tracked) |\n|---|---|---|")
for m, users in sorted(miss.items()):
    print(f"| {m} | {', '.join(local[m])} | {len(users)}: {', '.join(sorted(users))} |")
print(f"\nuntracked modules imported by tracked scripts: {len(miss)}" + (" -- none" if not miss else ""))
sys.exit(0 if not miss else 1)
