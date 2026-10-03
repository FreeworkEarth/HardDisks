#!/usr/bin/env python3
"""##CHRIS 2026-10-02: Task E1 -- inventory of every uncommitted modification and untracked path in the main tree, with
sizes, for Chris. Nothing is committed or changed. `owner` marks the paths that CC itself made or is about to commit in this
batch; everything else is listed as Chris's (or unknown) and is left alone.

usage: python3 hspist3/validation/worktree_inventory_261002.py
"""
import os, subprocess, sys
REPO = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip()
CC = {  # paths made by CC (this or earlier batches) or committed by CC in this batch
    "hspist3/00ALLINONE_5190846dirty_20261014": "CC: binary backup (Task B, previous batch)",
    "hspist3/00ALLINONE_pre_boxw_20261014": "CC: binary backup (box-width batch)",
    "hspist3/experiment_validation.c": "CC commits in E2 (build input, never tracked)",
    "hspist3/experiment_validation.h": "CC commits in E2 (build input, never tracked)",
    "hspist3/validation/paper1_figures_20261001.py": "CC (2026-10-01 batch); committed in G1",
    "hspist3/validation/paper2_figures_20261001.py": "CC (2026-10-01 batch); not in this batch",
    "hspist3/validation/worktree_inventory_261002.py": "CC: this script",
}
for f in ("estimator_floor", "massladder_line", "massladder_residuals", "slowmode"):
    for e in ("pdf", "png"):
        CC[f"0000_PLAN_OVERALL/paper1_speedofsound/experiments/final/261001_p1_{f}.{e}"] = "CC (2026-10-01 batch); committed in G1"
for f in ("level1_path", "level4_acf", "level4_bars", "zeta_tcut"):
    for e in ("pdf", "png"):
        CC[f"0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261001_p2_{f}.{e}"] = "CC (2026-10-01 batch); not in this batch"

def size(p):
    full = os.path.join(REPO, p)
    if os.path.islink(full) or os.path.isfile(full): return os.lstat(full).st_size
    if os.path.isdir(full):
        t = 0
        for root, _, files in os.walk(full):
            for f in files:
                try: t += os.lstat(os.path.join(root, f)).st_size
                except OSError: pass
        return t
    return None

def human(b):
    if b is None: return "(deleted)"
    for u in ("B", "KiB", "MiB", "GiB"):
        if b < 1024 or u == "GiB": return f"{b:.0f} {u}" if u == "B" else f"{b:.1f} {u}"
        b /= 1024

out = subprocess.run(["git", "status", "--porcelain=v1", "-z"], capture_output=True, cwd=REPO).stdout.decode().split("\0")
rows = []
for e in out:
    if not e: continue
    code, p = e[:2], e[3:]
    kind = {"??": "untracked", " M": "modified", " D": "deleted", "M ": "staged"}.get(code, code)
    if code == " M" and os.path.exists(os.path.join(REPO, p, ".git")): kind = "submodule: local edits"
    rows.append((p, kind, size(p), CC.get(p.rstrip("/"), "Chris / unknown")))
print("### E1: uncommitted and untracked paths in the main tree (`git status --porcelain`), sizes on disk\n")
print("| # | path | status | size | owner |"); print("|---|---|---|---|---|")
for i, (p, k, s, o) in enumerate(rows, 1):
    print(f"| {i} | `{p}` | {k} | {human(s)} | {o} |")
tot = lambda sel: sum(s for _, _, s, o in rows if s and sel(o))
print(f"\n{len(rows)} paths: {sum(1 for r in rows if r[1]=='modified')} modified, {sum(1 for r in rows if r[1]=='deleted')} deleted, "
      f"{sum(1 for r in rows if r[1]=='untracked')} untracked, {sum(1 for r in rows if r[1].startswith('submodule'))} submodules with local edits; "
      f"size of Chris/unknown paths {human(tot(lambda o: o.startswith('Chris')))}, of CC paths {human(tot(lambda o: o.startswith('CC')))}")
print("\n`git diff --stat` (tracked modifications only):\n")
print("    " + subprocess.run(["git", "diff", "--stat"], capture_output=True, text=True, cwd=REPO).stdout.replace("\n", "\n    "))
