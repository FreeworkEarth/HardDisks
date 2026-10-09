#!/usr/bin/env python3
"""##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the LOADER PROVENANCE GUARD.

Every loader that feeds a paper figure or table calls guard(path) on the data it reads. guard() raises AcceleratedRunError for
data whose run records say that the accelerated EDMD backend (edmd_core/edmd_accelerated.c, selected by --edmd-acc) wrote them:
that backend misses collisions (261012 sec. 4.7; known since 2026-08-26). It returns the path unchanged otherwise, so a call
can wrap the path a loader opens (pd.read_csv(guard(p))) and the loader's numbers cannot move.

PARSING: the provenance scan's own (validation/provenance_edmd_acc_261009.py): FLAG finds every '--edmd-acc' and its value,
BACKEND the binary's 'EDMD backend: accelerated|default' log line, klass() sorts files into run records, run logs and
summaries. RULE, stricter than the scan's classification: only the explicit default spellings 0, false, no and off pass. Every
other value refuses: =1 or any other non-zero digit, =yes, =true, =on, the bare flag, and any value the parser cannot read.

WHICH RECORDS GOVERN A DATA FILE: those in its own directory and in every ancestor directory below the experiments root --
run records (00_COMMAND*, 01_PLOT_COMMANDS.md, command.txt, *.command.txt, run_params.json), run logs (*.log) and CSV summaries
with a 'command' column (a campaign's summary CSV can be its only record, e.g. mass_sweep_eta02_N600). Never the experiments
root itself (a directory whose name starts with 'experiments'): it collects many campaigns -- e.g.
experiments_energy_transfer/energy_transfer_runs_2walls.csv logs unrelated runs -- and must not decide for any of them. A path
outside any experiments tree is checked up to the filesystem root. Each directory is scanned once per process (cache;
clear_cache() for tests). LIMIT: a copy of data kept outside its run's directory chain (e.g. an analysis/ copy of traces
whose 00_COMMAND.md sits in a sibling raw_simulations/ folder) cannot be traced; the campaign-level provenance scan
(provenance_edmd_acc_261009.py, 261012 sec. 4.7.1) covers that case, and no paper script reads such copies.

usage:  import edmd_acc_guard;  df = pd.read_csv(edmd_acc_guard.guard(path))      (validation/ must be on sys.path)
        python3 validation/edmd_acc_guard.py --list    prints every guard() call in hspist3/ (the wired loaders)
"""
import functools, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from provenance_edmd_acc_261009 import FLAG, BACKEND, klass   # the provenance scan's own parsing

DEFAULT_VALUES = frozenset({"0", "false", "no", "off"})


class AcceleratedRunError(RuntimeError):
    """Data written by a run of the accelerated EDMD backend (or a run record that does not say the default backend)."""


def refusing_value(v):
    """None if the flag value selects the default backend, else a reason. v is FLAG's captured value (None = bare flag)."""
    if v is None or v.strip() == "":
        return "bare --edmd-acc"
    s = v.strip("'\",;)`").lower()
    return None if s in DEFAULT_VALUES else f"--edmd-acc value '{v}'"


def text_refusal(text):
    """The first refusing occurrence in a text, or None."""
    if "edmd-acc" in text:
        for m in FLAG.finditer(text):
            r = refusing_value(m.group(1) if m.group(1) is not None else m.group(2))
            if r:
                return r
    if "EDMD backend" in text:
        for m in BACKEND.finditer(text):
            if m.group(1) == "accelerated":
                return "log line 'EDMD backend: accelerated'"
    return None


def _file_text(p, cls):
    try:
        with open(p, errors="ignore") as fh:
            if cls == "summary csv":
                head = fh.readline()
                if "command" not in head.lower():
                    return None
                return head + fh.read()
            return fh.read()
    except OSError:
        return None


@functools.lru_cache(maxsize=None)
def _dir_refusal(d):
    """The first refusing governing record directly in directory d, as 'path: reason', or None."""
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return None
    for n in names:
        cls = klass(n)
        if cls not in ("run record", "run log", "summary csv"):
            continue
        p = os.path.join(d, n)
        if not os.path.isfile(p):
            continue
        t = _file_text(p, cls)
        if t:
            r = text_refusal(t)
            if r:
                return f"{p}: {r}"
    return None


def guard(path):
    """Return path if no governing record says the accelerated backend wrote it; raise AcceleratedRunError otherwise."""
    p = os.path.abspath(os.fspath(path))
    d = p if os.path.isdir(p) else os.path.dirname(p)
    while True:
        if os.path.basename(d).startswith("experiments"):
            break                                    # the experiments root: a container, never scanned
        r = _dir_refusal(d)
        if r:
            raise AcceleratedRunError(f"refused {path}: written by the accelerated EDMD backend ({r}; 261012 sec. 4.7)")
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return path


def clear_cache():
    _dir_refusal.cache_clear()


INDIRECT = re.compile(r"\bT\.(_load|cell_runs|a1_leaf_table|health_of|assert_wall_thickness)\(")


def wired_calls(root, pattern=re.compile(r"\bedmd_acc_guard\.guard\(")):
    """Every guard() call (or, with pattern=INDIRECT, every call of a guarded loader of tests_20260913) in the Python
    sources under root, outside the data trees and the tests: (relpath, line number, line)."""
    out = []
    for dp, dn, fn in os.walk(root):
        dn[:] = [x for x in dn if not x.startswith("experiments") and x not in (".git", "__pycache__", "kissfft")]
        for f in sorted(fn):
            if not f.endswith(".py") or f == os.path.basename(__file__) or f.startswith("test_"):
                continue
            q = os.path.join(dp, f)
            for i, l in enumerate(open(q, errors="ignore"), 1):
                if pattern.search(l) and not l.lstrip().startswith("#"):
                    out.append((os.path.relpath(q, root), i, l.strip()))
    return sorted(out)


if __name__ == "__main__":
    if "--list" in sys.argv:
        HS = os.path.dirname(HERE)
        rows = wired_calls(HS)
        print("### Direct guard() calls (the loaders)\n\n| file | line | call |\n|---|---|---|")
        for f, i, l in rows:
            print(f"| {f} | {i} | `{l[:120]}` |")
        print(f"\n{len(rows)} guard() calls in {len({r[0] for r in rows})} files")
        ind = [r for r in wired_calls(HS, INDIRECT) if r[0] != "validation/tests_20260913.py"]
        print("\n### Loaders guarded through tests_20260913 (calls of T._load, T.cell_runs, T.a1_leaf_table, T.health_of, "
              "T.assert_wall_thickness, each of which calls guard())\n\n| file | line | call |\n|---|---|---|")
        for f, i, l in ind:
            print(f"| {f} | {i} | `{l[:120]}` |")
        print(f"\n{len(ind)} calls in {len({r[0] for r in ind})} files")
