#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.14, M3 amendment e): the LOADER GUARD FOR gen3 RUNS, in the style of edmd_acc_guard.

check(path) raises Gen3RunError for data written by a gen3 run (00ALLINONE --engine=gen3) whose run record is dirty or
missing; it returns the path unchanged otherwise. edmd_acc_guard.guard() calls it, so every loader wired to that guard (the 20
paper scripts, tests_20260913's guarded loaders) refuses such data without a change of its own, and gen2 data is untouched:
check() returns before it reads anything else when no governing record shows gen3.

THE RUN RECORD (00ALLINONE.c, g3_run_record; printed into the run log (stdout) at the end of EVERY gen3 trajectory, also when
all counters are 0, so a missing line is itself a finding):
  [EDMD3-HEALTH] <run id>: clean=<1|0> <every EDMD3_Health counter as name=value> validator_every=<K> hash=<16 hex digits>
      t_end=<internal units> engine=gen3 build="<git hash[-dirty]>" target=<build target> run_s=<s> engine_s=<s>
      driver_share=<fraction>
  run ids: speed of sound "L0=<%.1f> M=<int> run=<int> seed=<uint>"; energy transfer "energy-transfer seed=<uint>";
  replay "replay <name>". A run that never advanced prints "clean=0 no gen3 state (the run never advanced)".
  Each gen3 state prints one build line when it is built, "[EDMD3] built #<n> at t = ...": one per trajectory.

WHICH RECORDS GOVERN A DATA FILE: edmd_acc_guard's rule -- the run records and run logs in the file's own directory and in every
ancestor directory below the experiments root (never the root itself, a directory whose name starts with 'experiments');
outside an experiments tree up to the filesystem root. Record kinds are the provenance scan's (klass): run records
(00_COMMAND*, 01_PLOT_COMMANDS.md, command.txt, *.command.txt, run_params.json) and run logs (*.log).

GEN3 DATA: a governing run log holds a gen3 line ([EDMD3-HEALTH] or "[EDMD3] built"), or the nearest directory (own, then
upwards) whose run records name an engine (--engine=gen2|gen3) names gen3.

RULE for gen3 data, refusal on the first of:
  1. a governing run log holds a run record whose clean value is not exactly 1 (0, or unreadable);
  2. a governing run log holds more build lines than run records (a trajectory ended without its record), or no governing
     run log holds any run record at all;
  3. the data file is a per-run file (<name>_L0_<L0>_wallmassfactor_<M>_run<r>.<ext>) and no governing run record has the run
     id L0 = <L0>, M = <M>, run = <r>.
LIMIT: data whose run left no record and no log in its directory chain (stdout not kept) cannot be recognised as gen3; the run
scripts of the gen3 campaigns keep the run log next to the data (sec. 4.7.14).

usage:  import edmd3_health_guard;  edmd3_health_guard.check(path)        (validation/ must be on sys.path)
"""
import functools, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from provenance_edmd_acc_261009 import klass   # the provenance scan's record kinds, as edmd_acc_guard

RECORD = re.compile(r"\[EDMD3-HEALTH\] (.*?): clean=(\S+)")
BUILT = re.compile(r"\[EDMD3\] built #")
ENGINE = re.compile(r"--engine[= ]['\"]?(gen2|gen3)\b")
PER_RUN = re.compile(r"_L0_(\d+)_wallmassfactor_(\d+)_run(\d+)\.[A-Za-z0-9]+$")
SOS_ID = re.compile(r"^L0=([\d.]+) M=(\d+) run=(\d+) seed=\d+$")


class Gen3RunError(RuntimeError):
    """Data written by a gen3 run whose run record says clean=0, or that has no run record."""


def _read(p):
    try:
        with open(p, errors="ignore") as fh:
            return fh.read()
    except OSError:
        return None


@functools.lru_cache(maxsize=None)
def _dir_info(d):
    """For directory d: (engines named by its run records, [(log, n_build, [(id, clean)])] for its run logs with gen3 lines)."""
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return frozenset(), ()
    engines, logs = set(), []
    for n in names:
        p = os.path.join(d, n)
        cls = klass(n)
        if cls not in ("run record", "run log") or not os.path.isfile(p):
            continue
        t = _read(p)
        if not t:
            continue
        if cls == "run record":
            engines.update(m.group(1) for m in ENGINE.finditer(t))
        else:
            if "EDMD3" not in t:
                continue
            recs = tuple((m.group(1), m.group(2)) for m in RECORD.finditer(t))
            nb = len(BUILT.findall(t))
            if recs or nb:
                logs.append((p, nb, recs))
    return frozenset(engines), tuple(logs)


def _chain(p):
    d = p if os.path.isdir(p) else os.path.dirname(p)
    out = []
    while True:
        if os.path.basename(d).startswith("experiments"):
            break                                    # the experiments root: a container, never scanned
        out.append(d)
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return out


def check(path):
    """Return path unless it is gen3 data with a dirty or missing run record (raise Gen3RunError)."""
    p = os.path.abspath(os.fspath(path))
    chain = _chain(p)
    infos = [(d,) + _dir_info(d) for d in chain]
    logs = [lg for _, _, ls in infos for lg in ls]
    named = next((eng for _, eng, _ in infos if eng), frozenset())
    if not logs and "gen3" not in named:
        return path                                  # not gen3 data: untouched
    where = f"refused {path}: gen3 run"
    for lg, nb, recs in logs:
        for rid, clean in recs:
            if clean != "1":
                raise Gen3RunError(f"{where} '{rid}' has clean={clean} ({lg}; 261012 sec. 4.7.14)")
        if nb > len(recs):
            raise Gen3RunError(f"{where} without a run record: {lg} has {nb} gen3 builds and {len(recs)} run records "
                               f"(261012 sec. 4.7.14)")
    allrecs = [rid for _, _, recs in logs for rid, _ in recs]
    if not allrecs:
        raise Gen3RunError(f"{where} without a run record: no governing run log holds an [EDMD3-HEALTH] line (261012 sec. 4.7.14)")
    m = PER_RUN.search(os.path.basename(p))
    if m:
        L0, M, r = int(m.group(1)), int(m.group(2)), int(m.group(3))
        ok = False
        for rid in allrecs:
            s = SOS_ID.match(rid)
            if s and int(round(float(s.group(1)))) == L0 and int(s.group(2)) == M and int(s.group(3)) == r:
                ok = True
                break
        if not ok:
            raise Gen3RunError(f"{where} without a run record: no run record with L0={L0} M={M} run={r} (261012 sec. 4.7.14)")
    return path


def clear_cache():
    _dir_info.cache_clear()
