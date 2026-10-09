# ##CHRIS 2026-10-09 (261012 sec. 4.7.10, decision 5): which data files feed the papers. Python imports usercustomize at start-up
# in every process when its folder is on PYTHONPATH, also in multiprocessing workers (spawn) and subprocesses. Only when GUARD_LOG
# is set, this wraps edmd_acc_guard.guard() so that every path a loader hands to the guard is appended to $GUARD_LOG; the real
# guard is still called and its result returned, so the scripts behave as in a normal run. The module is imported from
# $GUARD_VALIDATION_DIR (the absolute validation/ path the scripts themselves put on sys.path), so later imports get this copy.
import os, sys
_log = os.environ.get("GUARD_LOG")
if _log:
    _vd = os.environ["GUARD_VALIDATION_DIR"]
    sys.path.insert(0, _vd)
    try:
        import edmd_acc_guard as _G
    finally:
        sys.path.remove(_vd)
    _real = _G.guard
    def _logged_guard(path, _real=_real, _log=_log):
        with open(_log, "a") as f:
            f.write(os.path.abspath(os.fspath(path)) + "\n")
        return _real(path)
    _G.guard = _logged_guard
