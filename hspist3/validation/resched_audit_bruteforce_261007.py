#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10, plan-author decision part D1): independent high-precision recomputation of the
earliest contact for every event the schedule-equivalence audit (--resched-audit) reports. The engine prints each mismatch
as an [EDMD-AUDIT-EV] line with the particles' state at %.17g (exact doubles); here the earliest contact from that state is
solved with Python's decimal module at 60 significant digits:
  disk-disk  |r + v tau| = 2R, earliest root, approaching pairs only (r.v < 0); an overlapping approaching pair is overdue (tau = 0)
  wall       gap / closing speed (overdue if the gap is <= 0 while closing)
  divider    face gap / closing speed in the divider frame (constant divider velocity between events, as the engine assumes)
  piston     likewise
and compared with the legacy and the heap times. For a disk-disk pair the relative discriminant (b^2 - v^2 c)/b^2 is printed:
near 0 the pair grazes and the root is ill-conditioned (a rounding-level change of the state moves it far).
usage (from hspist3/): python3 validation/resched_audit_bruteforce_261007.py <log> [<log> ...]
"""
import re, sys
from decimal import Decimal as Dec, getcontext
getcontext().prec = 60
EV = re.compile(r"\[EDMD-AUDIT-EV\] (\w+) type=(\w+) (.*)")


def parse(line):
    m = EV.search(line)
    if not m: return None
    kind, typ, rest = m.groups(); d = dict(kind=kind, type=typ)
    for tok in rest.split():
        k, v = tok.split("=", 1)
        d[k] = [Dec(u) for u in v.split(",")] if "," in v else (Dec(v) if k not in ("a", "b", "d") else int(v))
    return d


def true_time(d):
    """(absolute earliest contact time or None, note) from the dumped state."""
    now, R = d["now"], d["R"]; x, y, vx, vy = d["A"]
    t = d["type"]
    if t == "AB":
        bx, by, bvx, bvy = d["B"]; rx, ry = bx - x, by - y; ux, uy = bvx - vx, bvy - vy
        b = rx * ux + ry * uy; c = rx * rx + ry * ry - (2 * R) ** 2; vv = ux * ux + uy * uy
        if b >= 0: return None, "separating (r.v >= 0): no contact"
        if c < 0: return now, "overlapping and approaching: overdue (tau = 0)"
        disc = b * b - vv * c
        if disc <= 0: return None, f"miss (discriminant {disc:.3e} <= 0)"
        tau = (-b - disc.sqrt()) / vv
        return now + tau, f"relative discriminant {disc / (b * b):.3e}"
    if t in ("WL", "WR", "WB", "WT"):
        W, H = d["boxW"], d["boxH"]
        gap, sp = {"WL": (x - R, -vx), "WR": (W - R - x, vx), "WB": (y - R, -vy), "WT": (H - R - y, vy)}[t]
        if sp <= 0: return None, "not closing"
        return (now, "overdue (gap <= 0)") if gap <= 0 else (now + gap / sp, f"gap {gap:.3e} px")
    if t in ("DL", "DR"):
        cx, u, th = d["div_x"], d["div_vx"], d["div_th"]
        if t == "DL":
            gap, rel = (cx - th / 2) - (x + R), vx - u
        else:
            gap, rel = (x - R) - (cx + th / 2), u - vx
        if rel <= 0: return None, "not closing in the divider frame"
        return (now + gap / rel, f"gap {gap:.3e} px, closing speed {rel:.3e}") if gap > 0 else (None, f"gap {gap:.3e} <= 0")
    if t in ("PL", "PR"):
        gap, rel = ((x - R) - d["pistonL_x"], -vx) if t == "PL" else (d["pistonR_x"] - (x + R), vx)
        if rel <= 0: return None, "not closing"
        return now + gap / rel, f"gap {gap:.3e} px"
    return None, "unknown type"


def report(lines, out=print, table=True):
    """Per reported event: the true earliest contact, both schedule times against it, and the error RELATIVE to the
    prediction horizon (t_true - now): rounding of a double prediction gives ~1e-16 to 1e-13 of the horizon, a wrong event
    gives order 1. Summary: the largest relative errors, and the number of events no true contact supports."""
    rows = [p for p in (parse(l) for l in lines) if p]
    if not rows:
        out("(no [EDMD-AUDIT-EV] lines: no mismatch was reported)"); return 0
    if table:
        out("| kind | type | a | b/d | now | horizon t_true - now | t_legacy - t_true | t_heap - t_true | heap error / horizon | brute force |")
        out("|---|---|---|---|---|---|---|---|---|---|")
    worst_h = worst_l = 0.0; unsupported = 0; hmax = 0.0
    for d in rows:
        tt, note = true_time(d)
        if tt is None:
            unsupported += 1
            if table: out(f"| {d['kind']} | {d['type']} | {d['a']} | {d.get('b', d.get('d', '-'))} | {float(d['now']):.6f} | - | "
                          f"{d['t_legacy']} | {d['t_heap']} | - | NO TRUE CONTACT: {note} |")
            continue
        hz = tt - d["now"]; hmax = max(hmax, float(hz))
        el = None if str(d["t_legacy"]).lower() == "nan" else abs(d["t_legacy"] - tt)
        eh = None if str(d["t_heap"]).lower() == "nan" else abs(d["t_heap"] - tt)
        rel = lambda e: float(e / hz) if (e is not None and hz > 0) else float("nan")
        if eh is not None and hz > 0: worst_h = max(worst_h, rel(eh))
        if el is not None and hz > 0: worst_l = max(worst_l, rel(el))
        if table:
            out(f"| {d['kind']} | {d['type']} | {d['a']} | {d.get('b', d.get('d', '-'))} | {float(d['now']):.6f} | {float(hz):.4e} | "
                f"{'absent' if el is None else format(float(d['t_legacy'] - tt), '+.3e')} | {'absent' if eh is None else format(float(d['t_heap'] - tt), '+.3e')} | "
                f"{rel(eh):.2e} | {note} |")
    out(f"\nbrute-force summary over {len(rows)} reported events: largest |t_heap - t_true| / horizon = {worst_h:.2e}; "
        f"largest |t_legacy - t_true| / horizon = {worst_l:.2e}; longest horizon {hmax:.3e} px-time; events with no true contact: {unsupported}")
    return len(rows)


if __name__ == "__main__":
    ls = []
    for p in sys.argv[1:]: ls += open(p, errors="ignore").read().split("\n")
    report(ls)
