#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.4 amendment b, sec. 4.7.6; milestone M2): the generation-3 contact tolerances,
derived from the time resolution, printed with their scale beside the contact errors measured in M1 and M2.

THE READING UNDER TEST (plan author): the contact error of an executed event is set by the time resolution -- ulp of the
origin-relative time, at most ulp(2^13) = 2^-39 internal units below the next origin shift -- times a speed, not by
ulp(d^2). If so, the measured max |gap| divided by the time quantum of the run is a SPEED of the order of the contact speeds,
the same for short runs (small times, small quantum) and long runs (times up to 2^13).

DERIVATION (first order, unit roundoff u_t/2 of a time below 2^14, u_t = ulp(2^13)):
  a contact predicted at t_p for t_c is stored rounded:          gap error <= v_n u_t / 2          (v_n: closing speed)
  at execution each disk is evaluated at t_c from its stamp:      <= (|v_i| + |v_j|) u_t / 2        (t_c - tau rounded)
  at prediction the partner was evaluated the same way:           <= |v_j| u_t / 2
  products and sums in cell-local coordinates (<= 32 px):         a few ulp(32 px) = 7e-15 px each (negligible)
  -> one executed contact:  |gap| <= (v_n + |v_i| + 2 |v_j|) u_t / 2 <= 1.5 v_ref u_t
  a pair re-predicted inside the overlap its own contact left (a third disk changes one of them at nearly the same time),
  with both evaluated again:  |gap| <= 1.5 v_ref u_t + v_ref u_t = 2.5 v_ref u_t.
  v_ref bounds every relative speed: with E the mechanical energy (kinetic + spring) plus the work done on the system,
  two unit-mass disks have |v_i| + |v_j| <= 2 sqrt(E); a disk and a body of mass M >= 1: <= sqrt(2 E (1 + 1/M)); a body of
  mass 0 adds its prescribed speed. So v_ref = sqrt(2 E (1 + 1/m_min)) + max |v_driven|, m_min = min(1, finite masses).
  c = |r|^2 - d^2 = 2 d gap + gap^2  ->  c_tol = K * 2 d * v_ref * u_t,   K = 4 > 2.5.
  Faces (outer walls, divider faces, pistons) compare ABSOLUTE positions (up to the box width) too:
  tol_face = K * v_ref * u_t + 8 ulp(box width).
usage (from hspist3/):  python3 validation/gen3_tolerances_261009.py
"""
import math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
M1 = os.path.join(HS, "experiments_gen3_m1_261008", "m1_audit_output.txt")
M2 = os.path.join(HS, "experiments_gen3_m2_261009", "m2_audit_output.txt")
M2Q = os.path.join(HS, "experiments_gen3_m2_261009", "m2_audit_quick_output.txt")
R, D = 12.0, 24.0                 # px
SHIFT = 8192.0                    # EDMD3_ORIGIN_SHIFT
K = 4.0                           # EDMD3_TOL_K
U_T = math.ulp(SHIFT)             # 2^-39
C_TOL_M1 = 64.0 * D * D * 2.220446049250313e-16
TOL_PAIR = max(1e-7, 1e-6 * D)    # experiment_validation.c


def v_ref(E, m_min=1.0, v_drv=0.0):
    return math.sqrt(2.0 * E * (1.0 + 1.0 / m_min)) + v_drv


def quantum(T_sig, t0=0.0):
    """the ulp of the largest origin-relative event time of a run of T sigma-time (24 units each) that starts at t0: an origin
    shift at 2^13 keeps every time below 2^13 + one event, so a run that reaches 2^13 has times just above it (ulp(2^13))"""
    tmax = t0 + 24.0 * T_sig
    return math.ulp(SHIFT) if tmax >= SHIFT else math.ulp(tmax)


def parse_m1(path):
    rows = []
    if not os.path.exists(path):
        return rows
    cell = None
    for line in open(path):
        m = re.match(r"### Cell (\w+): N = (\d+), .*T = (\d+) sigma-time, cell width (\d+) px", line)
        if m:
            cell = {"run": "M1", "cell": m.group(1) + ("" if m.group(4) == "32" else f" ({m.group(4)} px cells)"), "N": int(m.group(2)), "T": float(m.group(3))}
            continue
        m = re.match(r"contact audit: gen3 (\d+) events, max \|gap\| pairs (\S+) px, walls (\S+) px; gen2 (\d+) events, max \|gap\| pairs (\S+) px, walls (\S+) px", line)
        if m and cell:
            cell.update(pair=float(m.group(2)), face=float(m.group(3)), gen2_pair=float(m.group(5)), gen2_face=float(m.group(6)),
                        E=float(cell["N"]), c_min=None)
            rows.append(cell); cell = None
    return rows


def parse_m2(path, label):
    rows = []
    if not os.path.exists(path):
        return rows
    cell = None
    for line in open(path):
        m = re.match(r"### Cell (\w+): ", line)
        if m:
            cell = {"run": label, "cell": m.group(1)}
            continue
        if cell is None:
            continue
        m = re.match(r"N = (\d+), box .*T = (\d+) sigma-time", line)
        if m:
            cell.update(N=int(m.group(1)), T=float(m.group(2)))
            m0 = re.search(r"gen3 loaded at t = (\S+) units", line)
            cell["t0"] = float(m0.group(1)) if m0 else 0.0
        m = re.match(r"\| gen3 \(A\) \| (\d+) \| (\S+) \| (\S+) \| (\S+) \| (\S+) \|", line)
        if m:
            cell.update(pair=float(m.group(2)), face=max(float(m.group(3)), float(m.group(4)), float(m.group(5))))
        m = re.match(r"\| gen2 \| (\d+) \| (\S+) \| (\S+) \| (\S+) \| (\S+) \|", line)
        if m:
            cell.update(gen2_pair=float(m.group(2)), gen2_face=max(float(m.group(3)), float(m.group(4)), float(m.group(5))))
        m = re.search(r"contact_c_min=(\S+)", line)
        if m:
            cell["c_min"] = float(m.group(1))
        m = re.search(r"E_bound = (\S+) kT, m_min = (\S+), v_ref = (\S+) px/unit, K = \S+ -> c_tol = (\S+) px\^2 .*tol_face = (\S+) px", line)
        if m:
            cell.update(E=float(m.group(1)), v_ref=float(m.group(3)), c_tol=float(m.group(4)), tol_face=float(m.group(5)))
            rows.append(cell); cell = None
    return rows


def main():
    print("# Generation-3 contact tolerances from the time resolution (261012 sec. 4.7.4 amendment b, sec. 4.7.6), "
          "printed by validation/gen3_tolerances_261009.py\n")
    print("## 1. The time resolution\n\n| origin-relative time [units] | ulp [units] |\n|---|---|")
    for t in (480.0, 960.0, 2400.0, 4096.0, math.nextafter(SHIFT, 0.0), SHIFT, 16383.0):
        print(f"| {t:.17g} | {math.ulp(t):.4g} |")
    print(f"\nu_t = ulp(EDMD3_ORIGIN_SHIFT) = ulp(2^13) = 2^-39 = {U_T:.6g} units: every executed time is below 2^13 + one event, "
          f"so below 2^14 (a time predicted more than 2^13 ahead keeps the ulp it was stored with: 3.6e-12 up to 2^15).")
    print(f"M1's threshold: c_tol = 64 ulp(d^2) = 64 x {math.ulp(D * D):.3g} = {C_TOL_M1:.4g} px^2, i.e. a gap of "
          f"{C_TOL_M1 / (2 * D):.3g} px; the rounding of c itself (ulp(d^2) = {math.ulp(D * D):.3g} px^2) is that scale, the time "
          f"rounding (v u_t ~ {U_T:.2g} px at v = 1 px/unit) is not.\n")

    print("## 2. The measured contact errors against the time quantum of each run\n")
    print("u_run = ulp of the largest event time the run reaches (the origin shift caps it at ulp(2^13)). If the reading holds, "
          "max |gap| / u_run is a speed [px/unit] of the order of the contact speeds (thermal speed sqrt(2) px/unit at kT = 1), "
          "whatever u_run is.\n")
    rows = parse_m1(M1) + parse_m2(M2Q, "M2 quick") + parse_m2(M2, "M2")
    print("| run | cell | N | T [sigma-time] (start) | u_run [units] | gen3 max pair gap [px] | / u_run [px/unit] | gen3 max face gap [px] | / u_run | "
          "gen2 max pair gap [px] | / u_run | 2 d gap [px^2] | / M1 c_tol | / M2 c_tol | most negative c of a contact_now [px^2] | its |c| / M1 c_tol (> 1: M1 would count an overlap repair) |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        if "pair" not in r:
            continue
        q = quantum(r["T"], r.get("t0", 0.0)); c = 2 * D * r["pair"]
        ct = r.get("c_tol") or K * 2 * D * v_ref(r.get("E", r["N"])) * U_T
        cm = r.get("c_min")
        print(f"| {r['run']} | {r['cell']} | {r['N']} | {r['T']:.0f}{' (t0 = ' + format(r.get('t0', 0.0), 'g') + ' units)' if r.get('t0') else ''} | {q:.3g} | {r['pair']:.3g} | {r['pair'] / q:.3g} | {r['face']:.3g} | "
              f"{r['face'] / q:.3g} | {r.get('gen2_pair', float('nan')):.3g} | {r.get('gen2_pair', float('nan')) / q:.3g} | {c:.3g} | "
              f"{c / C_TOL_M1:.3g} | {c / ct:.3g} | {'-' if cm is None else f'{cm:.3g}'} | {'-' if cm is None else f'{abs(cm) / C_TOL_M1:.3g}'} |")
    print("\n(gen2 is not cell-local: its positions are absolute (ulp(box) ~ 1e-13 px) and its time is not origin-relative, so its "
          "errors scale with the absolute time and box; printed for comparison only.)\n")

    print("## 3. The derived tolerances at kT = 1 (E = N kT in two dimensions, no driven body, m_min = 1)\n")
    print(f"| N | E [kT] | v_ref = 2 sqrt(E) [px/unit] | c_tol = {K:g} x 2 d v_ref u_t [px^2] | as a gap c_tol / 2d [px] | "
          f"tol_face at box 481.2 px [px] | at 3721.2 px [px] | validator: c at tol_pair / c_tol |\n|---|---|---|---|---|---|---|---|")
    for N in (100, 400, 900, 1600):
        E = float(N); v = v_ref(E); ct = K * 2 * D * v * U_T
        tf = lambda W: K * v * U_T + 8 * math.ulp(W)
        print(f"| {N} | {E:.0f} | {v:.4g} | {ct:.4g} | {ct / (2 * D):.4g} | {tf(481.2):.4g} | {tf(3721.2):.4g} | {2 * D * TOL_PAIR / ct:.3g} |")
    print(f"\nThe validator's pair tolerance {TOL_PAIR:.3g} px corresponds to c = {2 * D * TOL_PAIR:.3g} px^2: the derived c_tol lies "
          f"between the rounding scale (section 2, at most ~1e-10 px^2) and a missed collision (>= the validator scale) by more than "
          f"an order of magnitude on each side.")


if __name__ == "__main__":
    main()
