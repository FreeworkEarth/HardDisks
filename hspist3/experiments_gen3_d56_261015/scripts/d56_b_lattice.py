#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56, task 2 b; notes 261012 sec. 4.7.133): THE LATTICE START OF THE 14 STAGE 2 CELLS
(eta 0.69, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80 at N 900 and 1600; the stage 1 geometry: H = 10 sqrt(N/100), L0 = N pi / (8 H eta) per
side, the divider at L0, thickness 0.05) with the orientation pinned, in cluster/gen3_partA_261014/pa_orient.py's way:
  1. Per cell three starts of the binary (the held path: --experiment=energy_transfer, --gen3-exact-box, --gen3-seeding=lattice, 60 held
     steps; pa_orient.cmd unchanged) plus --gen3-dump-initial: the default start (no pin), --gen3-lattice-orientation=0 and =turned.
     Read from the log: the [EDMD3-RUN] arrangement per compartment (rows' orientation, equal or alternate rows, n_s x n_r, the spacings
     a_s and a_r, the wall clearance g, the smallest gap s, the vacancies), the [EDMD3-M4C] line (the orientation found per compartment),
     "[EDMD3-M4] INFEASIBLE", the exit code and the [EDMD3-HEALTH] record. Read from the dumped initial state: the smallest disk-disk
     gap per compartment, the smallest gap to the outer walls and to the divider's faces, the disks per side, the divider's centre
     against L0.
  2. A replica of g3_lattice_fit (00ALLINONE.c, M4/M4c) in Python, checked against the binary's arrangement and spacings for every start,
     then used for what the binary does not print: the anisotropic strain of the start, e = h / (h_ideal) - 1 with h the spacing between
     neighbouring close-packed lines and h_ideal = (sqrt 3 / 2) times the spacing along them (0 for a perfect triangular lattice;
     negative: the lines are pressed together), and the strain-free compartment width W* of the same arrangement at the same height
     (the width at which e = 0; the binary's own spacing rules, the clearance g solved as there), against the actual width W = L0 - t/2.
Writes into experiments_gen3_d56_261015/b_lattice/ (a new folder): runs/<cell>_<start>/ (stdout.log, initial.txt, the run's small
csv files), b_starts.csv and REPORT.md. At most 6 processes at once.
usage (from hspist3/): python3 experiments_gen3_d56_261015/scripts/d56_b_lattice.py [--out <a new folder instead>]
"""
import concurrent.futures as cf, math, os, re, subprocess, sys
import numpy as np, pandas as pd
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d56_261015")
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_partA_261014")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import pa_orient as PO                                     # cmd(): the start-only command of decision 53 item 2, unchanged
import tests_20260913 as T
BIN = os.path.join(HS, "experiments_gen3_d50_261014", "bin_0544b7e8-dirty", "00ALLINONE")
OUT = os.path.join(E, "b_lattice")
PX, TH = 24.0, 0.05                                        # px per sigma; the divider's thickness [sigma]
ETAS, NS = (0.69, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80), (900, 1600)
STARTS = (("default", None), ("pin0", "0"), ("pinturned", "turned"))
BASE_B = 20261017                                          # the seeds: run_seed(BASE_B, cell, start, 0); no other list holds this base


def cells():
    out = []
    for N in NS:
        H = 10.0 * math.sqrt(N / 100.0)
        for eta in ETAS:
            out.append(dict(cell=f"S2_N{N:04d}_eta{eta:.6f}", N=N, eta=eta, H=H, L0=N * math.pi / (8.0 * H * eta)))
    return out


# ---------------------------------------------------------------------------------------------- the replica of g3_lattice_fit
def lat_eval(ns, nr, alt, R, Lr, Ls, g):
    us, ur = Ls - 2.0 * (R + g), Lr - 2.0 * (R + g)
    if not (us >= 0.0) or not (ur >= 0.0): return -math.inf, 0.0, 0.0, 0.0
    if alt and (ns < 2 or nr < 2): return -math.inf, 0.0, 0.0, 0.0
    a_s = us / (ns - 1) if ns > 1 else 0.0
    if alt: a_r = ur / (nr - 1)
    else: a_r = ur / ((nr - 1) + (0.5 if ns > 1 else 0.0)) if nr > 1 else (2.0 * ur if ns > 1 else 0.0)
    d = math.inf
    if nr > 1: d = a_r
    if ns > 1: d = min(d, math.sqrt(a_s * a_s + 0.25 * a_r * a_r))
    if ns > 2 and 2.0 * a_s < d: d = 2.0 * a_s
    gw = g
    if ns == 1: gw = min(g, 0.5 * Ls - R)
    if ns == 1 and nr == 1: gw = min(0.5 * Ls - R, 0.5 * Lr - R)
    return d - 2.0 * R, a_s, a_r, gw


def sites(ns, nr, alt): return ((ns + 1) // 2) * nr + (ns // 2) * (nr - 1) if alt else ns * nr


def cryst(o, ns, nr, a_s, a_r):
    along = 1 if ns < 2 else (0 if nr < 2 else int(a_r < 2.0 * a_s)); return 0 if along == int(o == 1) else 1


def solve_g(ns, nr, alt, R, Lr, Ls):
    sv = lat_eval(ns, nr, alt, R, Lr, Ls, 0.0)[0]
    if not (sv > 0.0): return None
    g = 0.0
    if math.isfinite(sv):
        for _ in range(200):
            gn = 0.5 * lat_eval(ns, nr, alt, R, Lr, Ls, g)[0]
            if not (gn > 0.0): g = 0.0; break
            if abs(gn - g) <= 1e-15 * (1.0 + gn): g = gn; break
            g = gn
    else: g = min(0.5 * Ls - R, 0.5 * Lr - R)
    return g


def lattice_fit(n, R, Wx, Hy, cryst_req=-1):
    best = None
    for o in (0, 1):
        Lr, Ls = (Hy, Wx) if o == 0 else (Wx, Hy)
        for alt in (0, 1):
            for ns in range(1, n + 1):
                nr = (n + ns - 1) // ns
                if alt:
                    if ns < 2: continue
                    nr = 2
                    while sites(ns, nr, 1) < n: nr += 1
                g = solve_g(ns, nr, alt, R, Lr, Ls)
                if g is None: continue
                sv, a_s, a_r, gw = lat_eval(ns, nr, alt, R, Lr, Ls, g)
                if not (sv > 0.0) or not (gw > 0.0): continue
                cr = cryst(o, ns, nr, a_s, a_r)
                if cryst_req >= 0 and cr != cryst_req: continue
                nv = sites(ns, nr, alt) - n
                if best is None or sv > best["s"] or (sv == best["s"] and nv < best["nvac"]):
                    best = dict(orient=o, alt=alt, ns=ns, nr=nr, nvac=nv, a_s=a_s, a_r=a_r, g=g, s=sv, cryst=cr)
    return best


def strain(b):
    """h / h_ideal - 1: h the spacing between neighbouring close-packed lines, h_ideal = (sqrt 3 / 2) x the spacing along a line"""
    along = 1 if b["ns"] < 2 else (0 if b["nr"] < 2 else int(b["a_r"] < 2.0 * b["a_s"]))
    return b["a_s"] / (math.sqrt(3) / 2 * b["a_r"]) - 1 if along else (b["a_r"] / 2) / (math.sqrt(3) / 2 * 2 * b["a_s"]) - 1


def at_width(b, R, W, Hy):
    o = b["orient"]; Lr, Ls = (Hy, W) if o == 0 else (W, Hy)
    g = solve_g(b["ns"], b["nr"], b["alt"], R, Lr, Ls)
    if g is None: return None
    sv, a_s, a_r, gw = lat_eval(b["ns"], b["nr"], b["alt"], R, Lr, Ls, g)
    return dict(b, a_s=a_s, a_r=a_r, g=g, s=sv) if sv > 0 and gw > 0 else None


def w_star(b, R, W, Hy):
    """the compartment width at which this arrangement (n_s, n_r, rows, orientation) is strain free, H fixed: a scan of 0.75 W .. 1.35 W
    for the sign change of e where the arrangement exists (a positive gap), then bisection"""
    def e(w):
        x = at_width(b, R, w, Hy); return strain(x) if x else math.nan
    ws = np.linspace(0.75 * W, 1.35 * W, 241); es = [e(w) for w in ws]
    for k in range(len(ws) - 1):
        if math.isfinite(es[k]) and math.isfinite(es[k + 1]) and es[k] * es[k + 1] <= 0:
            lo, hi, el = ws[k], ws[k + 1], es[k]
            for _ in range(80):
                mid = 0.5 * (lo + hi); em = e(mid)
                if (em > 0) == (el > 0): lo, el = mid, em
                else: hi = mid
            return 0.5 * (lo + hi)
    return math.nan


# ---------------------------------------------------------------------------------------------- the binary's starts
SIDE = re.compile(r"(left|right)=\[orient (\w+) rows (\w+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) g (\S+) s (\S+) n_vac (\d+) jitter (\S+) px\]")


def run_start(c, ci, si, start, pin):
    d = os.path.join(OUT, "runs", f"{c['cell']}_{start}"); os.makedirs(d)
    seed = T.run_seed(BASE_B, ci, si, 0)
    a = PO.cmd(BIN, d, c["N"], f"{c['H']:.17g}", f"{c['L0']:.17g}", f"{c['L0']:.17g}", seed, pin) + [f"--gen3-dump-initial={d}/initial.txt"]
    with open(os.path.join(d, "stdout.log"), "w") as o:
        rc = subprocess.run(a, stdout=o, stderr=subprocess.STDOUT, env=PO.ENV, cwd=d).returncode
    log = open(os.path.join(d, "stdout.log"), errors="replace").read()
    r = dict(cell=c["cell"], N=c["N"], eta=c["eta"], H=c["H"], L0=c["L0"], start=start, seed=seed, rc=rc,
             infeasible=int("[EDMD3-M4] INFEASIBLE" in log), clean=int(bool(re.search(r"^\[EDMD3-HEALTH\] .*: clean=1 ", log, re.M))))
    m4c = re.search(r"^\[EDMD3-M4C\] lattice orientation requested: .*?; left (.+?), right (.+)$", log, re.M)
    r["m4c_left"], r["m4c_right"] = (m4c.group(1), m4c.group(2)) if m4c else ("", "")
    hdr = re.search(r"^\[EDMD3-RUN\].*$", log, re.M)
    for m in (SIDE.finditer(hdr.group(0)) if hdr else []):
        k = m.group(1)
        r.update({f"{k}_orient": m.group(2), f"{k}_rows": m.group(3), f"{k}_ns": int(m.group(4)), f"{k}_nr": int(m.group(5)),
                  f"{k}_a_s": float(m.group(6)) / PX, f"{k}_a_r": float(m.group(7)) / PX, f"{k}_g": float(m.group(8)) / PX,
                  f"{k}_s": float(m.group(9)) / PX, f"{k}_nvac": int(m.group(10))})
        o = 0 if m.group(2) == "parallel" else 1
        r[f"{k}_cryst"] = "0" if cryst(o, int(m.group(4)), int(m.group(5)), float(m.group(6)), float(m.group(7))) == 0 else "turned"
    ini = os.path.join(d, "initial.txt")
    if rc == 0 and os.path.exists(ini): r.update(check_initial(ini, c))
    return r


def check_initial(p, c):
    L = open(p).read().split("\n"); assert L[0] == "gen3-initial v1", L[0]
    N, NL = int(L[1].split()[1]), int(L[1].split()[3]); W = float.fromhex(L[2].split()[1]); Hh = float.fromhex(L[3].split()[1])
    R = float.fromhex(L[4].split()[1]); cen, th = (float.fromhex(v) for v in L[5].split()[1:3])
    xy = np.array([[float.fromhex(v) for v in l.split()[:2]] for l in L[6:6 + N]])
    left = xy[:, 0] < cen; out = dict(n_left=int(left.sum()), n_right=int((~left).sum()), divider_minus_L0=(cen / PX) - c["L0"],
                                      box_w=W / PX, box_h=Hh / PX)
    gaps = []
    for side in (left, ~left):
        q = xy[side]; dd = np.hypot(q[:, None, 0] - q[None, :, 0], q[:, None, 1] - q[None, :, 1]); np.fill_diagonal(dd, np.inf)
        gaps.append((dd.min() - 2 * R) / PX)
    cross = np.hypot(xy[left][:, None, 0] - xy[~left][None, :, 0], xy[left][:, None, 1] - xy[~left][None, :, 1]).min() - 2 * R
    wall = min((xy[:, 0] - R).min(), (W - xy[:, 0] - R).min(), (xy[:, 1] - R).min(), (Hh - xy[:, 1] - R).min()) / PX
    div = min((cen - th / 2 - xy[left, 0] - R).min(), (xy[~left, 0] - R - (cen + th / 2)).min()) / PX
    out.update(min_gap_left=gaps[0], min_gap_right=gaps[1], min_gap_across=cross / PX, min_gap_walls=wall, min_gap_divider=div,
               no_overlap=int(min(gaps[0], gaps[1], cross / PX, wall, div) > 0))
    return out


def main():
    global OUT
    if "--out" in sys.argv: OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    os.makedirs(os.path.join(OUT, "runs")); C = cells()
    jobs = [(c, ci, si, st, pin) for ci, c in enumerate(C) for si, (st, pin) in enumerate(STARTS)]
    with cf.ThreadPoolExecutor(6) as ex: res = list(ex.map(lambda j: run_start(*j), jobs))
    rows, mism = [], []
    for r in res:
        c = [x for x in C if x["cell"] == r["cell"]][0]; R = 0.5 * PX
        W = c["L0"] * PX - 0.5 * TH * PX; Hp = c["H"] * PX; req = {"default": -1, "pin0": 0, "pinturned": 1}[r["start"]]
        b = lattice_fit(c["N"] // 2, R, W, Hp, req)              # the replica in px, as the binary computes (W = c - t/2, c = the box's centre)
        if b is None:
            r.update(replica="infeasible"); rows.append(r); continue
        same = ("left_ns" in r and r["left_ns"] == b["ns"] and r["left_nr"] == b["nr"] and r["left_rows"] == ("alternate" if b["alt"] else "equal")
                and r["left_orient"] == ("parallel" if b["orient"] == 0 else "perpendicular") and abs(r["left_a_s"] - b["a_s"] / PX) < 1e-6
                and abs(r["left_a_r"] - b["a_r"] / PX) < 1e-6 and abs(r["left_s"] - b["s"] / PX) < 1e-6)
        if "left_ns" in r and not same: mism.append(f"{r['cell']} {r['start']}: the binary {r['left_ns']}x{r['left_nr']} {r['left_rows']}, the replica {b['ns']}x{b['nr']} alt {b['alt']}")
        ws = w_star(b, R, W, Hp) / PX; W /= PX
        r.update(replica="same as the binary" if same else ("the binary did not start" if "left_ns" not in r else "DIFFERS"),
                 strain_pct=100 * strain(b), W=W, W_star=ws, W_over_Wstar_pct=100 * (W / ws - 1) if math.isfinite(ws) else math.nan,
                 aspect=W / c["H"], aspect_star=ws / c["H"] if math.isfinite(ws) else math.nan, lattice_cryst="0" if b["cryst"] == 0 else "turned",
                 eta_star=c["N"] * math.pi / 4 / ((2 * ws + TH) * c["H"]) if math.isfinite(ws) else math.nan)
        rows.append(r)
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "b_starts.csv"), index=False, float_format="%.10g")
    rep = ["# d56_b_lattice.py: the lattice start of the 14 Stage 2 cells (decision 56 task 2 b)", "",
           f"binary: {os.path.relpath(BIN, HS)}; per cell three starts (default, pin 0, pin turned), the start only (pa_orient.cmd), "
           "plus --gen3-dump-initial. Geometry: H = 10 sqrt(N/100), L0 = N pi / (8 H eta) per side, the divider at L0.", "",
           "## Per start: the checks", "",
           "| cell | start | exit, health | M4C (left / right) | left = right | smallest gaps in the dumped start: disks, walls, divider [sigma] | N_L / N_R | divider - L0 [sigma] |",
           "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r.get("rc") != 0 or "left_ns" not in r:
            rep.append(f"| {r['cell']} | {r['start']} | exit {r.get('rc')}{' INFEASIBLE' if r.get('infeasible') else ''} | {r.get('m4c_left', '')} / {r.get('m4c_right', '')} | - | - | - | - |")
            continue
        sym = all(r[f"left_{k}"] == r[f"right_{k}"] for k in ("orient", "rows", "ns", "nr", "nvac", "cryst", "a_s", "a_r"))
        rep.append(f"| {r['cell']} | {r['start']} | {r['rc']}, {'clean' if r['clean'] else 'NOT CLEAN'} | {r['m4c_left'] or '-'} / {r['m4c_right'] or '-'} | "
                   f"{'yes' if sym else 'NO'} | {min(r['min_gap_left'], r['min_gap_right']):.4f}, {r['min_gap_walls']:.4f}, {r['min_gap_divider']:.4f} | "
                   f"{r['n_left']} / {r['n_right']} | {r['divider_minus_L0']:+.1e} |")
    nf = sum(1 for r in rows if r.get("rc") == 0 and r.get("clean") and r.get("no_overlap"))
    rep += ["", f"{nf} of {len(rows)} starts exit 0 with a clean health record and no overlap; the replica of g3_lattice_fit gives the binary's "
            f"arrangement and spacings in {sum(1 for r in rows if r.get('replica') == 'same as the binary')} of {len(rows)}" + "".join(f"\n- {x}" for x in mism), "",
            "## Per cell: the arrangement, its strain and where it would be strain free", "",
            "| cell | aspect W/H | default picks | pin 0: arrangement, vac, gap s, strain e, eta* | pin turned: arrangement, vac, gap s, strain e, eta* |",
            "|---|---|---|---|---|"]
    for c in C:
        rr = {r["start"]: r for r in rows if r["cell"] == c["cell"]}
        def ent(r): return (f"{r['left_orient'][:4]} {r['left_rows'][:3]} {r['left_ns']}x{r['left_nr']}, {r['left_nvac']}, {r['left_s']:.3f}, "
                            f"{r['strain_pct']:+.1f} %, {r['eta_star']:.4f}") if "left_ns" in r else "-"
        rep.append(f"| {c['cell']} | {rr['default'].get('aspect', math.nan):.4f} | {rr['default'].get('left_cryst', '-')} | {ent(rr['pin0'])} | {ent(rr['pinturned'])} |")
    rep += ["", "Entries: rows' orientation (para/perp), rows (equ/alt), n_s x n_r; crystal 0: close-packed lines along x; turned: along y "
            "(psi6 orientation 30 deg). vac: the vacancies (sites - N/2). gap s: the lattice's smallest surface gap before the jitter [sigma]. "
            "e: the spacing between neighbouring close-packed lines over (sqrt 3 / 2) x the spacing along them, minus 1 (negative: the lines "
            "pressed together, the crystal compressed across them and stretched along them). eta*: the packing fraction at which the same "
            "arrangement in the same height H is strain free (the compartment width W* that makes e = 0, by the binary's own spacing rules; "
            "eta* = N pi / (4 (2 W* + t) H)). W = L0 - t/2, t = 0.05."]
    open(os.path.join(OUT, "REPORT.md"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep))


if __name__ == "__main__":
    main()
