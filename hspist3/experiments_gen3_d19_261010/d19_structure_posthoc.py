#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.45; plan-author decision 19, item 1) -- POST HOC, existing snapshots and virial blocks only, no
runs; information. WHAT IS THE HIGHER-Z GROUP? For N = 1600 at eta 0.740 and 0.780 and N = 400 at eta 0.740, every compartment that
has a snapshot (seeds 0 and 1 of P6's turned and 0-deg starts, P4's lattice start A and compressed start B, P5's compressed starts):
the snapshots of the last quarter of the hold (1.5e4 <= t <= 2e4 sigma-time: t = 16000, 18000, 20000; for P4a, which held 1e5, the
same window), each read from its G3SNAP1 file (positions and the divider centre in sigma).
THE RULES (per compartment c: the disks with x < divider centre are c = 0):
  neighbours   the Delaunay triangulation of the compartment's disks, edges longer than 1.4 sigma removed (psi6's cutoff)
  orientation  psi6(j) = the mean of exp(6 i theta_jk) over j's neighbours; arg(sum over the compartment) / 6 in degrees, (-30, 30]
  lines        the close-packed axis is y for a turned compartment (|orientation| > 15 deg) and x otherwise; the lines are counted from
               the disks' coordinate across that axis (x for lines along y): a histogram in bins of 0.025 sigma, smoothed by a Gaussian
               of 0.1 sigma, its local maxima above a quarter of the largest; in four equal bands along the line axis (a band-to-band
               change of the count is an extra half-line, a dislocation) and over the whole length; each disk belongs to the nearest
               whole-length maximum, which gives the disks per line (min, max, histogram)
               [a first version counted chains of bonds within 15 deg of an axis; at eta 0.74 thermal motion broke the chains into
               fragments (25-32 "lines" for 22), so it was replaced before any number from it was used]
  coordination the 5- and 7-coordinated disks among the INTERIOR disks (centre >= 1.2 sigma from every wall and from the divider face)
  vacancies / interstitials  per line its capacity = the sites that fit along the line axis at the median in-line spacing a:
               floor((extent - 2R) / a) + 1 for a line whose first disk lies within a/4 of the wall, floor((extent - 2R - a/2) / a) + 1
               for a shifted line; vacancies = the sum of (capacity - disks) where positive, interstitials = the sum of (disks -
               capacity) where positive
  per quantity the median over the three snapshots (the counts of a snapshot fluctuate thermally)
  Z            the dt-weighted Z_c of the virial blocks 300-399 (Q4)
  Z group      per cell, among the compartments at |orientation| > 15 deg (turned), split at the largest gap of their Z (as
               p5_compartments_posthoc.py); compartments at 0 deg are their own group, "0 deg"
Last, per cell: for every structural number, whether its ranges in the lower and the higher group overlap; the numbers whose ranges
are disjoint separate the groups.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d19_261010/d19_structure_posthoc.py
  python3 experiments_gen3_d19_261010/d19_structure_posthoc.py --selftest
"""
import argparse, cmath, glob, math, os, struct, sys
import numpy as np
from scipy.spatial import Delaunay

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
P4D = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
P5D = os.path.join(HS, "experiments_gen3_p5_261010", "data_06e7248")
P6D = os.path.join(HS, "experiments_gen3_p6_261010", "data_f72d764")
CUT, ANG, INNER, TIMES = 1.4, 15.0, 1.2, (16000.0, 18000.0, 20000.0)
CELLS = ("N1600_eta0.740", "N1600_eta0.780", "N400_eta0.740")


def starts(cell):
    out = [("P6 turned lattice (clean)", os.path.join(P6D, cell, "turned"))]
    if os.path.isdir(os.path.join(P6D, cell, "zero")): out.append(("P6 0-deg lattice (addition)", os.path.join(P6D, cell, "zero")))
    out.append(("P4 A (present M4 family)", os.path.join(P4D, cell, "A")))
    out.append(("P4 B compressed u 0.01", os.path.join(P4D, cell, "B")))
    for u in ("0.001", "0.0001"):
        d = os.path.join(P5D, f"{cell}_u{u}", "B")
        if os.path.isdir(d): out.append((f"P5 B compressed u {u}", d))
    return out


def read_snaps(p):
    b = open(p, "rb").read()
    assert b[:8] == b"G3SNAP1\n", "magic"
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, Tsn = struct.unpack_from("<4d", b, 16); o = 48; out = []
    while o < len(b):
        t, = struct.unpack_from("<d", b, o); ph, k = struct.unpack_from("<ii", b, o + 8); o += 16
        dx = struct.unpack_from(f"<{k}d", b, o); o += 8 * k
        xy = np.array(struct.unpack_from(f"<{2 * N}d", b, o)).reshape(N, 2); o += 16 * N
        out.append((t, ph, dx, xy))
    return (N, nd, W, H, R), out


def zq4(d):
    import csv
    f = glob.glob(os.path.join(d, "virial_blocks_*_run0.csv"))
    if len(f) != 1: return None
    bl = [x for x in csv.DictReader(open(f[0])) if x["phase"] == "hold" and 300 <= int(x["block"]) < 400]
    if len(bl) != 100: return None
    W = sum(float(x["dt_sigma"]) for x in bl)
    return [sum(float(x[f"Z_{c}"]) * float(x["dt_sigma"]) for x in bl) / W for c in (0, 1)]


def edges(P):
    tri = Delaunay(P); E = set()
    for s in tri.simplices:
        for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
            a, b = (a, b) if a < b else (b, a)
            if np.hypot(*(P[a] - P[b])) <= CUT: E.add((a, b))
    return E


def components(n, E):
    par = list(range(n))
    def f(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for a, b in E: par[f(a)] = f(b)
    out = {}
    for i in range(n): out.setdefault(f(i), []).append(i)
    return list(out.values())


def peaks(v, lo, hi):
    """the local maxima of the smoothed histogram of v over [lo, hi]: bins 0.025 sigma, Gaussian 0.1 sigma, above a quarter of the largest"""
    nb = max(8, int(round((hi - lo) / 0.025))); h, e = np.histogram(v, bins=nb, range=(lo, hi))
    k = np.arange(-16, 17); g = np.exp(-0.5 * (k * 0.025 / 0.1) ** 2); hs = np.convolve(h, g / g.sum(), mode="same")
    c = 0.5 * (e[1:] + e[:-1]); top = hs.max() if len(hs) else 0.0
    return [c[i] for i in range(1, len(hs) - 1) if hs[i] > hs[i - 1] and hs[i] >= hs[i + 1] and hs[i] > 0.25 * top]


def analyse(P, x0, x1, H, R):
    """one compartment: P (n, 2) positions in sigma, the compartment spans x0..x1 (the walls/divider faces), height H"""
    n = len(P); E = edges(P); sig = 2.0 * R
    nb = [[] for _ in range(n)]
    for a, b in E: nb[a].append(b); nb[b].append(a)
    s = 0j
    for j in range(n):
        if nb[j]: s += sum(cmath.exp(6j * math.atan2(*(P[k] - P[j])[::-1])) for k in nb[j]) / len(nb[j])
    orient = math.degrees(cmath.phase(s)) / 6.0 if abs(s) > 0 else float("nan")
    cp = "y" if abs(orient) > 15 else "x"
    across, along = (P[:, 0], P[:, 1]) if cp == "y" else (P[:, 1], P[:, 0])
    alo, ahi = (x0, x1) if cp == "y" else (0.0, H)                 # the extent across the lines
    llo, lhi = (0.0, H) if cp == "y" else (x0, x1)                 # the extent along the lines
    bands = []
    for q in range(4):
        b0, b1 = llo + q * (lhi - llo) / 4, llo + (q + 1) * (lhi - llo) / 4
        m = (along >= b0) & (along < b1); bands.append(len(peaks(across[m], alo, ahi)))
    pk = np.array(peaks(across, alo, ahi)); L = len(pk)
    per = []; vac = inter = 0
    if L:
        lab = np.argmin(np.abs(across[:, None] - pk[None, :]), axis=1)
        sp = []
        for i in range(L):
            yy = np.sort(along[lab == i]); per.append(len(yy))
            if len(yy) > 1: sp.extend(np.diff(yy))
        a = float(np.median(sp)) if sp else float("nan")
        for i in range(L):
            yy = np.sort(along[lab == i])
            if not len(yy) or not math.isfinite(a): continue
            first = yy[0] - (llo + R)
            cap = (math.floor((lhi - llo - 2 * R) / a + 1e-9) + 1) if first < a / 4 else (math.floor((lhi - llo - 2 * R - a / 2) / a + 1e-9) + 1)
            d = cap - len(yy); vac += max(0, d); inter += max(0, -d)
    inner = [j for j in range(n) if P[j][0] - x0 >= INNER * sig and x1 - P[j][0] >= INNER * sig and P[j][1] >= INNER * sig and H - P[j][1] >= INNER * sig]
    c5 = sum(1 for j in inner if len(nb[j]) == 5); c7 = sum(1 for j in inner if len(nb[j]) == 7)
    return dict(n=n, orient=orient, cp=cp, L=L, bands=bands, lmin=min(per) if per else 0, lmax=max(per) if per else 0, hist=per, c5=c5, c7=c7,
                inner=len(inner), vac=vac, inter=inter, lx=0, ly=0)


def compartments(snap, head):
    (N, nd, W, H, R) = head; t, ph, dx, xy = snap
    c = dx[0]; th = 0.05                                        # the divider's thickness [sigma] (--wall-thickness=0.05)
    L = xy[xy[:, 0] < c]; Rr = xy[xy[:, 0] >= c]
    return [(L, 0.0, c - th / 2), (Rr, c + th / 2, W)], H, R


def run_one(d):
    f = glob.glob(os.path.join(d, "snapshots_*_run0.bin"))
    if len(f) != 1: return None
    head, snaps = read_snaps(f[0])
    sel = [s for s in snaps if s[1] == 0 and any(abs(s[0] - t) < 1.0 for t in TIMES)]
    if not sel: return None
    per = [[], []]
    for s in sel:
        comps, H, R = compartments(s, head)
        for k, (P, x0, x1) in enumerate(comps): per[k].append(analyse(P, x0, x1, H, R))
    out = []
    for k in (0, 1):
        rs = per[k]; med = lambda key: float(np.median([r[key] for r in rs]))
        med_min = lambda rr: float(np.median([min(r["bands"]) for r in rr])); med_max = lambda rr: float(np.median([max(r["bands"]) for r in rr]))
        last = rs[-1]
        z = complex(0)
        for r in rs: z += cmath.exp(6j * math.radians(r["orient"]))
        out.append(dict(times=[s[0] for s in sel], orient=math.degrees(cmath.phase(z)) / 6.0, cp=last["cp"], L=med("L"), bands=last["bands"],
                        bmin=med_min(rs), bmax=med_max(rs),
                        lmin=med("lmin"), lmax=med("lmax"), hist=last["hist"], c5=med("c5"), c7=med("c7"), inner=med("inner"),
                        vac=med("vac"), inter=med("inter"), n=last["n"]))
    return out


def hist_txt(h):
    from collections import Counter
    c = Counter(h); return " ".join(f"{k}x{c[k]}" for k in sorted(c))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: return selftest()
    print("# Decision 19, item 1 (POST HOC): what is the higher-Z group? -- printed by experiments_gen3_d19_261010/d19_structure_posthoc.py; "
          "existing snapshots and virial blocks only, no runs; information\n")
    print(f"Rules: neighbours = Delaunay edges <= {CUT} sigma (compartment's own disks); lines = chains of neighbour bonds within {ANG:.0f} deg of a "
          f"box axis, >= 3 disks; coordination over interior disks (>= {INNER} sigma from walls and the divider face); vacancies / interstitials "
          "against ceil(L/2) n + floor(L/2) (n - 1) sites (L lines along the close-packed axis, n the longest line); medians over the snapshots "
          "at t = 16000, 18000, 20000; Z = Q4 of the virial blocks\n")
    allrows = {}
    for cell in CELLS:
        rows = []
        for lab, sd in starts(cell):
            for r in (0, 1):
                d = os.path.join(sd, f"seed{r}"); res = run_one(d); z = zq4(d)
                if not res or not z: continue
                for k in (0, 1): rows.append(dict(start=lab, seed=r, c=k, z=z[k], **res[k]))
        turned = [x for x in rows if abs(x["orient"]) > 15]
        zs = sorted(x["z"] for x in turned)
        g, i = max((zs[j + 1] - zs[j], j) for j in range(len(zs) - 1)); cut = 0.5 * (zs[i] + zs[i + 1])
        for x in rows: x["group"] = "0 deg" if abs(x["orient"]) <= 15 else ("lower" if x["z"] < cut else "higher")
        allrows[cell] = (rows, cut, g)
        print(f"## {cell}\n\nZ groups among the turned compartments split at their largest gap {g:.4f} (cut {cut:.4f})\n")
        print("| start | seed | comp. | Z (Q4) | group | orientation [deg] | close-packed axis | lines in the 4 bands (last snapshot) | lines, whole length | disks per line "
              "min / max | histogram (last snapshot) | 5-coord. | 7-coord. | interior disks | vacancies | interstitials |\n" + "|---" * 16 + "|")
        for x in rows:
            print(f"| {x['start']} | {x['seed']} | {x['c']} | {x['z']:.4f} | {x['group']} | {x['orient']:+.1f} | {x['cp']} | {' '.join(str(v) for v in x['bands'])} | "
                  f"{x['L']:.0f} | {x['lmin']:.0f} / {x['lmax']:.0f} | {hist_txt(x['hist'])} | {x['c5']:.0f} | {x['c7']:.0f} | {x['inner']:.0f} | "
                  f"{x['vac']:.0f} | {x['inter']:.0f} |")
        print("\nper structural number, the ranges in the lower and the higher group (turned compartments):\n")
        print("| number | lower group: min - max | higher group: min - max | separates (disjoint ranges) |\n|---|---|---|---|")
        lo = [x for x in rows if x["group"] == "lower"]; hi = [x for x in rows if x["group"] == "higher"]
        for key, name in (("L", "lines along the close-packed axis (whole length)"), ("bmin", "fewest lines in a band"), ("bmax", "most lines in a band"), ("lmin", "shortest line"), ("lmax", "longest line"), ("c5", "5-coordinated"),
                          ("c7", "7-coordinated"), ("vac", "vacancies"), ("inter", "interstitials")):
            if not lo or not hi: continue
            a0, a1 = min(x[key] for x in lo), max(x[key] for x in lo); b0, b1 = min(x[key] for x in hi), max(x[key] for x in hi)
            sep = a1 < b0 or b1 < a0
            print(f"| {name} | {a0:.0f} - {a1:.0f} | {b0:.0f} - {b1:.0f} | {'YES' if sep else 'no'} |")
        print()
    print("(POST HOC, information only; no verdict.)")


def selftest():
    """a perfect turned lattice (lines along y) and the same lattice with one line removed: lines, coordination, orientation"""
    ok = True
    d = 1.0; a_r, a_s = d, d * math.sqrt(3) / 2
    for L, label in ((10, "perfect 10 lines"), (9, "9 lines")):
        P = []
        for k in range(L):
            for j in range(20 - (k % 2)):
                P.append((0.5 + k * a_s, 0.5 + j * a_r + (0.5 * a_r if k % 2 else 0.0)))
        P = np.array(P); W = 0.5 + (L - 1) * a_s + 0.5; H = 0.5 + 19 * a_r + 0.5
        r = analyse(P, 0.0, W, H, 0.5)
        e = r["cp"] == "y" and r["L"] == L and r["lmax"] == 20 and r["lmin"] == 19 and abs(abs(r["orient"]) - 30.0) < 0.5 and r["c5"] == 0 and r["c7"] == 0 \
            and r["vac"] == 0 and r["inter"] == 0
        ok &= e
        print(f"selftest {label}: cp {r['cp']}, lines {r['L']}, per line {r['lmin']}-{r['lmax']}, orientation {r['orient']:+.2f}, 5/7 {r['c5']}/{r['c7']}, "
              f"vac {r['vac']}, inter {r['inter']}: {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
