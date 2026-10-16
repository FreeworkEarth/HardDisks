#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE LATTICE TOOLS OF THE STAGE 2 PILOT P2 (no simulation;
imported by gen_p2_koa_261016.py, p2_reduce.py and the Mac analysis).
  1. THE REPLICA of g3_lattice_fit (00ALLINONE.c, M4 and M4c) in the binary's own units (px, 24 per sigma; the disk radius 12 px; a
     compartment of width W = c - t/2 where c = the box centre L0, t = 0.05 sigma = 1.2 px), as in experiments_gen3_d56_261015/scripts/
     d56_b_lattice.py, which checked it against the binary's arrangement and spacings in 42 of 42 starts (Report D56-1 (b)):
       lattice_fit(n, R, W, H, cryst_req)  the binary's choice (the largest smallest gap; the fewest vacancies on a tie), cryst_req -1
                                           (default), 0 or 1 (turned);
       strain(arr)                         e = h / ((sqrt 3 / 2) d) - 1: d the spacing along the close-packed lines, h the spacing between
                                           them (0 for a perfect triangular lattice; negative: the lines pressed together);
       at_width(arr, R, W, H), w_star(arr, R, W, H)   the same arrangement at another width; the width at which e = 0 (H fixed).
  2. READERS: the initial state of --gen3-dump-initial ("gen3-initial v1": box coordinates in px, C hex floats) and the position
     snapshots of --gen3-snapshots ("G3SNAP1": sigma, box coordinates); both returned in sigma.
  3. THE STRAIN OF A CONFIGURATION, per compartment (config_strain): the disks whose centres lie at least `margin` (2 sigma) from every
     wall and from the divider's faces (the bulk); per disk the local psi6 over its neighbours within rc = 1.35 a on its own side
     (a = sqrt(2 / (sqrt 3 rho)), rho the compartment's number density, the divider's position as recorded); the compartment's
     orientation theta0 = arg(the mean local psi6 of its ordered bulk disks, |psi6| >= 0.7) / 6; its line direction theta_L = the bond
     direction theta0 + k 60 deg closest to the x or the y axis (0 deg: lines along x; 90 deg: along y, the turned crystal);
     the in-grain disks: ordered bulk disks whose own orientation lies within 7.5 deg of theta0 (mod 60);
     their bonds within 15 deg of a lattice direction: d = the mean length of the bonds along theta_L, l = of the four others;
     e = sqrt((4 l^2 / d^2 - 1) / 3) - 1, which equals strain() for the binary's arrangement (the along bonds are the in-row spacing a_r
     when a_r < 2 a_s, else 2 a_s; the others sqrt(a_s^2 + a_r^2 / 4)). Thermal motion enters d and l alike.
usage: imported; python3 p2_lattice.py --selftest checks the estimator on the binary's lattices (perfect and jittered).
"""
import gzip, math, struct, sys
import numpy as np

PX = 24.0                    # px per sigma (PIXELS_PER_SIGMA)
TH = 0.05                    # the divider's thickness [sigma]


# ---------------------------------------------------------------------------------------------- 1. the replica of g3_lattice_fit
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
    along = 1 if b["ns"] < 2 else (0 if b["nr"] < 2 else int(b["a_r"] < 2.0 * b["a_s"]))
    return b["a_s"] / (math.sqrt(3) / 2 * b["a_r"]) - 1 if along else (b["a_r"] / 2) / (math.sqrt(3) / 2 * 2 * b["a_s"]) - 1


def at_width(b, R, W, Hy):
    o = b["orient"]; Lr, Ls = (Hy, W) if o == 0 else (W, Hy)
    g = solve_g(b["ns"], b["nr"], b["alt"], R, Lr, Ls)
    if g is None: return None
    sv, a_s, a_r, gw = lat_eval(b["ns"], b["nr"], b["alt"], R, Lr, Ls, g)
    return dict(b, a_s=a_s, a_r=a_r, g=g, s=sv) if sv > 0 and gw > 0 else None


def w_star(b, R, W, Hy):
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


def fmt_arr(b):
    return f"{'para' if b['orient'] == 0 else 'perp'} {'alt' if b['alt'] else 'equ'} {b['ns']}x{b['nr']} {'0' if b['cryst'] == 0 else 'turned'}"


# ---------------------------------------------------------------------------------------------- 2. readers
def opener(p, mode="rt"): return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode.replace("t", "") if "b" in mode else mode)


def read_dump(p):
    """--gen3-dump-initial: 'gen3-initial v1', N, N_L, boxW, boxH, radius, divider centre and thickness (px, C hex floats), then
    x y vx vy per disk; returned in sigma"""
    with opener(p) as fh: L = fh.read().split("\n")
    if L[0] != "gen3-initial v1": raise ValueError(f"{p}: not a gen3-initial v1 file")
    N, NL = int(L[1].split()[1]), int(L[1].split()[3])
    W, H, R = (float.fromhex(L[k].split()[1]) / PX for k in (2, 3, 4))
    cen, th = (float.fromhex(v) / PX for v in L[5].split()[1:3])
    xy = np.array([[float.fromhex(v) for v in l.split()[:2]] for l in L[6:6 + N]]) / PX
    return dict(N=N, NL=NL, W=W, H=H, R=R, cen=cen, th=th, xy=xy)


def read_snapshots(p):
    """--gen3-snapshots: 'G3SNAP1', N, dividers, box width and height, radius, interval (sigma); per snapshot t, phase (0 hold, 1 record),
    the divider centres, N x (x, y)"""
    with opener(p, "rb") as fh: b = fh.read()
    if b[:8] != b"G3SNAP1\n": raise ValueError(f"{p}: not a G3SNAP1 file")
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, every = struct.unpack_from("<dddd", b, 16); off = 48; frames = []
    while off < len(b):
        t, ph, ndv = struct.unpack_from("<dii", b, off); off += 16
        cen = struct.unpack_from(f"<{ndv}d", b, off); off += 8 * ndv
        xy = np.frombuffer(b, dtype="<f8", count=2 * N, offset=off).reshape(N, 2).copy(); off += 16 * N
        frames.append(dict(t=t, phase=ph, cen=cen[0] if ndv else math.nan, xy=xy))
    return dict(N=N, W=W, H=H, R=R, every=every, frames=frames)


# ---------------------------------------------------------------------------------------------- 3. the strain of a configuration
def _wrap60(a): return (np.asarray(a) + 30.0) % 60.0 - 30.0


def config_strain(xy, cen, W, H, th=TH, margin=2.0, rc_fac=1.35, psi_min=0.7, grain_tol=7.5, ang_tol=15.0):
    """per compartment (0 left, 1 right) a dict: e, theta0 [deg, mod 60 in [0, 60)], line [deg: 0 or 90 for lines along x or y], d, l,
    n (disks), n_bulk, n_used, psi6_local_bulk (mean |psi6| of the bulk), ordered_bulk (fraction >= psi_min), n_along, n_across"""
    out = []
    for k, side in enumerate((xy[:, 0] < cen, xy[:, 0] > cen)):
        q = xy[side]; n = len(q)
        x0, x1 = (0.0, cen - th / 2) if k == 0 else (cen + th / 2, W)
        rho = n / ((x1 - x0) * H); a = math.sqrt(2.0 / (math.sqrt(3.0) * rho)); rc = rc_fac * a
        d = q[None, :, :] - q[:, None, :]; r = np.hypot(d[..., 0], d[..., 1]); nb = (r < rc) & (r > 0)
        ang = np.arctan2(d[..., 1], d[..., 0])
        z = np.where(nb, np.exp(6j * ang), 0).sum(axis=1); cnt = nb.sum(axis=1); psi = np.where(cnt > 0, z / np.maximum(cnt, 1), 0)
        bulk = (q[:, 0] - x0 >= margin) & (x1 - q[:, 0] >= margin) & (q[:, 1] >= margin) & (H - q[:, 1] >= margin)
        ordered = bulk & (np.abs(psi) >= psi_min)
        res = dict(n=n, n_bulk=int(bulk.sum()), psi6_local_bulk=float(np.abs(psi[bulk]).mean()) if bulk.any() else math.nan,
                   ordered_bulk=float(ordered.sum() / max(1, bulk.sum())))
        if ordered.sum() < 10:
            res.update(e=math.nan, theta0=math.nan, line=math.nan, d=math.nan, l=math.nan, n_used=int(ordered.sum()), n_along=0, n_across=0)
            out.append(res); continue
        th0 = math.degrees(np.angle(psi[ordered].mean())) / 6.0
        own = np.degrees(np.angle(psi)) / 6.0
        grain = ordered & (np.abs(_wrap60(own - th0)) <= grain_tol)
        dirs = [(th0 + 60.0 * m) % 180.0 for m in range(3)]
        dist_axis = [min(abs(t - 0.0), abs(t - 90.0), abs(t - 180.0)) for t in dirs]
        tl = dirs[int(np.argmin(dist_axis))]
        bang = np.degrees(ang) % 180.0
        rel = (bang - tl + 90.0) % 180.0 - 90.0                      # bond angle relative to the line direction, in [-90, 90)
        use = grain[:, None] & nb
        along = use & (np.abs(rel) <= ang_tol)
        across = use & (np.abs(np.abs(rel) - 60.0) <= ang_tol)
        dd, ll = float(r[along].mean()) if along.any() else math.nan, float(r[across].mean()) if across.any() else math.nan
        e = math.sqrt(max(0.0, (4.0 * ll * ll / (dd * dd) - 1.0) / 3.0)) - 1.0 if math.isfinite(dd) and math.isfinite(ll) else math.nan
        res.update(e=e, theta0=float(0.0 if th0 % 60.0 > 60.0 - 1e-6 else th0 % 60.0), line=float(0.0 if min(abs(tl), abs(tl - 180.0)) < 45.0 else 90.0), d=dd, l=ll,
                   n_used=int(grain.sum()), n_along=int(along.sum()), n_across=int(across.sum()))
        out.append(res)
    return out


def selftest():
    """the estimator on the binary's own lattices: perfect and with the default jitter (f = 0.25), against strain()"""
    rng = np.random.default_rng(1); R = 0.5 * PX; worst = 0.0; lines = []
    for N, H, eta, req in ((900, 30.0, 0.76, 1), (900, 30.0, 0.72, 0), (1600, 40.0, 0.80, 1), (900, 30.0, 0.78, 0)):
        L0 = N * math.pi / (8.0 * H * eta); W = L0 * PX - 0.5 * TH * PX; Hp = H * PX
        b = lattice_fit(N // 2, R, W, Hp, req)
        for f in (0.0, 0.25):
            pts = []
            for side in (0, 1):
                x0 = 0.0 if side == 0 else L0 * PX + 0.5 * TH * PX
                Lr, Ls = (Hp, W) if b["orient"] == 0 else (W, Hp)
                jit = f * min(b["s"], 2.0 * b["g"])
                for kk in range(b["ns"]):
                    for j in range(b["nr"] - 1 if (b["alt"] and kk % 2) else b["nr"]):
                        S = R + b["g"] + kk * b["a_s"]; T = R + b["g"] + j * b["a_r"] + (0.5 * b["a_r"] if (b["ns"] > 1 and kk % 2) else 0.0)
                        rr, ph = jit * math.sqrt(rng.random()), 2 * math.pi * rng.random()
                        x, y = (S, T) if b["orient"] == 0 else (T, S)
                        pts.append((x0 + x + rr * math.cos(ph), y + rr * math.sin(ph)))
            xy = np.array(pts) / PX
            res = config_strain(xy, L0, 2 * L0, H)
            e_rep = strain(b)
            for k in (0, 1):
                worst = max(worst, abs(res[k]["e"] - e_rep)) if f == 0.0 else worst
            lines.append(f"N {N} eta {eta} pin {req} jitter {f}: {fmt_arr(b)}, e replica {100 * e_rep:+.3f} %, estimator "
                         f"{100 * res[0]['e']:+.3f} / {100 * res[1]['e']:+.3f} %, theta0 {res[0]['theta0']:.2f}, line {res[0]['line']:.0f}, "
                         f"used {res[0]['n_used']} of {res[0]['n_bulk']} bulk")
    return worst, lines


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        w, ls = selftest(); print("\n".join(ls)); print(f"largest |e_estimator - e_replica| on the perfect lattices (vacancies included): {w:.2e}")
