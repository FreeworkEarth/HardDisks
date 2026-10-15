#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 50 item 3 d; 261012 sec. 4.7.113; registration sec. 4.8.4): THE PART A REDUCER, option (ii):
one held trajectory's run folder -> one row.

It reads the record lines of the run's log:
  [EDMD3-HELD-GEOMETRY]   the geometry (gate C2's quantities: the box, the divider's centre and thickness, the compartments, L_eff, eta);
  [EDMD3-RUN]             the lattice of each compartment (the start), T_eq;
  [EDMD3-HELD] t0, t1     the divider's per-face impulse counters at T_eq and at the last step of the hold (exact values in %a);
  [EDMD3-HELD] window     F_L, F_R over (t0, t1];
  [EDMD3-HEALTH]          clean, the event hash, t_end.
F_L = (J_L(t1) - J_L(t0)) / T_w and F_R = -(J_R(t1) - J_R(t0)) / T_w, T_w = (t1 - t0) / 24, are recomputed from the exact values and must
equal the window line's. The temperatures come from the released tail's trace as reduce_AF.py reads them (the KE_gas_left and
KE_gas_right means over the counts of SegCounts). With an event log (the check subset), it also forms reduce_AF.py's F_L, F_R, n_L,
n_R, u_wall_max and W_div over [t0, t1) and their differences from the readout.

FAILS CLOSED (exit 1, nothing written) when:
  - a record line is missing or doubled, or the run is not clean;
  - the readout's steps are not T_eq and T_eq + the window;
  - the counts are not N_s each, or the geometry is off by more than 1e-9 sigma;
  - T_w differs from the window steps x dt by more than 1e-9 (relative);
  - the held divider moved or worked (v_div, dx_div or dW_div not exactly 0);
  - with a log: the counts differ, the forces differ by more than 1e-9 (relative; the log prints 12 digits), or u_wall or W_div is not 0.
usage: pa_reduce.py <run folder> <out csv> <rel> <j> <s> <seed> <N_s> <T_eq steps> <window steps>
"""
import os, re, sys
import numpy as np, pandas as pd

PPS = 24.0                                   # internal units per sigma (and per sigma-time)
DT = float(np.float32(0.4)) / PPS            # the driver's step, --fixed-dt=0.4 as a float, in sigma-time
COLS = ["rel", "cell", "j", "s", "seed", "N", "N_s", "L0", "x_div_requested", "box", "H", "x_div", "t_div", "R", "w_L", "w_R",
        "L_eff_L", "L_eff_R", "eta", "dW_box", "dx_requested", "lattice_left", "lattice_right", "T_eq_header", "t0", "t1", "T_w",
        "F_L", "F_R", "n_L", "n_R", "v_div_t0", "v_div_t1", "dx_div", "dW_div", "T_L", "T_R", "clean", "hash", "t_end",
        "F_L_log", "F_R_log", "n_L_log", "n_R_log", "u_wall_max_log", "W_div_log", "dF_L_rel_log", "dF_R_rel_log"]


def die(msg):
    sys.exit(f"ABORT pa_reduce {sys.argv[1]}: {msg}")


def one(lines, prefix):
    L = [l for l in lines if l.startswith(prefix)]
    if len(L) != 1:
        die(f"{len(L)} lines starting '{prefix}' (exactly one expected)")
    return L[0]


def kv(text):
    return dict(re.findall(r"(\w+)=(\S+)", text))


def held(lines, which, step):
    l = one(lines, f"[EDMD3-HELD] {which} ")
    dec, ex = l.split(" exact ", 1)
    a, x = kv(dec), kv(ex)
    if int(a["step"]) != step:
        die(f"{which} read at step {a['step']}, expected {step}")
    return dict(step=int(a["step"]), t=float.fromhex(x["t"]), J_L=float.fromhex(x["J_L"]), J_R=float.fromhex(x["J_R"]),
                x_div=float.fromhex(x["x_div"]), n_L=int(a["n_L"]), n_R=int(a["n_R"]), v_div=float(a["v_div"]), W_div=float(a["W_div"]))


def main():
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    d, out, rel, j, s, seed, NS, teq, win = sys.argv[1:10]
    j, s, NS, teq, win = int(j), int(s), int(NS), int(teq), int(win)
    lines = open(os.path.join(d, "stdout.log"), errors="replace").read().splitlines()

    geo = one(lines, "[EDMD3-HELD-GEOMETRY] ")
    m = re.search(r"box=(\S+) x (\S+) sigma \(2 L0 (\S+), \|dW\| (\S+)\) x_div=(\S+) sigma \(requested (\S+), \|dx\| (\S+)\)", geo)
    if not m:
        die("the geometry line does not parse")
    box, H, twoL0, dWb, xd, xreq, dxr = (float(v) for v in m.groups())
    g = kv(geo[geo.index(") t=") + 2:])
    N, NL, NR = int(g["N"]), int(g["N_L"]), int(g["N_R"])
    if N != 2 * NS or NL != NS or NR != NS:
        die(f"N {N}, N_L {NL}, N_R {NR}; expected {2 * NS} = {NS} + {NS}")
    if not (dWb <= 1e-9 and dxr <= 1e-9):
        die(f"the geometry is off: |box - 2 L0| {dWb}, |x_div - requested| {dxr} sigma")

    hdr = one(lines, "[EDMD3-RUN] held ")
    lat = re.search(r" (left=\[.*?\]) (right=\[.*?\])", hdr)
    teq_h = re.search(r" T_eq=(\S+) sigma-time", hdr)
    L0 = float(re.search(r"held L0=(\S+) ", hdr).group(1))

    a, b = held(lines, "t0", teq), held(lines, "t1", teq + win)
    w = kv(one(lines, "[EDMD3-HELD] window "))
    Tw = (b["t"] - a["t"]) / PPS
    FL, FR = (b["J_L"] - a["J_L"]) / Tw, -(b["J_R"] - a["J_R"]) / Tw
    for k, v in (("T_w", Tw), ("F_L", FL), ("F_R", FR)):
        if abs(float(w[k]) - v) > 1e-15 * abs(v):
            die(f"the window line's {k} {w[k]} is not the recomputed {v!r}")
    if int(w["n_L"]) != b["n_L"] - a["n_L"] or int(w["n_R"]) != b["n_R"] - a["n_R"]:
        die("the window line's counts are not t1 - t0")
    if abs(Tw / (win * DT) - 1.0) > 1e-9:
        die(f"T_w {Tw!r} is not the window of {win} steps ({win * DT!r} sigma-time)")
    if not (a["v_div"] == 0.0 and b["v_div"] == 0.0 and float(w["dx_div"]) == 0.0 and float(w["dW_div"]) == 0.0 and b["x_div"] == a["x_div"]):
        die(f"the held divider moved or worked: v_div {a['v_div']}, {b['v_div']}, dx_div {w['dx_div']}, dW_div {w['dW_div']}")

    hl = one(lines, "[EDMD3-HEALTH] ")
    h = kv(hl.split(": ", 1)[1])
    if h.get("clean") != "1":
        die(f"the run record says clean={h.get('clean')}")

    tr = pd.read_csv(os.path.join(d, "tr.csv"), usecols=["Time", "KE_gas_left", "KE_gas_right", "SegCounts"], low_memory=False)
    if len(tr) == 0:
        die("the tail's trace has no rows")
    cnt = tr["SegCounts"].astype(str).str.split(";").iloc[0]
    nl, nr = int(cnt[0]), int(cnt[-1])
    TL, TR = tr["KE_gas_left"].mean() / nl, tr["KE_gas_right"].mean() / nr

    row = dict(rel=rel, cell=rel.rstrip("/").split("/")[-2], j=j, s=s, seed=int(seed), N=N, N_s=NS, L0=L0, x_div_requested=xreq, box=box,
               H=H, x_div=xd, t_div=float(g["t"]), R=float(g["R"]), w_L=float(g["w_L"]), w_R=float(g["w_R"]), L_eff_L=float(g["L_eff_L"]),
               L_eff_R=float(g["L_eff_R"]), eta=float(g["eta"]), dW_box=dWb, dx_requested=dxr,
               lattice_left=lat.group(1) if lat else "", lattice_right=lat.group(2) if lat else "",
               T_eq_header=float(teq_h.group(1)) if teq_h else np.nan, t0=a["t"] / PPS, t1=b["t"] / PPS, T_w=Tw, F_L=FL, F_R=FR,
               n_L=b["n_L"] - a["n_L"], n_R=b["n_R"] - a["n_R"], v_div_t0=a["v_div"], v_div_t1=b["v_div"], dx_div=float(w["dx_div"]),
               dW_div=float(w["dW_div"]), T_L=TL, T_R=TR, clean=1, hash=h.get("hash", ""), t_end=float(h.get("t_end", "nan")))
    ev = os.path.join(d, "ev.csv")
    if os.path.exists(ev):
        e = pd.read_csv(ev, usecols=["t_sigma", "kind", "u_wall", "dE", "dp"])
        dd = e[(e["kind"] == "D0") & (e["t_sigma"] >= a["t"] / PPS) & (e["t_sigma"] < b["t"] / PPS)]
        dp = dd["dp"].to_numpy()
        if (dp == 0).any():
            die("divider events with dp == 0 in the log (face undefined)")
        fl, fr = -dp[dp < 0].sum() / Tw, dp[dp > 0].sum() / Tw
        row.update(F_L_log=fl, F_R_log=fr, n_L_log=int((dp < 0).sum()), n_R_log=int((dp > 0).sum()),
                   u_wall_max_log=float(np.abs(dd["u_wall"]).max()) if len(dd) else 0.0, W_div_log=float(dd["dE"].sum()),
                   dF_L_rel_log=FL / fl - 1.0, dF_R_rel_log=FR / fr - 1.0)
        if row["n_L_log"] != row["n_L"] or row["n_R_log"] != row["n_R"]:
            die(f"the log's counts {row['n_L_log']}, {row['n_R_log']} differ from the readout's {row['n_L']}, {row['n_R']}")
        if abs(row["dF_L_rel_log"]) > 1e-9 or abs(row["dF_R_rel_log"]) > 1e-9:
            die(f"the log's forces differ from the readout's by {row['dF_L_rel_log']:.3g}, {row['dF_R_rel_log']:.3g}")
        if row["u_wall_max_log"] != 0.0 or row["W_div_log"] != 0.0:
            die(f"the log has u_wall_max {row['u_wall_max_log']}, W_div {row['W_div_log']}")
    pd.DataFrame([row], columns=COLS).to_csv(out, index=False, float_format="%.17g")


if __name__ == "__main__":
    main()
