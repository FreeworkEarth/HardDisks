#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 38 item 6; 261012 sec. 4.7.80): THE M5 SOLID-CELL JOB, the runner. The registered M5 rule
(sec. 4.7, gate plan item 5: "the double run's contact maxima and energy drift within the 80-bit run's by the floor of item 3, at one
fluid and one solid cell; decorrelation time reported") was read on K2-2's part (d), whose cells are fluid (N 100 at pi/8, N 400 at
eta 0.70). This job adds the solid cell of decision 38 item 6: N 400, eta 0.74, M 500, four seeds per build, the koa (double) and the
koa-ld (80-bit long double) builds of a643bc7 that ~/harddisks_gen3 holds since runsheet step 2 (not rebuilt).
The protocol is part (d)'s, unchanged: g3_rate_ld_koa.py's sos() (cluster/gen3_koa_261009, imported, not copied): the lattice start
(M4) and the exact box, 2000 steps held, 200 predicted periods recorded, HD_CONTACT_AUDIT=1, the koa and the koa-ld run of a seed with
the same command and seed; its record() reads every run's record. Seeds: run_seed(20261219, 2, 400, r), r = 0..3 -- part (d)'s base,
a stream it did not use (part (c) used stream 0, part (d) stream 1); checked distinct from those 20 seeds before anything runs.
The 8 runs at once (one core each). Per run it prints the row part (d) printed: clean, events, the event hash, the largest contact gap,
u_t, gap / u_t, argmax nu. The rule itself is read on the Mac after the fetch (m5_solid_verdict_261012.py).
usage (inside the job, from ~/harddisks_gen3_k2/hspist3): python3 cluster/gen3_m5_solid_261012/m5_solid_koa.py --bin B --bin-ld B_LD --out OUT
"""
import argparse, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import g3_rate_ld_koa as RL                                              # sos(), run(), record(), argmax_nu() as K2-2 used them
import tests_20260913 as T
N, ETA, M, SEEDS, STREAM, HOLD = 400, 0.74, 500, 4, 2, 2000


def seeds():
    old = {T.run_seed(RL.BASE, 0, k, 0) for k in range(12)} | {T.run_seed(RL.BASE, 1, n, r) for n in (100, 400) for r in range(4)}
    new = [T.run_seed(RL.BASE, STREAM, N, r) for r in range(SEEDS)]
    if len(set(new)) != SEEDS or set(new) & old: sys.exit("STOP: the solid cell's seeds are not distinct from K2-2's 20")
    return new


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--bin-ld", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    O = os.path.abspath(a.out); S = seeds()
    H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * ETA)
    c0 = RL.sos(a.bin, None, N, ETA, M, S[0], HOLD)
    stride = int(next(x for x in c0 if x.startswith("--speed-sound-log-stride=")).split("=")[1])
    nu_true = T.kr_cs(ETA) * T.k_root(M / N) / (2 * math.pi * T.l_eff(L0))
    print(f"== M5, the solid cell: N {N}, eta {ETA}, M {M} (alpha {M / N:g}), H {H:g}, L0 {L0:.17g}; {SEEDS} seeds per build "
          f"(run_seed({RL.BASE}, {STREAM}, {N}, r)): {', '.join(str(s) for s in S)}; held {HOLD} steps, 200 predicted periods")
    print(f"   trace stride {stride} steps (sos(): KR nu with tests_20260913.N_SIDE = {T.N_SIDE}); about "
          f"{1 / (nu_true * T.DT_SIGMA * stride):.1f} samples per KR period at alpha {M / N:g} (information)\n", flush=True)
    jobs, procs = [], []
    for r, s in enumerate(S):
        for tag, b in (("koa", a.bin), ("koa-ld", a.bin_ld)):
            d = os.path.join(O, "ld", f"N{N}_M{M}_eta{ETA:.2f}", f"seed{r}_{tag}")
            c = [x if not x.startswith("--speed-sound-run-dir=") else f"--speed-sound-run-dir={d}" for x in RL.sos(b, None, N, ETA, M, s, HOLD)]
            jobs.append((r, tag, d)); procs.append(RL.run(c, d))
    t0 = time.time()
    for p in procs: p.wait()
    print(f"the 8 runs ended after {time.time() - t0:.0f} s; exit codes {', '.join(str(p.returncode) for p in procs)}\n")
    print("| N | M | eta | seed | build | clean | events | event hash | max contact gap [px] | u_t | gap / u_t [px/unit] | argmax nu |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
    bad = 0
    for r, tag, d in jobs:
        q = RL.record(d)
        if not q: print(f"| {N} | {M} | {ETA} | {r} | {tag} | no run record | | | | | | |"); bad += 1; continue
        bad += q["clean"] != "1"
        print(f"| {N} | {M} | {ETA} | {r} | {tag} | {q['clean']} | {q['events']} | {q['hash']} | {q['gap']:.3g} | {q['u_t']:.3g} | "
              f"{q['gap'] / q['u_t']:.3g} | {RL.argmax_nu(d):.6f} |")
    if bad: print(f"\nSTOP: {bad} of the 8 runs have no clean run record"); sys.exit(1)


if __name__ == "__main__":
    main()
