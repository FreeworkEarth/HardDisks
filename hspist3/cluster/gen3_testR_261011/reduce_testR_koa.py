#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.57): THE KOA REDUCTION OF ONE TEST R CELL FOLDER m_<M>, after
reduce_B.py (merge_koa.py runs both): the per-trajectory tables the Mac needs, so the fetch carries no traces.
  nud_runs.csv     per trajectory: run, seed, nu_arg (red_nu.csv's registered argmax), nu_d, fit_ok, nu_ref, dt, n. nu_d is REGISTERED
                   (decision 28): fit_acf exactly as testG_koa_posthoc_261011.py at 60909289 computes it -- the imported
                   validation/resched_testT_followup_261007.fit_acf on each trajectory's ACF (acf_runs.npz, reduce_B.py's), the ACFs
                   truncated to the shortest of the folder, nu_d = omega / 2 pi, a fit that raises gives NaN (fit_ok 0) -- with the start
                   value of decision 28: the mass's pooled argmax mean of its OWN engine, i.e. the mean of red_nu.csv's nu over this folder
                   (the post hoc script pooled both engines; the registration names the engine's own).
  runs_health.csv  per ##RUN section of run.log: run, seed, node, seconds, build line count, gen3 run-record count and clean flag,
                   engine line, gen2's divider-policy line, health-pattern lines (testG_worker_koa.sh's HEALTH; for gen3 without its own
                   [EDMD3-HEALTH] record), the contact-audit line and its largest |contact distance| [px], the end-of-record energies
                   (E_div, E_tot of HD_KE_TRACE=1's "[EDMD-ENERGY] 2 end of record").
  failed_runs.csv  per .failed_run<r>_* folder (a run the worker refused; after the merge renamed per chunk): its name, run, the health-
                   pattern lines and the last line of its stdout.log.
Prints one line: the counts. Exit 0, or 2 with a STOP line (an input missing or inconsistent).
usage: python3 cluster/gen3_testR_261011/reduce_testR_koa.py <.../m_<M>>
"""
import glob, math, os, re, sys, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import resched_testT_followup_261007 as FU                      # fit_acf: the sec. 4.4.13 item 4a model (as the post hoc script)
warnings.filterwarnings("ignore")
HEALTH = re.compile(r"EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue")      # testG_worker_koa.sh's pattern
CONTACT = re.compile(r"^\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), outer walls (\S+), "
                     r"divider (\S+), pistons (\S+)$", re.M)
ENERGY = re.compile(r"^\[EDMD-ENERGY\] 2 end of record\s+t=\S+ E_gas=\S+ E_div=(\S+) E_tot=(\S+)", re.M)


def nud_table(d):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
    if red["dt"].nunique() != 1: sys.exit(f"STOP: more than one dt in {d}/red_nu.csv")
    z = np.load(os.path.join(d, "acf_runs.npz"))
    if sorted(z.files) != sorted(f"run{r}" for r in red["run"]): sys.exit(f"STOP: acf_runs.npz and red_nu.csv list other runs in {d}")
    nl = min(len(z[f"run{r}"]) for r in red["run"])
    dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    nd, ok = [], []
    for r in red["run"]:
        a = np.asarray(z[f"run{r}"][:nl], dtype=float)
        try: nd.append(FU.fit_acf(a, dt, nu_ref)[4] / (2 * math.pi)); ok.append(1)
        except Exception: nd.append(math.nan); ok.append(0)
    out = pd.DataFrame(dict(run=red["run"], seed=red["seed"], nu_arg=red["nu"], nu_d=nd, fit_ok=ok, nu_ref=nu_ref, dt=dt, n=red["n"]))
    out.to_csv(os.path.join(d, "nud_runs.csv"), index=False, float_format="%.17g")
    return len(out), int(len(out) - sum(ok))


def health_table(d):
    log = open(os.path.join(d, "run.log"), errors="replace").read()
    rows = []
    for sec in re.split(r"^##RUN ", log, flags=re.M)[1:]:
        h = re.match(r"\S+ \S+ run (\d+) seed (\d+) node (\S+) \((\d+) s\)", sec)
        if not h: continue
        recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", sec, re.M)
        other = sum(1 for l in sec.splitlines() if HEALTH.search(l) and not l.startswith("[EDMD3-HEALTH]"))
        c = CONTACT.search(sec); e = ENERGY.search(sec)
        rows.append(dict(run=int(h.group(1)), seed=int(h.group(2)), node=h.group(3), seconds=int(h.group(4)),
                         built_lines=len(re.findall(r"^\[EDMD3\] built #", sec, re.M)), gen3_records=len(recs),
                         gen3_clean=int(recs == ["1"]), engine_gen3_line=int(" engine=gen3 " in sec),
                         gen2_policy_line=int("[EDMD-RESCHED] divider events: minimal (default)" in sec),
                         edmd3_lines=len(re.findall(r"^\[EDMD3", sec, re.M)), health_lines=other,
                         contact_line=int(bool(c)), contact_max_px=(max(float(v) for v in c.groups()[1:]) if c else math.nan),
                         E_div=(float(e.group(1)) if e else math.nan), E_tot=(float(e.group(2)) if e else math.nan)))
    pd.DataFrame(rows).to_csv(os.path.join(d, "runs_health.csv"), index=False, float_format="%.17g")
    return len(rows)


def failed_list(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, ".failed_run*"))):
        m = re.match(r"^\.failed_run(\d+)_", os.path.basename(f)); so = os.path.join(f, "stdout.log")
        txt = open(so, errors="replace").read() if os.path.exists(so) else ""
        rows.append(dict(folder=os.path.basename(f), run=int(m.group(1)) if m else -1,
                         health_lines=sum(1 for l in txt.splitlines() if HEALTH.search(l)),
                         last_line=(txt.strip().splitlines() or ["(no stdout.log)"])[-1][:200]))
    pd.DataFrame(rows, columns=["folder", "run", "health_lines", "last_line"]).to_csv(os.path.join(d, "failed_runs.csv"), index=False)
    return len(rows)


def main():
    if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]): sys.exit("usage: reduce_testR_koa.py <.../m_<M>>")
    d = sys.argv[1]
    for f in ("red_nu.csv", "acf_runs.npz", "run.log"):
        if not os.path.exists(os.path.join(d, f)): print(f"STOP: {d} has no {f}"); sys.exit(2)
    n, nf = nud_table(d); nh = health_table(d); nfl = failed_list(d)
    print(f"{n} trajectories: nud_runs.csv ({nf} fit failures), runs_health.csv ({nh} ##RUN sections), failed_runs.csv ({nfl} failed runs)")


if __name__ == "__main__":
    main()
