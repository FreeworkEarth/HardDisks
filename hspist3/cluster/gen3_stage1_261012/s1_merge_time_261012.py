#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 45 item 5 e, addendum 1 item 2): THE MERGE JOB'S TIME AND WHERE IT WENT, on the Mac, from the
fetched files' modification times (nothing on KOA is read; nothing is written into the fetched folder). fetch_s1_koa.sh's rsync -a keeps
them, files and folders alike; KOA's Lustre keeps whole seconds. merge_s1_koa.py and s1_merge.sbatch write, in this order:
  scan     the job's checks and the scan (reads only); its end: <root>/merge_plan_<UTC time>.tsv;
  moves    the plan carried out; its end: merged/.merge_complete (its first line is that time);
  then per cell folder, in sorted order:
  reduce_B    reduce_B.py: red_nu.csv and acf_runs.npz in each m_<M>; its end: the latest of them;
  reduce_s1   reduce_s1_koa.py per m_<M>: nud_runs.csv, runs_health.csv, failed_runs.csv, acf_mean_fit.csv, s1_runs.csv, snap_runs.csv;
              its end: the latest of them in the cell;
  then per m_<M>:
  sha         the SHA-256 of each trace, appended to .sha256_uncompressed; its end: that file's time;
  gzip        gzip -9 of each trace; its end: the m_<M> folder's own time (its last trace replaced by the .gz);
  end         s1_merge.sbatch copies the job's log into job_logs/ as its last command: the log's time is the job's end.
The job's start is in no fetched file: --start 'YYYY-mm-dd HH:MM:SS' (UTC, from sacct) gives it; without it the scan is "unknown". The
order of the cells and folders is checked against the merge log's own reduce_B / reduce_s1 / compressed lines, and every phase must end
after the one before it (else the line is marked).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_merge_time_261012.py <group> [--root R] [--start 'YYYY-mm-dd HH:MM:SS']
"""
import argparse, datetime as DT, glob, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
S1_OUT = ("nud_runs.csv", "runs_health.csv", "failed_runs.csv", "acf_mean_fit.csv", "s1_runs.csv", "snap_runs.csv")
UTC = DT.timezone.utc
ZONES = {"UTC": 0, "GMT": 0, "HST": -10}                        # the zone .merge_complete names (merge_s1_koa.py's %Z); the plan's name uses the same clock


def mt(p): return int(os.stat(p).st_mtime)                     # whole seconds (Lustre)
def ts(t): return DT.datetime.fromtimestamp(t, UTC).strftime("%Y-%m-%d %H:%M:%S")
def hms(s): return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("group"); ap.add_argument("--root", default=os.path.join(HS, "experiments_gen3_stage1_koa_261012"))
    ap.add_argument("--start"); a = ap.parse_args()
    G = os.path.join(a.root, a.group); M = os.path.join(G, "merged"); base = os.path.join(M, "stage1", "gen3")
    if not os.path.isdir(base): sys.exit(f"STOP: no fetched folder {base}")
    print(f"# INTERIM -- group {a.group}: the merge job's time and where it went (s1_merge_time_261012.py, from the fetched files' times; UTC)\n")
    plans = sorted(glob.glob(os.path.join(G, "merge_plan_*.tsv"))); logs = sorted(glob.glob(os.path.join(G, "job_logs", f"gen3-s1-merge-{a.group}_*.out")))
    if len(plans) != 1 or len(logs) != 1: sys.exit(f"STOP: expected one plan and one merge log, found {len(plans)} and {len(logs)}")
    first = open(os.path.join(M, ".merge_complete")).readline().strip(); t_mc = mt(os.path.join(M, ".merge_complete"))
    zone = first[20:].strip(); tz = DT.timezone(DT.timedelta(hours=ZONES.get(zone, 0)))      # KOA writes UTC; the Mac rehearsal HST
    t_plan = mt(plans[0]); stamp = DT.datetime.strptime(re.search(r"merge_plan_(\d{8}_\d{6})\.tsv$", plans[0]).group(1), "%Y%m%d_%H%M%S").replace(tzinfo=tz)
    rec = DT.datetime.strptime(first[:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=tz)
    t_end = mt(logs[0]); job = re.search(r"_(\d+)\.out$", logs[0]).group(1)
    print(f"- the files keep the merge's times: the plan's name says {stamp.astimezone(UTC):%H:%M:%S} UTC, its file {ts(t_plan)[11:]}; .merge_complete "
          f"says '{first}', its file {ts(t_mc)[11:]} UTC (file minus record: {t_plan - int(stamp.timestamp()):+d} s and {t_mc - int(rec.timestamp()):+d} s)"
          + ("" if zone in ZONES else f" -- the zone '{zone}' is not known here; read as UTC"))
    # the merge log's order of the per-cell steps
    order = []
    for l in open(logs[0], errors="replace"):
        m = re.match(r"^reduce_B stage1/gen3/(\S+): exit (\d+)", l) or re.match(r"^reduce_s1 stage1/gen3/(\S+): exit (\d+)", l) \
            or re.match(r"^  stage1/gen3/(\S+): (\d+) traces compressed$", l)
        if m: order.append((l.split()[0] if not l.startswith("  ") else "gzip", m.group(1), int(m.group(2))))
    rows, prev, bad, seq = [], t_mc, [], []
    cells = [c for c in sorted(glob.glob(os.path.join(base, "*"))) if glob.glob(os.path.join(c, "m_*"))]
    for c in cells:
        cid = os.path.basename(c); ms = sorted(glob.glob(os.path.join(c, "m_*")))
        seq += [("reduce_B", cid)] + [("reduce_s1", f"{cid}/{os.path.basename(m)}") for m in ms] + [("gzip", f"{cid}/{os.path.basename(m)}") for m in ms]
        rb = max(mt(os.path.join(m, f)) for m in ms for f in ("red_nu.csv", "acf_runs.npz") if os.path.exists(os.path.join(m, f)))
        rs = max(mt(os.path.join(m, f)) for m in ms for f in S1_OUT if os.path.exists(os.path.join(m, f)))
        sha = gz = 0; p = rs; n_tr = 0
        for m in ms:
            s_end = mt(os.path.join(m, ".sha256_uncompressed")); g_end = max(mt(m), s_end)
            if not (p <= s_end <= g_end): bad.append(f"{cid}/{os.path.basename(m)}: sha or gzip out of order")
            sha += s_end - p; gz += g_end - s_end; p = g_end
        if not (prev <= rb <= rs): bad.append(f"{cid}: reduce_B or reduce_s1 out of order")
        rows.append(dict(cid=cid, n=len(ms), rB=rb - prev, rs=rs - rb, sha=sha, gz=gz, tot=p - prev)); prev = p
    tail = t_end - prev
    got = [(k, w) for k, w, _ in order]
    print(f"- the merge log ({os.path.relpath(logs[0], G)}, job {job}): {len(order)} step lines, "
          + ("in the order assumed here (cells sorted; reduce_B, then reduce_s1 per m_<M>, then sha and gzip per m_<M>)" if got == seq
             else f"NOT in the order assumed here ({len(seq)} expected) -- the split below is not valid")
          + f"; exit codes other than 0 in reduce_B / reduce_s1: {sum(1 for k, w, x in order if k != 'gzip' and x != 0)}")
    if bad: print("- OUT OF ORDER: " + "; ".join(bad))
    t_start = None
    if a.start: t_start = int(DT.datetime.strptime(a.start, "%Y-%m-%d %H:%M:%S").replace(tzinfo=UTC).timestamp())
    R = lambda k: sum(r[k] for r in rows)
    ph = [("scan (the job's checks and the scan; to the plan)", t_start, t_plan), ("moves (to .merge_complete)", t_plan, t_mc),
          ("reduce_B (reduce_B.py, 16 cells)" if len(rows) == 16 else f"reduce_B ({len(rows)} cells)", None, R("rB")),
          (f"reduce_s1 (reduce_s1_koa.py, {R('n')} folders)", None, R("rs")), ("sha (the SHA-256 lines of the traces)", None, R("sha")),
          ("gzip -9 of the traces", None, R("gz")), ("after the last folder (to the job log's copy)", None, tail)]
    total = (t_end - t_start) if t_start is not None else None
    den = total if total else (t_end - t_plan)                    # without the start: the shares of the span from the plan to the end
    print(f"\n| phase | from | to | seconds | h:mm:ss | share of {'the job' if total else 'the span from the plan to the end'} |\n|---|---|---|---|---|---|")
    for lab, x, y in ph:
        if lab.startswith(("scan", "moves")):
            d = (y - x) if x is not None else None
            fr = ts(x)[11:] if x is not None else "unknown"; to = ts(y)[11:]
        else:
            d, fr, to = y, "", ""
        sh = f"{100 * d / den:.1f} %" if (d is not None and den and (total or not lab.startswith("scan"))) else "-"
        print(f"| {lab} | {fr} | {to} | {d if d is not None else '-'} | {hms(d) if d is not None else '-'} | {sh} |")
    print(f"| **the job** | {ts(t_start)[11:] if t_start is not None else 'unknown'} | {ts(t_end)[11:]} | {total if total is not None else '-'} | "
          f"{hms(total) if total is not None else '-'} | |")
    print(f"\nfrom the plan to the end: {t_end - t_plan} s ({hms(t_end - t_plan)}); the job's start "
          + (f"{a.start} UTC (given)" if a.start else "is in no fetched file (sacct gives it)") + f"; the end {ts(t_end)} UTC (the merge log's copy)")
    print(f"\nper cell (seconds):\n\n| cell | m_<M> folders | reduce_B | reduce_s1 | sha | gzip | total |\n|---|---|---|---|---|---|---|")
    for r in rows: print(f"| {r['cid']} | {r['n']} | {r['rB']} | {r['rs']} | {r['sha']} | {r['gz']} | {r['tot']} |")
    print(f"| all | {R('n')} | {R('rB')} | {R('rs')} | {R('sha')} | {R('gz')} | {R('tot')} |")


if __name__ == "__main__":
    main()
