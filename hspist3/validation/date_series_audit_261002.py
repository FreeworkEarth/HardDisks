#!/usr/bin/env python3
"""##CHRIS 2026-10-02 (Task I): the "series" dates written into file names, ##CHRIS headers and STATUS timestamps,
against the git commit dates, for every commit since 2026-09-20. Git commit timestamps are authoritative; nothing is
renamed. From 2026-10-02 on, new names and STATUS lines use the machine date.

Per commit, the series date is the LATEST date written into what the commit ADDED:
  - STATUS timestamps   `- **YYYY-MM-DD HH:MM:SS HST**` in added lines of 260913_tests_STATUS.md
  - code/doc headers    `##CHRIS YYYY-MM-DD` in added lines of text files (.py .sh .md .tex .c .h .sbatch)
  - file names          `26MMDD_...` or `..._2026MMDD...` in the names of files the commit added
offset = series date - commit date (calendar days, commit date in HST = the repo's -1000).
usage: python3 hspist3/validation/date_series_audit_261002.py
"""
import datetime as dt, os, re, subprocess
REPO = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip()
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True, cwd=REPO, errors="replace").stdout
TEXT = ["*.py", "*.sh", "*.md", "*.tex", "*.c", "*.h", "*.sbatch"]
RX_STATUS = re.compile(r"^\+- \*\*(20\d\d-\d\d-\d\d) \d\d:\d\d:\d\d")
RX_CHRIS = re.compile(r"##CHRIS (20\d\d-\d\d-\d\d)")
RX_NAME6 = re.compile(r"(?<![\d])(2[56][01]\d[0-3]\d)_")       # 26MMDD_ (and 25MMDD_)
RX_NAME8 = re.compile(r"(?<![\d])(202[56][01]\d[0-3]\d)(?![\d])")

def day(s):
    try:
        return dt.date(int(s[:4]), int(s[5:7]), int(s[8:10])) if "-" in s else \
               (dt.date(2000 + int(s[:2]), int(s[2:4]), int(s[4:6])) if len(s) == 6 else dt.date(int(s[:4]), int(s[4:6]), int(s[6:8])))
    except ValueError:
        return None

def scan(*rng):
  rows = []
  for line in git("log", *rng, "--reverse", "--date=format:%Y-%m-%d %H:%M", "--format=%h%x09%cd%x09%s").splitlines():
      h, cd, subj = line.split("\t", 2)
      cdate = dt.date.fromisoformat(cd[:10]); found = {}
      diff = git("show", "--format=", "--unified=0", "--no-color", h, "--", *TEXT)
      for l in diff.splitlines():
          if not l.startswith("+") or l.startswith("+++"): continue
          m = RX_STATUS.match(l)
          if m and (d := day(m.group(1))): found.setdefault("STATUS", set()).add(d)
          for m in RX_CHRIS.finditer(l):
              if (d := day(m.group(1))): found.setdefault("##CHRIS", set()).add(d)
      for name in git("show", "--format=", "--name-only", "--diff-filter=A", h).splitlines():
          base = os.path.basename(name)
          for rx in (RX_NAME6, RX_NAME8):
              for m in rx.finditer(base):
                  d = day(m.group(1))
                  if d and dt.date(2026, 8, 1) <= d <= dt.date(2026, 12, 31): found.setdefault("file name", set()).add(d)
      if found:
          src, sd = max(((k, max(v)) for k, v in found.items()), key=lambda kv: kv[1])
          rows.append((h, cdate, cd, sd, (sd - cdate).days, src, subj))
      else:
          rows.append((h, cdate, cd, None, None, "-", subj))
  return rows

rows = scan("--since=2026-09-20")

print("### Series date written into each commit vs its git commit date (every commit since 2026-09-20, oldest first)\n")
print("| # | hash | commit date (HST) | series date | offset [days] | latest series date found in | subject |")
print("|---|---|---|---|---|---|---|")
for i, (h, c, cd, sd, off, src, subj) in enumerate(rows, 1):
    print(f"| {i} | {h} | {cd} | {sd or '-'} | {'' if off is None else f'{off:+d}'} | {src} | {subj[:90]} |")
aff = [r for r in rows if r[4] is not None and r[4] >= 2]
print(f"\ncommits: {len(rows)}; carrying a series date: {sum(r[3] is not None for r in rows)}; "
      f"series date >= 2 days ahead of the commit: {len(aff)}")
if aff:
    print(f"first affected: {aff[0][0]} ({aff[0][2]}, series {aff[0][3]}, {aff[0][4]:+d} d); "
          f"last affected: {aff[-1][0]} ({aff[-1][2]}, series {aff[-1][3]}, {aff[-1][4]:+d} d)")
    by = {}
    for r in aff: by.setdefault(r[1], []).append(r[4])
    print("offset by commit day: " + "; ".join(f"{d}: {min(v):+d} to {max(v):+d}" for d, v in sorted(by.items())))

early = [r for r in scan("--since=2026-08-01", "--until=2026-09-20") if r[4] is not None and r[4] >= 2]
print(f"before 2026-09-20 (commits since 2026-08-01): {len(early)} with a series date >= 2 days ahead; "
      + (f"earliest {early[0][0]} ({early[0][2]}, series {early[0][3]}, {early[0][4]:+d} d)" if early else "none"))
