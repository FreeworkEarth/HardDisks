#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 31, item 2; 261012 sec. 4.7.61): record the stage 1 registration text
(261010_prereg_v1_stage1_fluid.md, the file Chris names) VERBATIM as sec. 4.8 of the notes, binary-safe, with no rewording. Amendments A1
(decision 28) and A2 (decision 30, POST HOC) follow it with their dates.
HOW (binary-safe). The section is appended to the end of the notes file in binary mode, as `cat >>` does: nothing above it changes, so
the notes' 15 raw CR bytes stay as they are. The file's bytes stand unchanged between two fence lines of F backticks, F being one more
than the longest run of backticks in the file (at least 3), so no line of the file can close the fence. If the file does not end with a
newline, one newline is added before the closing fence. The header line of the section says so, and --check takes it off again.
THE RECORDED TEXT, for the SHA-256: the bytes after the opening fence line's newline, up to the closing fence line, less that added
newline if there is one. They must equal the file byte for byte, and both SHA-256 are printed.
A1 and A2 are not retyped. Each is cut from the verbatim block of the decision that carries it, in the committed notes: A1 from sec.
4.7.57 (decision 28), A2 from sec. 4.7.59 (decision 30). The cut runs from the amendment's first line to its last line, as named below,
and is printed in its own fence.
usage (from the repository root):
  python3 hspist3/validation/record_stage1_prereg_261011.py --selftest
  python3 hspist3/validation/record_stage1_prereg_261011.py --file <path> --dry-run   build the section into the scratch file named by
                                                                                      --out, print the SHA-256s; the notes are untouched
  python3 hspist3/validation/record_stage1_prereg_261011.py --file <path> --append    append it to the notes
  python3 hspist3/validation/record_stage1_prereg_261011.py --file <path> --check <commit>
                                                                                      the recorded text in <commit>'s notes = the file?
"""
import argparse, hashlib, os, re, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
NOTES = "0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md"
NAME = "261010_prereg_v1_stage1_fluid.md"
BEGIN, END = f"<!-- BEGIN VERBATIM: {NAME} -->", f"<!-- END VERBATIM: {NAME} -->"
# the amendments: (label, the section that carries the decision verbatim, first line, last line), cut from the committed notes
A1 = ("4.7.57", " AMENDMENT A1 to the stage 1 registration (sec. 4.8.10, gate G6 a),",
      "  Test G and Test R are both reported wherever either is cited.\"")
A2 = ("4.7.59", "    AMENDMENT A2 to amendment A1 (POST HOC: written after C1's",
      "    is bounded separately by the contact audit (max 4.5e-12 px).")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def fence_for(b):
    runs = [len(m) for m in re.findall(rb"`+", b)]
    return "`" * max(3, (max(runs) + 1) if runs else 3)


def section_of(notes_text, num):
    m = re.search(rf"^### {re.escape(num)} .*?(?=^### |\Z)", notes_text, re.M | re.S)
    if not m: sys.exit(f"STOP: sec. {num} not found in the notes")
    return m.group(0)


def cut(notes_text, spec):
    num, first, last = spec
    lines = section_of(notes_text, num).split("\n")
    i = [k for k, l in enumerate(lines) if l == first]; j = [k for k, l in enumerate(lines) if l == last]
    if len(i) != 1 or len(j) != 1 or j[0] < i[0]: sys.exit(f"STOP: sec. {num}: the amendment's first/last line found {len(i)}/{len(j)} times")
    return "\n".join(lines[i[0]:j[0] + 1])


def build(file_bytes, path_given, notes_text, when, script_commit):
    f = fence_for(file_bytes); nl = file_bytes.endswith(b"\n")
    a1, a2 = cut(notes_text, A1), cut(notes_text, A2)
    head = (f"\n### 4.8 The stage 1 registration, verbatim: {NAME}, with amendment A1 (decision 28) and amendment A2 (decision 30, "
            f"POST HOC) (recorded {when}, machine date) [SOURCE]\n\n"
            f"**How this section was made** (plan-author decision 31, item 2; sec. 4.7.61).\n"
            f"- The file Chris named: `{path_given}`; {len(file_bytes)} bytes; SHA-256 `{sha(file_bytes)}`.\n"
            f"- Its bytes stand below unchanged, by `hspist3/validation/record_stage1_prereg_261011.py` (main {script_commit}). They are "
            f"appended in binary mode between the fence line after the comment line BEGIN VERBATIM and the fence line before the comment "
            f"line END VERBATIM; the fence is {len(f)} backticks.\n"
            f"- The file ends with a newline: {'yes' if nl else 'no (one newline was added before the closing fence; it is not part of the text)'}.\n"
            f"- To check: `python3 hspist3/validation/record_stage1_prereg_261011.py --file <the file> --check <commit>` cuts the recorded "
            f"text out of that commit's notes and compares it with the file, byte for byte and by SHA-256.\n\n"
            f"{BEGIN}\n\n{f}\n").encode("utf-8")
    tail = ((b"" if nl else b"\n") + f"{f}\n\n{END}\n\n".encode("utf-8") +
            (f"**Amendment A1 to the stage 1 registration** (plan-author decision 28, 2026-10-11 ~09:45 HST; dated in it \"2026-10-11, "
             f"before any Test R or C1 data\"; recorded in sec. 4.7.57), verbatim from that section:\n\n```\n{a1}\n```\n\n"
             f"**Amendment A2 to amendment A1 -- POST HOC** (written after C1's data, and to be flagged so wherever it is cited; plan-author "
             f"decision 30, 2026-10-11 ~11:30 HST; recorded in sec. 4.7.59), verbatim from that section:\n\n```\n{a2}\n```\n\n"
             f"(Not part of the registration, for orientation only: A2's condition was printed on C1's data as MET, POST HOC, in sec. "
             f"4.7.60. The registered C1 verdict NOT DECIDED stands. Test R's verdict is Report R1's.)\n").encode("utf-8"))
    return head + file_bytes + tail


def recorded_text(notes_bytes):
    """the bytes between the opening fence line after BEGIN and the closing fence line before END, less an added newline"""
    b = notes_bytes; bl, el = ("\n" + BEGIN + "\n").encode(), ("\n" + END + "\n").encode()
    if b.count(bl) != 1 or b.count(el) != 1: sys.exit("STOP: the BEGIN/END lines: not exactly one each")
    i = b.find(bl) + 1; j = b.find(el) + 1                      # i, j: the starts of the two comment lines
    seg = b[i + len(BEGIN):j]                                   # "\n\n" + fence + "\n" + text [+ "\n"] + fence + "\n\n"
    m = re.match(rb"\n\n(`{3,})\n", seg)
    if not m: sys.exit("STOP: no opening fence after BEGIN")
    f = m.group(1); body = seg[m.end():]
    closing = b"\n" + f + b"\n\n"
    if not body.endswith(closing): sys.exit("STOP: no closing fence before END")
    text = body[:-len(closing) + 1]                             # keeps the text's own final newline
    hdr = b[:i].rsplit(b"\n### 4.8 ", 1)[-1]
    if b"The file ends with a newline: no" in hdr: text = text[:-1]
    return text


def git_show(commit, path):
    return subprocess.run(["git", "-C", ROOT, "show", f"{commit}:{path}"], check=True, capture_output=True).stdout


def selftest():
    ok = {}
    notes = open(os.path.join(ROOT, NOTES), "rb").read().decode("utf-8")
    cases = {"plain, ends with a newline": b"# Stage 1\n\n## 4.8.10 Gate G6\n(a) text\n",
             "no final newline": b"# Stage 1\nlast line without newline",
             "backtick runs and a fence inside": b"text ```` four\n```python\nx = 1\n```\n",
             "CR bytes, tabs, UTF-8 and a BOM": "﻿# Stufe 1\r\n\tEinrückung – ± η\r\n".encode("utf-8"),
             "empty lines at the end": b"a\n\n\n"}
    with tempfile.TemporaryDirectory() as td:
        for lab, fb in cases.items():
            base = b"# notes\nline with a raw CR\r\n" * 3
            sec = build(fb, "/tmp/x.md", notes, "2026-10-11 00:00 HST", "selftest")
            got = recorded_text(base + sec)
            ok[f"{lab}: recorded text = the file (bytes and SHA-256)"] = got == fb and sha(got) == sha(fb)
            ok[f"{lab}: the notes above are unchanged"] = (base + sec).startswith(base)
    a1, a2 = cut(notes, A1), cut(notes, A2)
    ok["A1 cut from sec. 4.7.57: 6 lines, from 'AMENDMENT A1' to '...wherever either is cited.\"'"] = (
        a1.count("\n") == 5 and a1.lstrip().startswith("AMENDMENT A1") and "Test R prints ENGINES AGREE and C1 prints PASS" in a1)
    ok["A2 cut from sec. 4.7.59: 12 lines (its reasons included), from 'AMENDMENT A2' to '...(max 4.5e-12 px).'"] = (
        a2.count("\n") == 11 and a2.lstrip().startswith("AMENDMENT A2") and "(|dt| > 1e-9 only) is <= 1e-10" in a2)
    print("| check | result |\n|---|---|")
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else 'FAIL'} |")
    n = sum(ok.values()); print(f"\nSELF-TEST: {'PASS' if n == len(ok) else 'FAIL'} ({n} of {len(ok)})")
    return n == len(ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--file"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out"); ap.add_argument("--append", action="store_true"); ap.add_argument("--check")
    a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    if not a.file: sys.exit("--file is required")
    fb = open(a.file, "rb").read()
    print(f"the file given: {a.file}\n  bytes: {len(fb)}\n  SHA-256: {sha(fb)}\n  ends with a newline: {fb.endswith(b'\\n')}; "
          f"CR bytes: {fb.count(b'\\r')}; longest backtick run: {max([len(m) for m in re.findall(rb'`+', fb)] or [0])}")
    if a.check:
        rt = recorded_text(git_show(a.check, NOTES))
        print(f"the recorded text in {a.check}:{NOTES}\n  bytes: {len(rt)}\n  SHA-256: {sha(rt)}\n"
              f"  identical to the file, byte for byte: {'YES' if rt == fb else 'NO'}")
        sys.exit(0 if rt == fb else 1)
    me = os.path.relpath(os.path.abspath(__file__), ROOT)
    if subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", me], capture_output=True, text=True).stdout.strip():
        sys.exit("STOP: commit this script first")
    commit = subprocess.run(["git", "-C", ROOT, "log", "-1", "--format=%h", "--", me], capture_output=True, text=True).stdout.strip()
    notes_path = os.path.join(ROOT, NOTES); nb = open(notes_path, "rb").read()
    if b"\n### 4.8 " in nb: sys.exit("STOP: the notes already hold a sec. 4.8")
    sec = build(fb, os.path.abspath(a.file), nb.decode("utf-8"), time.strftime("%Y-%m-%d %H:%M %Z"), commit)
    if a.dry_run:
        out = a.out or sys.exit("--dry-run needs --out")
        open(out, "wb").write(nb + sec); rt = recorded_text(nb + sec)
    elif a.append:
        with open(notes_path, "ab") as h: h.write(sec)
        rt = recorded_text(open(notes_path, "rb").read())
    else:
        sys.exit("one of --dry-run, --append, --check")
    print(f"the recorded text ({'dry run: ' + a.out if a.dry_run else NOTES})\n  bytes: {len(rt)}\n  SHA-256: {sha(rt)}\n"
          f"  identical to the file, byte for byte: {'YES' if rt == fb else 'NO'}\n  the section appended: {len(sec)} bytes")
    sys.exit(0 if rt == fb else 1)


if __name__ == "__main__":
    main()
