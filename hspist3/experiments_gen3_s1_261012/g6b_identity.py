#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 b; 261012 sec. 4.7.81): THE PROOF FOR GATE G6 (b). The M5 job on the a643bc7
double and long-double pair is the gate PROVIDED every source compiled into the binary except 00ALLINONE.c is byte-identical between
a643bc7 and f72d764, and the parts of 00ALLINONE.c that differ are quoted. Read-only (git show / git diff of the two commits).
  1. The compiled sources: the Makefile's link line ($(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS)) at each commit, and
     every file they #include "..." (recursively, resolved next to the including file and in hspist3/); kissfft is a submodule, so
     its files are read from its own repository at the gitlink commit of each tree.
  2. Per file: the blob's SHA-256 at a643bc7 and at f72d764; git diff --stat over the same paths.
  3. 00ALLINONE.c: every hunk of the diff with its line ranges, its enclosing function (git's hunk header) and the ##CHRIS tag that
     names its change; the full diff is written to g6b_00ALLINONE_diff.txt.
usage (from hspist3/): python3 experiments_gen3_s1_261012/g6b_identity.py > experiments_gen3_s1_261012/g6b_identity_output.txt"""
import hashlib, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); ROOT = os.path.dirname(HS)
A, F = "a643bc7", "f72d764"


def git(*a, cwd=ROOT, binary=False):
    r = subprocess.run(["git", "-C", cwd] + list(a), capture_output=True, check=True)
    return r.stdout if binary else r.stdout.decode("utf-8", "replace")


def blob(commit, path):
    """bytes of hspist3/<path> at commit; kissfft/* from the submodule's own repository at the tree's gitlink"""
    if path.startswith("kissfft/"):
        link = git("ls-tree", commit, "hspist3/kissfft").split()[2]
        return git("show", f"{link}:{path[len('kissfft/'):]}", cwd=os.path.join(HS, "kissfft"), binary=True), link
    return git("show", f"{commit}:hspist3/{path}", binary=True), None


def exists(commit, path):
    try: blob(commit, path); return True
    except subprocess.CalledProcessError: return False


def sources(commit):
    mk = git("show", f"{commit}:hspist3/Makefile")
    var = lambda n: re.search(rf"^{n}\s*\??:?=\s*(.*)$", mk, re.M).group(1).split()
    link = re.search(r"^\t\$\(CC\) \$\(LDFLAGS\) -o \$@ (.*)$", mk, re.M).group(1)
    assert "$(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS)" in link, link
    top = [var("TG")[0] + ".c"] + var("KISSFFT_SRCS") + var("EDMD_SRCS") + var("VALIDATION_SRCS")
    seen, todo = [], list(top)
    while todo:
        p = todo.pop(0)
        if p in seen: continue
        seen.append(p)
        text = blob(commit, p)[0].decode("utf-8", "replace")
        for inc in re.findall(r'^\s*#\s*include\s+"([^"]+)"', text, re.M):
            for cand in (os.path.normpath(os.path.join(os.path.dirname(p), inc)), os.path.normpath(inc)):
                if exists(commit, cand): todo.append(cand); break
            else:
                print(f"   (#include \"{inc}\" in {p} at {commit}: not in the tree; a system header or generated)")
    return top, seen


def main():
    print(f"# Gate G6 (b): the sources compiled into 00ALLINONE at {A} and {F} (decision 39 item 5 b), printed by experiments_gen3_s1_261012/g6b_identity.py\n")
    ta, sa = sources(A); tf, sf = sources(F)
    print(f"the Makefile's link line, both commits: $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS) -> {' '.join(tf)}")
    print(f"compiled sources and their #include closure: {A} {len(sa)} files, {F} {len(sf)} files; the same list: {sorted(sa) == sorted(sf)}\n")
    allp = sorted(set(sa) | set(sf) | {"Makefile"})
    print("| file | SHA-256 at a643bc7 | SHA-256 at f72d764 | identical |\n|---|---|---|---|")
    differ = []
    for p in allp:
        ba, la = blob(A, p) if exists(A, p) else (None, None); bf, lf = blob(F, p) if exists(F, p) else (None, None)
        ha = hashlib.sha256(ba).hexdigest() if ba is not None else "absent"; hf = hashlib.sha256(bf).hexdigest() if bf is not None else "absent"
        same = ha == hf
        if not same: differ.append(p)
        print(f"| {p}{' (submodule ' + la[:8] + ')' if la else ''} | {ha} | {hf} | {'yes' if same else '**NO**'} |")
    print(f"\nfiles that differ: {', '.join(differ) if differ else 'none'}")
    print(f"kissfft gitlink: {A} {git('ls-tree', A, 'hspist3/kissfft').split()[2]}, {F} {git('ls-tree', F, 'hspist3/kissfft').split()[2]}")
    paths = ["hspist3/" + p for p in allp if not p.startswith("kissfft/")] + ["hspist3/kissfft"]
    print(f"\ngit diff --stat {A} {F} -- <the same paths>:\n")
    print("    " + git("diff", "--stat", A, F, "--", *paths).rstrip().replace("\n", "\n    "))
    ok = differ == ["00ALLINONE.c"]
    print(f"\nEVERY COMPILED SOURCE EXCEPT 00ALLINONE.c IS BYTE-IDENTICAL BETWEEN {A} AND {F}: {'YES' if ok else 'NO'}")
    print("(the build also passes -DBUILD_GIT=<commit> (Makefile PROVENANCE): the version line, \"git a643bc7\" or \"git f72d764\", differs by design)")
    # ------------------------------------------------------------ the parts of 00ALLINONE.c that differ
    d = git("diff", "-U0", A, F, "--", "hspist3/00ALLINONE.c")
    open(os.path.join(HERE, "g6b_00ALLINONE_diff.txt"), "w").write(git("diff", A, F, "--", "hspist3/00ALLINONE.c"))
    hunks = re.split(r"(?m)^(?=@@ )", d)[1:]
    print(f"\n## The parts of 00ALLINONE.c that differ: {len(hunks)} hunks (git diff -U0 {A} {F}; the full diff with context in g6b_00ALLINONE_diff.txt)\n")
    print(f"| # | {A} lines | {F} lines | - / + | enclosing function (git's hunk header) | the change, by its ##CHRIS tag or first changed line |\n|---|---|---|---|---|---|")
    for k, h in enumerate(hunks, 1):
        m = re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@ ?(.*)", h)
        a0, an, f0, fn, ctx = int(m.group(1)), int(m.group(2) or 1), int(m.group(3)), int(m.group(4) or 1), m.group(5).strip()
        body = h.split("\n")[1:]
        minus = [l[1:] for l in body if l.startswith("-")]; plus = [l[1:] for l in body if l.startswith("+")]
        tag = next((re.sub(r"\s+", " ", l.strip())[:150] for l in plus + minus if "##CHRIS" in l), None)
        first = next((re.sub(r"\s+", " ", l.strip())[:150] for l in plus + minus if l.strip()), "")
        what = tag or first
        print(f"| {k} | {a0}-{a0 + an - 1 if an else a0} ({an}) | {f0}-{f0 + fn - 1 if fn else f0} ({fn}) | {len(minus)} / {len(plus)} | "
              f"{ctx[:70].replace('|', chr(92) + '|')} | {what.replace('|', chr(92) + '|')} |")   # a pipe in the quoted code is written \| (markdown)
    st = git("diff", "--numstat", A, F, "--", "hspist3/00ALLINONE.c").split()
    print(f"\n00ALLINONE.c: +{st[0]} -{st[1]} lines; the commits: " +
          "; ".join(l[:110] for l in git("log", "--format=%h %s", f"{A}..{F}", "--", "hspist3/00ALLINONE.c").strip().split("\n")))


if __name__ == "__main__":
    main()
