#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.33; plan-author decision 15, item 1b): EVERY FILE THE KOA STEPS OF THE RUNSHEET EXECUTE OR READ,
with its git blob in the gate commit a643bc7 (what KOA builds) and in the runsheet's commit f8e763e, and whether the runsheet's sparse
checkout puts it on KOA. The list is traced, not typed, from the entry points the runsheet names step by step:
  - shell scripts (.sh, .sbatch): the repository paths their code lines name (comment lines dropped), "$HERE"/"$CF" paths, `make`
    targets (the Makefile and its sources);
  - Python scripts: every local import (ast, at any depth), resolved through the sys.path inserts the scripts make (evaluated from their
    own os.path expressions, in import order), recursively;
  - C: the Makefile's source lists and the harness builds, with the closure of their #include "..." lines;
  - data the scripts read by a fixed name (task lists, SHA-256 files), from the shell lines and a short list of Python reads stated below.
Each file is printed with its step, where it runs (KOA or the Mac) and its role. A path a code line names that is missing in a643bc7 is
printed with the reason it can or cannot matter. Read-only (git objects only).
usage (from anywhere in the repository): python3 hspist3/cluster/gen3_koa_261009/koa_files_check.py
"""
import ast, os, re, subprocess, sys

GATE_C, RUN_C = "a643bc7edb1b41574c0f9af0c5fea199dc797680", "f8e763e4e9e449a762d71bfb452828f88e4ae9a1"
SPARSE = ("/hspist3/*.c", "/hspist3/*.h", "/hspist3/*.py", "/hspist3/Makefile", "/hspist3/edmd_core/", "/hspist3/cluster/", "/hspist3/validation/",
          "/hspist3/kissfft")
K = "hspist3/cluster/gen3_koa_261009/"
ENTRY = [   # (step, where, role, path): the files the runsheet itself names
    ("2", "KOA", "executed (bash)", K + "build_gen3_koa.sh"),
    ("3", "KOA", "executed (sbatch)", K + "g3_gate_koa.sbatch"),
    ("4", "KOA", "executed (sbatch)", K + "g3_det_xnode.sbatch"),
    ("4", "KOA", "executed (sbatch)", K + "g3_xnode_compare.sbatch"),
    ("5", "KOA", "executed (sbatch)", K + "testG_koa.sbatch"),
    ("5", "KOA", "executed (sbatch)", K + "merge_testG_koa.sbatch"),
    ("ext", "Mac", "executed (python3; only if the verdict asks for an extension)", K + "gen_testG_koa_261009.py"),
    ("6", "Mac", "executed (bash)", K + "fetch_gen3_koa.sh"),
    ("6", "Mac", "executed (python3)", K + "testG_koa_verdict_261009.py"),
    ("all", "Mac", "read by Chris", K + "runsheet_gen3_koa_261009.txt"),
]
PY_READS = {   # files Python code opens by a fixed name (not visible to the import trace)
    "hspist3/cluster/gen3_gate_261009/gen_testG_261009.py": ["hspist3/cluster/confinement_20261013/tasks_B_epi8_H_H10_L10.txt"],
    K + "testG_koa_verdict_261009.py": [K + "tasks_testG_koa_261009.txt"],
}
NOT_TAKEN = {   # paths named on a code path the KOA steps never take
    K + "af_stream_reduce.py": "testG_worker_koa.sh calls it only on the FIFO path (KEEP_EV unset); testG_koa.sbatch sets KEEP_EV=1 for every "
                               "AF line, so KOA never runs it (the Mac's copy is cluster/gen3_gate_261009/af_stream_reduce.py)",
}
THIRD = {"numpy", "pandas", "scipy", "matplotlib"}
ROOT = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip()


def git(*a): return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True)
def blob(c, p): r = git("rev-parse", f"{c}:{p}"); return r.stdout.strip() if r.returncode == 0 else None
def show(c, p): r = git("show", f"{c}:{p}"); return r.stdout if r.returncode == 0 else None
TREE = set(git("ls-tree", "-r", "--name-only", GATE_C, "--", "hspist3").stdout.split("\n")) - {""}


def in_sparse(p):
    q = "/" + p
    for pat in SPARSE:
        if pat.endswith("/"):
            if q.startswith(pat): return True
        elif "*" in pat:
            d, ext = pat.rsplit("/", 1)[0], pat.rsplit("*", 1)[1]
            if q.startswith(d + "/") and "/" not in q[len(d) + 1:] and q.endswith(ext): return True
        elif q == pat or q.startswith(pat + "/"): return True
    return False


def code_lines(text):
    out = []
    for l in text.split("\n"):
        s = l.strip()
        if not s or s.startswith("#"): continue
        m = re.search(r"\s#[^!]", l)            # a trailing comment, if the quotes before it are balanced
        if m and l[:m.start()].count('"') % 2 == 0 and l[:m.start()].count("'") % 2 == 0: l = l[:m.start()]
        out.append(l)
    return out


def shell_refs(path, text):
    """repository paths a shell script's code lines name: (path, role)"""
    here, hs = os.path.dirname(path), "hspist3"
    refs = []
    for l in code_lines(text):
        for m in re.finditer(r"(?<![\w/$.-])((?:cluster|validation|edmd_core)/[A-Za-z0-9_./-]*[A-Za-z0-9_])", l):
            p = hs + "/" + m.group(1)
            pre = l[:m.start()]
            role = ("sourced" if re.search(r"(^|\s)(source|\.)\s+$", pre) else "executed (python3)" if re.search(r"python3\s+$", pre)
                    else "executed (bash)" if re.search(r"bash\s+$", pre) else "compiled" if p.endswith((".c", ".h")) else "read")
            refs.append((p, role))
        for m in re.finditer(r"\$(HERE|CF)/([A-Za-z0-9_.-]+)", l):
            base = here if m.group(1) == "HERE" else hs + "/cluster/confinement_20261013"
            pre = l[:m.start()].rstrip('"')
            refs.append((base + "/" + m.group(2), "executed (python3)" if pre.rstrip().endswith("python3") else "read"))
        if re.search(r"\bmake\b.*\bkoa", l):
            refs.append(("hspist3/Makefile", "read (make)"))
    return refs


def makefile_sources(c):
    t = show(c, "hspist3/Makefile") or ""
    var = {m.group(1): m.group(2).split() for m in re.finditer(r"^(\w+)\s*[?:]?=\s*(.*)$", t, re.M)}
    rule = re.search(r"^\$\(TG\):\s*(.*)$", t, re.M)
    def expand(tok, depth=0):   # every $(VAR) inside a token, recursively (a list variable gives several tokens)
        m = re.search(r"\$\((\w+)\)", tok)
        if not m or depth > 5: return [tok]
        return [y for v in (var.get(m.group(1)) or [""]) for y in expand(tok[:m.start()] + v + tok[m.end():], depth + 1)]
    srcs = [y for tok in (rule.group(1).split() if rule else []) for y in expand(tok)]
    return ["hspist3/" + s for s in srcs if s.endswith((".c", ".h"))]


def includes(path, text):
    out = []
    for m in re.finditer(r'^\s*#\s*include\s+"([^"]+)"', text, re.M):
        out.append(os.path.normpath(os.path.join(os.path.dirname(path), m.group(1))))
    return out


def py_paths(path, text):
    """the sys.path inserts a Python file makes at its top level, evaluated from its own os.path expressions"""
    env, ins = {"__file__": "/" + path}, []
    def ev(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, str): return n.value
        if isinstance(n, ast.Name): return env.get(n.id)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute):
            f = ast.unparse(n.func); a = [ev(x) for x in n.args]
            if any(x is None for x in a): return None
            if f == "os.path.dirname": return os.path.dirname(a[0])
            if f == "os.path.abspath": return a[0]
            if f == "os.path.join": return os.path.join(*a)
        return None
    try: tree = ast.parse(text)
    except SyntaxError: return []
    for st in tree.body:
        if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name):
            v = ev(st.value)
            if v is not None: env[st.targets[0].id] = v
        for node in ast.walk(st) if isinstance(st, ast.Expr) else []:
            if isinstance(node, ast.Call) and ast.unparse(node.func) == "sys.path.insert" and len(node.args) == 2:
                v = ev(node.args[1])
                if v: ins.append(os.path.normpath(v).lstrip("/"))
    return ins


def py_imports(text):
    try: tree = ast.parse(text)
    except SyntaxError: return []
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import): out += [a.name.split(".")[0] for a in n.names]
        elif isinstance(n, ast.ImportFrom) and n.module and n.level == 0: out.append(n.module.split(".")[0])
    return out


def main():
    print("# The files the KOA runsheet's steps execute or read: gate commit a643bc7 against the runsheet's commit f8e763e (261012 sec. 4.7.33),"
          " printed by hspist3/cluster/gen3_koa_261009/koa_files_check.py\n")
    seen = {}   # path -> [(step, where, role, via)]
    missing = []
    def add(p, step, where, role, via):
        seen.setdefault(p, []).append((step, where, role, via))
    stdlib = set(getattr(sys, "stdlib_module_names", ())) | THIRD
    for step, where, role, p in ENTRY:
        add(p, step, where, role, "the runsheet")
        stack = [(p, [os.path.dirname(p)])]
        done = set()
        while stack:
            f, spath = stack.pop()
            if f in done: continue
            done.add(f)
            text = show(GATE_C, f)
            if text is None: continue
            if f.endswith((".sh", ".sbatch")):
                for q, r in shell_refs(f, text):
                    if q not in TREE and any(t.startswith(q + "/") for t in TREE): continue     # a directory (a data root, logs/)
                    add(q, step, where, r, f)
                    if q not in TREE and q not in NOT_TAKEN: missing.append((q, f))
                    if q.endswith((".sh", ".sbatch", ".py", ".c", ".h")): stack.append((q, [os.path.dirname(q)]))
                    if q == "hspist3/Makefile":
                        for s in makefile_sources(GATE_C): add(s, step, where, "compiled (make)", q); stack.append((s, []))
            elif f.endswith(".py"):
                sp = list(spath)
                for d in py_paths(f, text): sp.insert(0, d)
                for mod in py_imports(text):
                    if mod in stdlib: continue
                    hit = next((d + "/" + mod + ".py" for d in sp if d + "/" + mod + ".py" in TREE), None)
                    if hit is None:
                        cands = sorted(t for t in TREE if t.endswith("/" + mod + ".py"))
                        add(f"(import {mod})", step, where, "import NOT RESOLVED in a643bc7" + (f"; by name: {cands}" if cands else ""), f)
                        missing.append((f"import {mod}", f)); continue
                    add(hit, step, where, "imported", f); stack.append((hit, sp))
                for q in PY_READS.get(f, []): add(q, step, where, "read (by name)", f)
            elif f.endswith((".c", ".h")):
                for q in includes(f, text):
                    if q in TREE: add(q, step, where, "compiled (#include)", f); stack.append((q, []))
                    elif not q.startswith("hspist3/kissfft/"): missing.append((q, f))
    print("| file | step | runs on | role (first reference) | blob in a643bc7 | blob in f8e763e | same | on KOA (sparse checkout) |\n|---|---|---|---|---|---|---|---|")
    differ, absent, outside = [], [], []
    for p in sorted(seen, key=lambda p: (min((s for s, *_ in seen[p]), key=lambda s: (s not in "23456", s)), p)):
        steps = sorted({s for s, *_ in seen[p]}); wh = sorted({w for _, w, *_ in seen[p]}); role = seen[p][0][2]
        if p.startswith("("): print(f"| {p} | {', '.join(steps)} | {', '.join(wh)} | **{role}** (from {seen[p][0][3]}) | | | | |"); continue
        a, b = blob(GATE_C, p), blob(RUN_C, p)
        if p.startswith("hspist3/kissfft/"): a = b = "(the separate kissfft clone, febd4ca)"
        same = "yes" if a == b else ("**NO**" if a and b else "")
        koa = "yes" if in_sparse(p) else "**no**"
        note = f" -- {NOT_TAKEN[p]}" if p in NOT_TAKEN else ""
        if a is None: absent.append((p, wh)); a = "**missing**"
        elif same == "**NO**": differ.append((p, wh))
        if "KOA" in wh and koa != "yes" and not p.startswith("hspist3/kissfft/"): outside.append(p)
        print(f"| {p} | {', '.join(steps)} | {', '.join(wh)} | {role}{note} | {a} | {b or '**missing**'} | {same} | {koa} |")
    print(f"\nthe gitlink hspist3/kissfft: a643bc7 {git('ls-tree', GATE_C, 'hspist3/kissfft').stdout.split()[2] if git('ls-tree', GATE_C, 'hspist3/kissfft').stdout else '-'}; "
          f"f8e763e {git('ls-tree', RUN_C, 'hspist3/kissfft').stdout.split()[2] if git('ls-tree', RUN_C, 'hspist3/kissfft').stdout else '-'} (step 1 clones kissfft "
          "separately at febd4caeed32e33ad8b2e0bb5ea77542c40f18ec)")
    koa_diff = [p for p, wh in differ if "KOA" in wh]; koa_abs = [p for p, wh in absent if "KOA" in wh and p not in NOT_TAKEN]
    print(f"\nfiles traced: {len([p for p in seen if not p.startswith('(')])}; differing between a643bc7 and f8e763e: {len(differ)} "
          f"({'; '.join(f'{p} [{chr(44).join(w)}]' for p, w in differ) or 'none'}); missing in a643bc7: {len(absent)} "
          f"({'; '.join(f'{p} [{chr(44).join(w)}]' for p, w in absent) or 'none'}); used on KOA but outside the sparse checkout: {outside or 'none'}; "
          f"unresolved references: {[m for m in missing if not m[0].startswith(K + 'af_stream')] or 'none'}")
    ok = not koa_diff and not koa_abs and not outside and not [m for m in missing if not m[0].startswith(K + "af_stream")]
    print(f"VERDICT (information): every file the KOA steps execute or read is present in a643bc7, identical in f8e763e and inside the sparse checkout: "
          f"{'YES -- KOA builds a643bc7 as the runsheet says; the files that differ are read or run on the Mac only' if ok else '**NO** -- see the rows above'}")


if __name__ == "__main__":
    main()
