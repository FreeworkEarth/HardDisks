#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.1, item 0): PROVENANCE of the accelerated EDMD backend. Does any recorded run, and so any
result, rest on edmd_core/edmd_accelerated.c (selected only by the command-line flag --edmd-acc, 00ALLINONE.c:1229 / :4608-4620;
no other route sets g_edmd_is_acc)?

PART A -- the data trees. Scans the whole working tree (tracked and untracked data; not .git, not kissfft) for
  run records   00_COMMAND*.md, 01_PLOT_COMMANDS.md, command.txt, *.command.txt, run_params.json
  run logs      run.log, run_*.log, stdout.log, every *.log (the binary prints 'EDMD backend: accelerated|default' when not --quiet)
  summaries     every *.csv whose header has a 'command' column (the energy-transfer and failure summaries carry the full command)
  launchers     *.sh, *.sbatch, *.py, *.command (scripts that start runs)
  notes         *.md, *.tex, *.txt not in the classes above (mentions only)
and classifies every occurrence of '--edmd-acc' by its value (0/false/no/off = the default backend, stated explicitly;
1/true/yes/on or no value = ACCELERATED) and every 'EDMD backend: accelerated' log line as ACCELERATED. Then, per campaign, which
analysis script, notebook or LaTeX source names the campaign's directory.

PART B -- the figure provenance of the two drafts (paper1_draft.tex, paper2_draft.tex), traced from the paper side:
  1. every \\includegraphics -> the scripts that name the figure's file stem (*.py, *.sh outside the data trees);
  2. those scripts plus every script a '% TODO-source' or '% FIGURES' comment of a draft names, and their local imports
     (recursively) = the paper's script closure;
  3. every string constant of the closure (f-string constant parts included) is tested against the directory names that occur
     ONLY on paths of accelerated runs ("distinctive" names): a script that names one could read accelerated data;
  4. every directory name the closure's constants name that exists in the data trees is a data root: its run records and logs
     are counted by backend; a root that holds accelerated runs (a shared container) is resolved by testing the names of its
     accelerated children against every wildcard component of the closure's constants;
  5. every directory-listing call of the closure (glob, listdir, walk, iterdir, rglob, scandir) is printed; those with a bare
     '*' or '**' argument are listed for a check by eye of their base;
  6. figures with no producing script on disk are listed with what they are and the data trees they come from (counted).
usage (from hspist3/): python3 validation/provenance_edmd_acc_261009.py
output of 2026-10-08 (HST): hspist3/experiments_gen3_design_261009/provenance_edmd_acc_261009_output.txt
"""
import ast, fnmatch, os, re, sys
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); REPO = os.path.dirname(HS)
SKIP_DIRS = {".git", "kissfft", "node_modules", ".venv", "__pycache__"}
FLAG = re.compile(r"--edmd-acc(?:=(\S+)|\s+(?=-)|\s*$|\s+(?!-)(\S+))?")
BACKEND = re.compile(r"EDMD backend:\s*(accelerated|default)")
DIAG = "experiments_gen3_design_261009/eventcounts/acc_"          # the sec. 4.7 diagnostic runs (--edmd-acc=1 on purpose)
SELF = ("provenance_edmd_acc_261009.py", "gen3_design_numbers_261009.py")
OUTNAME = "provenance_edmd_acc_261009_output.txt"                 # this script's saved output: not scanned (stable re-runs)
DRAFTS = ("0000_PLAN_OVERALL/paper1_speedofsound/writeup/paper1_draft.tex",
          "0000_PLAN_OVERALL/paper2_energytransfer/writeup/paper2_draft.tex")
# figures of the drafts with no producing script on disk: what they are (with the source of that statement) and their data trees
NO_SCRIPT = {
    "260922_apparatus_paper": ("the apparatus drawn by the simulation itself (caption, paper1_draft.tex:80-81): a picture of the "
                               "geometry, no measured quantity", []),
    "260923_level3_settled_comparison": ("made for the 2026-09-23 figure pack (260913_tests_STATUS.md, 2026-09-22 12:34:44 HST; "
                                         "commit ba83144) from the Level-3 cells; no script on disk", ["hspist3/experiments_energy_transfer/level3_*"]),
}
LISTING = re.compile(r"\b(glob|iglob|listdir|scandir|walk|iterdir|rglob)\s*\(")
GENERIC = re.compile(r"(EDMD|mode\d_.*|plots|analysis.*|raw_simulations|eta_[0-9p.]+|m_\d+|N\d+|r\d+|run\d*|\.run\d*|logs?|final|archive"
                     r"|experiments|data|figures?|\.mplconfig|legacy|minimal|c\d+)$")


def klass(name):
    n = name.lower()
    if n.startswith("00_command") or n == "01_plot_commands.md" or n == "command.txt" or n.endswith(".command.txt") or n == "run_params.json":
        return "run record"
    if n.endswith(".log"): return "run log"
    if n.endswith(".csv"): return "summary csv"
    if n.endswith((".sh", ".sbatch", ".py", ".command")): return "launcher"
    if n.endswith((".md", ".tex", ".txt")): return "notes"
    return None


def value_kind(v):
    if v is None or v == "": return "ACCELERATED"
    v = v.strip("'\",;)`").lower()
    if v in ("0", "false", "no", "off"): return "default (explicit)"
    if v in ("1", "true", "yes", "on"): return "ACCELERATED"
    return f"other ({v})"


def scan_text(path, cls):
    """Occurrences (kind, evidence) in one file. CSVs: only files with a 'command' header column."""
    hits = []
    try:
        size = os.path.getsize(path)
        with open(path, errors="ignore") as fh:
            if cls == "summary csv":
                head = fh.readline()
                if "command" not in head.lower(): return None
                lines = [head] + fh.readlines() if size < 200e6 else [head]
            else:
                lines = fh.readlines() if size < 200e6 else []
    except (OSError, UnicodeError):
        return None
    for l in lines:
        if "edmd-acc" in l:
            for m in FLAG.finditer(l):
                v = m.group(1) if m.group(1) is not None else m.group(2)
                hits.append((value_kind(v), l.strip()[:160]))
        if "EDMD backend" in l:
            m = BACKEND.search(l)
            if m: hits.append(("ACCELERATED" if m.group(1) == "accelerated" else "default (log line)", l.strip()[:160]))
        if cls in ("launcher", "notes") and "edmd_acc" in l and "edmd-acc" not in l:
            hits.append(("code/doc mention", l.strip()[:160]))
    return hits


def scan_tree():
    """Every classified file of the working tree: list of (relpath, class, hits)."""
    rows = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for f in files:
            cls = klass(f)
            if cls is None or f == OUTNAME: continue
            p = os.path.join(root, f); h = scan_text(p, cls)
            if h is None: continue
            rows.append((os.path.relpath(p, REPO), cls, h))
    return rows


def sources(exts):
    """Scripts / sources outside the data trees (hspist3/experiments_*): list of (relpath, text)."""
    out = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith("experiments_")]
        for f in files:
            if f.endswith(exts) and f not in SELF:
                try: out.append((os.path.relpath(os.path.join(root, f), REPO), open(os.path.join(root, f), errors="ignore").read()))
                except OSError: pass
    return out


def campaign(p):
    parts = p.split(os.sep)
    for k, q in enumerate(parts):
        if q.startswith("experiments_"):
            tail = parts[k + 1:-1]
            keep = [t for t in tail if not re.fullmatch(r"(raw_simulations|eta_0p\d+|m_\d+|N\d+|\.run\d*|plots|analysis.*)", t)]
            if q == "experiments_speed_of_sound": keep = keep[2:] if len(keep) > 2 else keep    # EDMD/mode*/...
            if not keep: return p          # a summary file directly in the experiments root: the file is its own key
            return os.sep.join(parts[:k + 1] + keep[:2])
    return os.path.dirname(p)


def only_in_docstring(rel, needle):
    """True if every occurrence of needle in the Python file rel is inside its module docstring."""
    if not rel.endswith(".py"): return False
    txt = open(os.path.join(REPO, rel), errors="ignore").read()
    try: doc = ast.get_docstring(ast.parse(txt), clean=False) or ""
    except SyntaxError: return False
    return txt.count(needle) == doc.count(needle) > 0


def part_a(rows, readers_src, closure):
    scanned = Counter(); with_hits = Counter(); kinds = defaultdict(Counter); acc = []
    for rel, cls, h in rows:
        scanned[cls] += 1
        if h:
            with_hits[cls] += 1
            for k, ev in h:
                kinds[cls][k] += 1
                if k == "ACCELERATED": acc.append((rel, cls, ev))
    print("# Provenance of the accelerated EDMD backend (261012 sec. 4.7.1, item 0)\n")
    print("## PART A -- the data trees\n")
    print("| file class | files scanned | files mentioning the flag/backend | occurrences by kind |\n|---|---|---|---|")
    for cls in ("run record", "run log", "summary csv", "launcher", "notes"):
        k = ", ".join(f"{a}: {b}" for a, b in sorted(kinds[cls].items())) or "-"
        print(f"| {cls} | {scanned[cls]} | {with_hits[cls]} | {k} |")
    runs = [a for a in acc if a[1] in ("run record", "run log", "summary csv")]
    camps = defaultdict(list)
    for p, cls, ev in runs: camps[campaign(p)].append((p, cls, ev))
    notes_src = []
    for root, dirs, files in os.walk(os.path.join(REPO, "0000_PLAN_OVERALL")):
        for f in files:
            if f.endswith(".md"):
                try: notes_src.append((os.path.relpath(os.path.join(root, f), REPO), open(os.path.join(root, f), errors="ignore").read()))
                except OSError: pass
    print("\n## Campaigns whose recorded runs used the ACCELERATED backend, and who reads them\n")
    print("(readers: every *.py, *.sh, *.tex and *.ipynb outside the data trees that names the campaign's directory; notes: 0000_PLAN_OVERALL/*.md)\n")
    print("| campaign directory | files with the flag | example evidence | read by scripts / notebooks / LaTeX | named in notes |\n|---|---|---|---|---|")
    flagged = []
    for c in sorted(camps):
        leaf = os.path.basename(c)
        rd = [n for n, txt in readers_src if leaf in txt]
        nt = [n for n, txt in notes_src if leaf in txt]
        diag = DIAG.split("/")[0] in c
        if rd and not diag: flagged += [(c, leaf, r) for r in rd]
        ev = camps[c][0][2].replace("|", "/")[:90]
        print(f"| {c}{' (sec. 4.7 diagnostic)' if diag else ''} | {len({x[0] for x in camps[c]})} | `{ev}` | {len(rd)}{': ' + ', '.join(rd[:3]) if rd else ''} | "
              f"{len(nt)}{': ' + ', '.join(os.path.basename(n) for n in nt[:3]) if nt else ''} |")
    launch = sorted({a[0] for a in acc if a[1] == "launcher"} - {"hspist3/validation/" + s for s in SELF})
    print("\n## Launchers that can select it\n")
    for l in launch:
        ev = [a[2] for a in acc if a[0] == l][0]
        print(f"- {l}: `{ev[:110]}`")
    nondiag = [c for c in camps if DIAG.split("/")[0] not in c]
    print(f"\ncampaigns with accelerated-backend runs (outside the sec. 4.7 diagnostic): {len(nondiag)}; "
          f"of them named by any analysis script, notebook or LaTeX source: {len({c for c, _, _ in flagged})}")
    print("\n## The readers flagged above, resolved\n")
    print("| campaign | reader | how it names the campaign | reader in the drafts' script closure (PART B) |\n|---|---|---|---|")
    bad = []
    for c, leaf, r in flagged:
        how = "only in its module docstring (an example layout)" if only_in_docstring(r, leaf) else "in code: READS IT"
        inpaper = r in closure
        if inpaper and "READS" in how: bad.append((c, r))
        print(f"| {c} | {r} | {how} | {'YES' if inpaper else 'no'} |")
    return acc, bad


def local_imports(p):
    try: tree = ast.parse(open(p, errors="ignore").read())
    except SyntaxError: return []
    mods = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import): mods += [a.name for a in n.names]
        elif isinstance(n, ast.ImportFrom) and n.module: mods.append(n.module)
    out = []
    for m in mods:
        for base in (os.path.dirname(p), HERE, HS):
            c = os.path.join(base, m.replace(".", "/") + ".py")
            if os.path.isfile(c): out.append(c); break
    return out


def constants(rel):
    """String constants of a Python file (f-string constant parts included, as ast visits them), or the words of a shell file."""
    txt = open(os.path.join(REPO, rel), errors="ignore").read()
    if rel.endswith(".sh"):
        return [(w, 0) for w in re.split(r"[\s\"'=]+", txt) if w]
    return [(n.value, n.lineno) for n in ast.walk(ast.parse(txt)) if isinstance(n, ast.Constant) and isinstance(n.value, str)]


def names_stem(rel, txt, stem):
    """The script names the figure stem literally, or (shell) through a file-name template such as name${g}_${mode}.bmp."""
    if stem in txt: return True
    if not rel.endswith(".sh"): return False
    for w in re.findall(r"[\w${}./-]*\$[\w${}./-]*", txt):
        base = os.path.splitext(os.path.basename(w))[0]
        parts = [t for t in re.split(r"(\$\{\w+\}|\$\w+)", base) if t]
        if "$" not in base or not parts or parts[0].startswith("$") or len(parts[0]) < 8: continue   # a literal prefix of >= 8 characters
        rx = "".join("[A-Za-z0-9_.-]+" if re.fullmatch(r"\$\{?\w+\}?", t) else re.escape(t) for t in parts)
        if re.fullmatch(rx, stem): return True
    return False


def closure_of_drafts(scripts):
    figs, seeds = [], set()
    for d in DRAFTS:
        s = open(os.path.join(REPO, d)).read()
        for st in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", s):
            stem = os.path.splitext(st)[0]
            pr = sorted(n for n, t in scripts if names_stem(n, t, stem))
            figs.append((os.path.basename(d), stem, pr)); seeds.update(pr)
        for sc in re.findall(r"%.*?([\w./-]+\.py)", s):
            for cand in (os.path.join("hspist3", sc), os.path.join("hspist3/validation", os.path.basename(sc)), sc):
                if os.path.isfile(os.path.join(REPO, cand)): seeds.add(os.path.normpath(cand)); break
    clo, todo = set(), [os.path.join(REPO, s) for s in seeds if s.endswith(".py")]
    while todo:
        p = todo.pop()
        if p in clo: continue
        clo.add(p); todo += local_imports(p)
    return figs, sorted(os.path.relpath(p, REPO) for p in clo) + sorted(s for s in seeds if s.endswith(".sh"))


def part_b(rows, figs, closure):
    print("\n## PART B -- the figure provenance of the two drafts\n")
    print("| draft | figure | scripts that name its file stem | note |\n|---|---|---|---|")
    for d, stem, pr in figs:
        note = NO_SCRIPT[stem][0] if stem in NO_SCRIPT else ("" if pr else "NO SCRIPT FOUND AND NOT EXPLAINED")
        for p in pr:
            if p.endswith(".sh"):      # a launcher: which backend its run commands select
                ks = sorted({k for k, _ in scan_text(os.path.join(REPO, p), "launcher") or [] if k != "code/doc mention"})
                note += f"{'; ' if note else ''}{os.path.basename(p)} runs the binary with: {', '.join(ks) or 'no --edmd-acc (default backend)'}"
        print(f"| {d} | {stem} | {', '.join(pr) or '-'} | {note} |")
    print(f"\nscript closure of the drafts (figure scripts, '% TODO-source' / '% FIGURES' scripts, their local imports): {len(closure)} files")
    for c in closure: print(f"- {c}")
    accf = {rel for rel, cls, h in rows if cls in ("run record", "run log", "summary csv") and DIAG not in rel
            and any(k == "ACCELERATED" for k, _ in h)}
    gen = set()
    for rel, cls, h in rows:
        if cls in ("run record", "run log", "summary csv") and rel not in accf: gen.update(rel.split(os.sep)[:-1])
    dist = set()
    for a in accf: dist.update(c for c in a.split(os.sep)[:-1] if c not in gen)
    consts = [(c, rel, ln) for rel in closure for c, ln in constants(rel)]
    t1 = sorted({(rel, ln, d) for c, rel, ln in consts for d in dist if d in c})
    print(f"\n### Test B1: constants of the closure that name a directory occurring only on accelerated-run paths\n")
    print(f"accelerated run files outside the sec. 4.7 diagnostic: {len(accf)}; distinctive directory names: {len(dist)}; "
          f"string constants in the closure: {len(consts)}; constants naming a distinctive directory: {len(t1)}")
    for rel, ln, d in t1: print(f"- {rel}:{ln} names {d}")
    # B2: data roots named by the closure
    idx = defaultdict(list)
    bases = [os.path.join(HS, d) for d in os.listdir(HS) if d.startswith("experiments_")]
    bases += [os.path.join(REPO, "0000_PLAN_OVERALL", p, "experiments") for p in ("paper1_speedofsound", "paper2_energytransfer")]
    for b in bases:
        for root, dirs, files in os.walk(b):
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
            for d in dirs: idx[d].append(os.path.relpath(os.path.join(root, d), REPO))
    named = sorted({p for c, _, _ in consts for p in re.split(r"[/*]", c) if p in idx and not GENERIC.match(p)})
    runrows = [(rel, cls, h) for rel, cls, h in rows if cls in ("run record", "run log")]
    wild = sorted({p for c, _, _ in consts for p in c.split("/") if re.search(r"[*?\[]", p) and p not in ("*", "**")})
    print(f"\n### Test B2: data roots named by the closure (directory names in its constants that exist in the data trees)\n")
    print("| directory name | instances | run records | --edmd-acc=0 explicit | ACCELERATED records | run logs | logs 'backend: default' "
          "| logs 'backend: accelerated' | accelerated children: matched by a wildcard of the closure |\n|---|---|---|---|---|---|---|---|---|")
    t2 = []
    for name in named:
        inst = idx[name]; rec = ex0 = accr = lg = lgd = lga = 0; kids = set()
        for rel, cls, h in runrows:
            if not any(rel.startswith(i + os.sep) for i in inst): continue
            ks = {k for k, _ in h}
            if cls == "run record":
                rec += 1; ex0 += "default (explicit)" in ks
                if "ACCELERATED" in ks:
                    accr += 1
                    i = next(i for i in inst if rel.startswith(i + os.sep)); kids.add(rel[len(i) + 1:].split(os.sep)[0])
            else:
                lg += 1; lgd += "default (log line)" in ks; lga += "ACCELERATED" in ks
        km = []
        for k in sorted(kids):
            m = [w for w in wild if fnmatch.fnmatchcase(k, w)]
            km.append(f"{k}: {', '.join(m) if m else 'none'}"); t2 += [(name, k, w) for w in m]
        print(f"| {name} | {len(inst)} | {rec} | {ex0} | {accr} | {lg} | {lgd} | {lga} | {'; '.join(km) or '-'} |")
    print(f"\nwildcard components in the closure's constants: {len(wild)}; accelerated children matched by one: {len(t2)}")
    for t in t2: print(f"- {t}")
    print("\n### Test B3: every directory-listing call of the closure (bases of the bare-wildcard calls to be checked by eye)\n")
    nbare = 0
    for rel in closure:
        if not rel.endswith(".py"): continue
        for i, l in enumerate(open(os.path.join(REPO, rel), errors="ignore").read().splitlines(), 1):
            s = l.strip()
            if s.startswith("#") or s.startswith("import") or not LISTING.search(s): continue
            bare = bool(re.search(r"[\"']\*{1,2}[\"']", s)); nbare += bare
            print(f"- {'BARE WILDCARD ' if bare else ''}{rel}:{i}: `{s[:150]}`")
    print(f"\nlisting calls with a bare '*' or '**' argument: {nbare}")
    print("\n### Test B4: data trees of the figures without a producing script\n")
    t4 = 0
    for stem, (what, trees) in NO_SCRIPT.items():
        for t in trees:
            sel = [(rel, cls, h) for rel, cls, h in runrows if fnmatch.fnmatchcase(rel.split(os.sep)[0] + os.sep + os.sep.join(rel.split(os.sep)[1:3]), t)]
            a = sum(1 for _, _, h in sel if any(k == "ACCELERATED" for k, _ in h)); t4 += a
            print(f"- {stem}: {t}: {sum(1 for x in sel if x[1] == 'run record')} run records, {sum(1 for x in sel if x[1] == 'run log')} run logs, "
                  f"ACCELERATED {a}")
        if not trees: print(f"- {stem}: {what}")
    return len(t1), len(t2), t4, nbare


def main():
    rows = scan_tree()
    scripts = sources((".py", ".sh"))
    figs, closure = closure_of_drafts(scripts)
    readers_src = [(n, t) for n, t in sources((".py", ".sh", ".tex", ".ipynb"))]
    acc, bad = part_a(rows, readers_src, set(closure))
    t1, t2, t4, nbare = part_b(rows, figs, closure)
    unexplained = [s for _, s, pr in figs if not pr and s not in NO_SCRIPT]
    ok = not bad and t1 == 0 and t2 == 0 and t4 == 0 and not unexplained
    print("\nVERDICT: " + (f"no figure of either draft, and no script that feeds one, reads a run of the accelerated backend (A: flagged readers "
                           f"outside the closure or docstring-only; B1 {t1}, B2 {t2}, B4 {t4}); {nbare} bare-wildcard listing calls are listed "
                           f"in B3 for a check of their bases by eye" if ok else
                           f"TRACE NEEDED -- readers in the closure {bad}, B1 {t1}, B2 {t2}, B4 {t4}, figures unexplained {unexplained}"))


if __name__ == "__main__":
    main()
