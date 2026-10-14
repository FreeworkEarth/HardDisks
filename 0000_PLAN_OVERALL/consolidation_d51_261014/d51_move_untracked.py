#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 51 item 3 b; 261012 sec. 4.7.114): THE STEP BEFORE `git worktree remove`. The engine-gen3
worktree (Repo/HardDisks_gen3) holds files outside git: untracked files (git worktree remove refuses them), ignored files (git worktree
remove deletes them with the worktree, silently: tested with git 2.50.1 on this Mac), the nested kissfft clone at the gitlink
hspist3/kissfft (git worktree remove refuses a populated gitlink) and tracked files modified by runs. This script moves every one of
them to the same relative path in the main checkout (Repo/HardDisks). It never deletes and never overwrites: a file whose target
exists there, the kissfft clone and the modified tracked files go to Repo/HardDisks_d51_hold/<same relative path> instead. Then the
tracked versions are restored in the worktree (git checkout -- <path>) and the empty gitlink folder is recreated, so that the worktree
is clean and `git worktree remove` can run without --force.
  default (--list): read-only. Prints the counts and sizes per folder and per destination, the files that go to the hold folder, the
                    modified tracked files, the processes whose working directory is inside the worktree, and the verdict.
  --go:             refuses unless the verdict of --list is READY; moves (os.rename: the same volume, nothing copied); writes the
                    manifest Repo/HardDisks_d51_hold/MANIFEST_d51_move.tsv (relative path, destination, bytes); checks that
                    `git status --porcelain --ignored --untracked-files=all` of the worktree is empty.
usage (from Repo/HardDisks): python3 0000_PLAN_OVERALL/consolidation_d51_261014/d51_move_untracked.py [--go]
"""
import collections, os, subprocess, sys

REPO = os.environ.get("D51_REPO", "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo")   # D51_REPO: a rehearsal layout only (d51_mover_rehearsal.sh)
G, M, HOLD = os.path.join(REPO, "HardDisks_gen3"), os.path.join(REPO, "HardDisks"), os.path.join(REPO, "HardDisks_d51_hold")
GITLINK = "hspist3/kissfft"
RUN_OUTPUTS = {"hspist3/energy_log.csv", "hspist3/wall_position.csv"}   # tracked files the viewers rewrite in their working directory


def git(*a, cwd=G):
    return subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True, check=True).stdout


def status():
    """(xy, path) of every entry, NUL-separated porcelain v1: untracked and ignored files one by one, ignored folders walked"""
    out = git("status", "--porcelain=v1", "-z", "--ignored=matching", "--untracked-files=all", "--ignore-submodules=none")
    ent = []
    for rec in out.split("\0"):
        if not rec: continue
        xy, p = rec[:2], rec[3:]
        if xy in ("??", "!!") and p.endswith("/"):
            for dp, dn, fn in os.walk(os.path.join(G, p)):
                for f in fn: ent.append((xy, os.path.relpath(os.path.join(dp, f), G)))
                for d in dn:
                    if os.path.islink(os.path.join(dp, d)): ent.append((xy, os.path.relpath(os.path.join(dp, d), G)))
        else:
            ent.append((xy, p))
    return ent


def cwd_procs():
    """processes whose working directory is inside the worktree (lsof -d cwd)"""
    r = subprocess.run(["lsof", "-d", "cwd", "-Fpcn"], capture_output=True, text=True).stdout.splitlines()
    out, pid, cmd = [], None, None
    for l in r:
        if l.startswith("p"): pid = l[1:]
        elif l.startswith("c"): cmd = l[1:]
        elif l.startswith("n") and (l[1:] == G or l[1:].startswith(G + "/")): out.append((pid, cmd, l[1:]))
    return out


def blocked(target):
    """the target exists, or a parent of it exists and is not a folder"""
    if os.path.lexists(target): return "exists"
    d = os.path.dirname(target)
    while d and d != M and d != HOLD and d != "/":
        if os.path.lexists(d) and not os.path.isdir(d): return f"parent {os.path.relpath(d, REPO)} is a file"
        d = os.path.dirname(d)
    return None


def plan():
    ent = status()
    moves, hold, problems = [], [], []
    for xy, p in ent:
        if p == GITLINK or p.startswith(GITLINK + "/"):
            continue
        if xy in ("??", "!!"):
            why = blocked(os.path.join(M, p))
            (hold if why else moves).append((xy, p, why or ""))
        elif xy.strip() == "M" and xy[0] == " " and p in RUN_OUTPUTS:
            hold.append((xy, p, "tracked file rewritten by runs; its committed version is restored in the worktree"))
        else:
            problems.append(f"{xy} {p}: a tracked change (commit it on engine-gen3 first)")
    gl = os.path.join(G, GITLINK)
    nested = os.path.isdir(os.path.join(gl, ".git"))
    for xy, p, _ in hold:
        if os.path.lexists(os.path.join(HOLD, p)): problems.append(f"hold target exists: HardDisks_d51_hold/{p}")
    if nested and os.path.lexists(os.path.join(HOLD, GITLINK)): problems.append(f"hold target exists: HardDisks_d51_hold/{GITLINK}")
    if git("diff", "--cached", "--name-only").strip(): problems.append("staged changes in the worktree (commit them first)")
    return moves, hold, nested, problems, cwd_procs()


def size(p):
    try: return os.lstat(os.path.join(G, p)).st_size
    except OSError: return 0


def report(moves, hold, nested, problems, procs):
    by = collections.defaultdict(lambda: [0, 0])
    for xy, p, _ in moves:
        parts = p.split("/"); k = "/".join(parts[:2]) if len(parts) > 2 else (parts[0] + "/" if len(parts) > 1 else "(top level)")
        by[(xy, k)][0] += 1; by[(xy, k)][1] += size(p)
    print(f"# d51_move_untracked.py {'--go: the list first' if '--go' in sys.argv else '--list (read-only)'}: {G} -> {M}\n")
    print("## to the main checkout, same relative path (?? untracked, !! ignored)\n")
    print("| kind | folder | files | MB |\n|---|---|---|---|")
    for (xy, k), (n, b) in sorted(by.items(), key=lambda kv: -kv[1][1]):
        print(f"| {xy} | {k} | {n} | {b / 1e6:.1f} |")
    tm = sum(size(p) for _, p, _ in moves)
    print(f"\n{len(moves)} files, {tm / 1e6:.1f} MB to the main checkout\n")
    print("## to Repo/HardDisks_d51_hold (the target exists in the main checkout, or a run output)\n")
    for xy, p, why in hold:
        print(f"- {xy} {p} ({size(p)} bytes): {why}")
    if nested:
        print(f"- the nested clone {GITLINK} (whole folder, {sum(os.lstat(os.path.join(dp, f)).st_size for dp, dn, fn in os.walk(os.path.join(G, GITLINK)) for f in fn) / 1e6:.1f} MB): "
              "a populated gitlink; the empty folder is recreated after the move")
    print(f"\n{len(hold) + int(nested)} entries to the hold folder\n")
    print("## processes whose working directory is inside the worktree (quit them before --go)\n")
    for pid, cmd, cwd in procs: print(f"- PID {pid} {cmd}: {cwd}")
    if not procs: print("- none")
    print("\n## problems\n")
    for x in problems: print(f"- {x}")
    if not problems: print("- none")
    ready = not problems and not procs
    print(f"\nVERDICT: {'READY' if ready else 'NOT READY'}")
    return ready


def go(moves, hold, nested):
    rows = []
    for xy, p, _ in moves:
        dst = os.path.join(M, p); os.makedirs(os.path.dirname(dst), exist_ok=True)
        b = size(p); os.rename(os.path.join(G, p), dst); rows.append((p, "main", b))
    restore = []
    for xy, p, why in hold:
        dst = os.path.join(HOLD, p); os.makedirs(os.path.dirname(dst), exist_ok=True)
        b = size(p); os.rename(os.path.join(G, p), dst); rows.append((p, "hold", b))
        if xy[0] == " ": restore.append(p)
    if nested:
        dst = os.path.join(HOLD, GITLINK); os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.rename(os.path.join(G, GITLINK), dst); os.makedirs(os.path.join(G, GITLINK)); rows.append((GITLINK + "/", "hold", 0))
    if restore: git("checkout", "--", *restore)
    os.makedirs(HOLD, exist_ok=True)
    with open(os.path.join(HOLD, "MANIFEST_d51_move.tsv"), "w") as f:
        f.write("path\tdestination\tbytes\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
    left, empty = [], 0   # an ignored folder emptied by the moves is still listed by git status --ignored; it holds nothing
    for rec in git("status", "--porcelain=v1", "-z", "--ignored=matching", "--untracked-files=all", "--ignore-submodules=none").split("\0"):
        if not rec: continue
        full = os.path.join(G, rec[3:])
        if rec[:2] == "!!" and rec.endswith("/") and os.path.isdir(full) and not any(
                fn or any(os.path.islink(os.path.join(dp, d)) for d in dn) for dp, dn, fn in os.walk(full)):
            empty += 1
        else:
            left.append(rec)
    print(f"moved {sum(r[1] == 'main' for r in rows)} files to the main checkout and {sum(r[1] == 'hold' for r in rows)} entries to "
          f"{HOLD}; restored {len(restore)} tracked files; manifest {HOLD}/MANIFEST_d51_move.tsv")
    print(f"GOOD: the worktree is clean ({empty} emptied ignored folders, which go with it); next: git worktree remove ../HardDisks_gen3"
          if not left else "NOT CLEAN, left in the worktree:\n" + "\n".join(left))


def main():
    moves, hold, nested, problems, procs = plan()
    ready = report(moves, hold, nested, problems, procs)
    if "--go" in sys.argv:
        if not ready: sys.exit("REFUSED: the verdict is NOT READY")
        go(moves, hold, nested)


if __name__ == "__main__":
    main()
