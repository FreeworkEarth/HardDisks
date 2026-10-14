#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 51 item 3 b; 261012 sec. 4.7.114): what `git worktree remove` does on this Mac's git with
# (1) a clean worktree, (2) a populated gitlink (a nested clone, as hspist3/kissfft), (2b) the clone moved out, (2c) its empty folder
# recreated, (3) only an ignored file. A throwaway repository in a NEW scratch folder; nothing outside it is touched.
#   usage: bash d51_worktree_remove_test.sh <new scratch folder>
set -u
T=${1:?a new scratch folder}; [ -e "$T" ] && { echo "REFUSED: $T exists"; exit 1; }
mkdir -p "$T" && cd "$T" && git init -q repo && cd repo && git config user.email t@t && git config user.name t
echo a > a.txt && printf 'ignored_*\n' > .gitignore && git add a.txt .gitignore && git commit -qm init
(mkdir subsrc && cd subsrc && git init -q && git config user.email t@t && git config user.name t && echo s > s.txt && git add s.txt && git commit -qm s)
git update-index --add --cacheinfo 160000,$(git -C subsrc rev-parse HEAD),sub && git commit -qm gitlink
git branch b1; git branch b2; git branch b3
echo "# $(git --version)"
git worktree add -q ../wt1 b1; echo "(1) clean worktree, gitlink not populated:"; git worktree remove ../wt1; echo "    exit $?, folder left: $([ -e ../wt1 ] && echo yes || echo no)"
git worktree add -q ../wt2 b2; git clone -q subsrc ../wt2/sub
echo "(2) populated gitlink (a nested clone):"; git worktree remove ../wt2 2>&1 | sed 's/^/    /'; echo "    folder left: $([ -e ../wt2 ] && echo yes || echo no)"
mv ../wt2/sub ../sub_moved_out; echo "(2b) the nested clone moved out:"; git -C ../wt2 status --porcelain --ignore-submodules=none | sed 's/^/    status: /'; git worktree remove ../wt2 2>&1 | sed 's/^/    /'
mkdir ../wt2/sub; echo "(2c) its empty folder recreated:"; git worktree remove ../wt2 2>&1 | sed 's/^/    /'; echo "    folder left: $([ -e ../wt2 ] && echo yes || echo no)"
git worktree add -q ../wt3 b3; echo x > ../wt3/ignored_data.txt; echo "(3) only an ignored file:"; git -C ../wt3 status --porcelain --ignored | sed 's/^/    status: /'
git worktree remove ../wt3 2>&1 | sed 's/^/    /'; echo "    exit ${PIPESTATUS[0]}, folder left: $([ -e ../wt3 ] && echo yes || echo no) (the ignored file went with it)"
