#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 51 item 3 b; 261012 sec. 4.7.114): A REHEARSAL OF THE MOVER (d51_move_untracked.py --go) on a
# synthetic layout in a NEW scratch folder: <T>/HardDisks (a main checkout) and <T>/HardDisks_gen3 (an engine-gen3 worktree of it) with
# every kind the real worktree has: tracked files, two run outputs rewritten (hspist3/energy_log.csv, hspist3/wall_position.csv), a
# populated gitlink hspist3/kissfft (a nested clone), untracked files in tracked and in new folders, ignored files and an ignored folder,
# untracked and ignored files whose path exists in the main checkout too, a symlink. Then: the list, --go, git worktree remove, and the
# checks (every file arrived once, nothing lost: the byte counts and SHA-256 of all files before and after, worktree gone, branch kept).
#   usage: bash d51_mover_rehearsal.sh <new scratch folder> <path of d51_move_untracked.py>
set -u
T=${1:?folder}; MV=${2:?mover}; [ -e "$T" ] && { echo "REFUSED: $T exists"; exit 1; }
mkdir -p "$T" && cd "$T" || exit 1
git init -q HardDisks && cd HardDisks && git config user.email t@t && git config user.name t
mkdir -p hspist3/cluster hspist3/exp_old && printf '*.pyc\n__pycache__/\nign_*\nbigdata/\n' > .gitignore
echo 'int main(void){return 0;}' > hspist3/00ALLINONE.c; echo t0 > hspist3/energy_log.csv; echo w0 > hspist3/wall_position.csv
echo c > hspist3/cluster/fetch.sh; echo e > hspist3/exp_old/kept.txt
(mkdir ../kiss && cd ../kiss && git init -q && git config user.email t@t && git config user.name t && echo k > k.c && git add k.c && git commit -qm k)
git add .gitignore hspist3 && git update-index --add --cacheinfo 160000,$(git -C ../kiss rev-parse HEAD),hspist3/kissfft && git commit -qm main0
git clone -q ../kiss hspist3/kissfft                                  # main's own populated kissfft, as on the Mac
echo "main's untracked" > hspist3/00_COMMAND.md; mkdir -p hspist3/__pycache__; echo pyc-main > hspist3/__pycache__/x.cpython-313.pyc
git branch engine-gen3 && git worktree add -q ../HardDisks_gen3 engine-gen3
G=../HardDisks_gen3
git -C $G config user.email t@t; git -C $G config user.name t
echo 'int main(void){return 1;}' > $G/hspist3/00ALLINONE.c && git -C $G commit -qam "gen3 work"     # committed (the mover refuses uncommitted sources)
mkdir -p $G/hspist3/experiments_gen3_x/tests/run1 $G/hspist3/experiments_gen3_x/bin_abc $G/hspist3/bigdata/sub $G/hspist3/__pycache__
echo t1 > $G/hspist3/energy_log.csv; echo w1 > $G/hspist3/wall_position.csv                       # the viewers' outputs
echo r > $G/hspist3/experiments_gen3_x/tests/run1/stdout.log; head -c 200000 /dev/urandom > $G/hspist3/experiments_gen3_x/bin_abc/00ALLINONE
echo d > $G/hspist3/exp_old/new_in_tracked_folder.txt; echo "gen3's untracked" > $G/hspist3/00_COMMAND.md
echo i > $G/hspist3/ign_note.txt; head -c 100000 /dev/urandom > $G/hspist3/bigdata/sub/raw.bin; echo pyc-gen3 > $G/hspist3/__pycache__/x.cpython-313.pyc
ln -s ../exp_old/kept.txt $G/hspist3/experiments_gen3_x/link_to_kept
git clone -q ../kiss $G/hspist3/kissfft; echo mod >> $G/hspist3/kissfft/k.c
( cd $G && git status --porcelain=v1 --ignored=matching --untracked-files=all --ignore-submodules=none | grep -v ' hspist3/kissfft$' | cut -c4- \
    | while read -r p; do if [ -d "$p" ]; then find "$p" -type f; else echo "$p"; fi; done | sort > ../mover_domain.txt
  xargs shasum -a 256 < ../mover_domain.txt | awk '{print $1}' | sort > ../before_sha.txt
  find hspist3/kissfft -path '*/.git' -prune -o -type f -print | xargs shasum -a 256 | awk '{print $1}' | sort >> ../before_sha.txt )
echo "the mover's domain: $(wc -l < ../mover_domain.txt | tr -d ' ') files plus the nested clone's working files" 
echo "== the list"; D51_REPO="$T" python3 "$MV" | sed 's/^/    /'
echo "== --go (its last two lines)"; D51_REPO="$T" python3 "$MV" --go | tail -2 | sed 's/^/    /'
echo "== git worktree remove ../HardDisks_gen3"; git worktree remove ../HardDisks_gen3; echo "    exit $?; worktree folder left: $([ -e ../HardDisks_gen3 ] && echo yes || echo no)"
echo "== git worktree list"; git worktree list | sed 's/^/    /'; echo "    branch engine-gen3: $(git rev-parse --short engine-gen3)"
echo "== checks"
for f in hspist3/experiments_gen3_x/tests/run1/stdout.log hspist3/experiments_gen3_x/bin_abc/00ALLINONE hspist3/exp_old/new_in_tracked_folder.txt \
         hspist3/ign_note.txt hspist3/bigdata/sub/raw.bin hspist3/experiments_gen3_x/link_to_kept; do echo "    main checkout: $f $([ -e "$f" ] || [ -L "$f" ] && echo present || echo MISSING)"; done
for f in hspist3/energy_log.csv hspist3/wall_position.csv hspist3/00_COMMAND.md hspist3/__pycache__/x.cpython-313.pyc hspist3/kissfft/k.c MANIFEST_d51_move.tsv; do
  echo "    hold: $f $([ -e "../HardDisks_d51_hold/$f" ] && echo "present ($(head -c 40 "../HardDisks_d51_hold/$f" | tr '\n' ' '))" || echo MISSING)"; done
echo "    main's own files untouched: 00_COMMAND.md '$(cat hspist3/00_COMMAND.md)', pyc '$(cat hspist3/__pycache__/x.cpython-313.pyc)', energy_log.csv '$(cat hspist3/energy_log.csv)'"
( cd .. && cat HardDisks/hspist3/experiments_gen3_x/tests/run1/stdout.log HardDisks/hspist3/experiments_gen3_x/bin_abc/00ALLINONE HardDisks/hspist3/exp_old/new_in_tracked_folder.txt \
    HardDisks/hspist3/ign_note.txt HardDisks/hspist3/bigdata/sub/raw.bin > /dev/null && echo "    the moved files read back" )
( cd .. && find HardDisks HardDisks_d51_hold -path '*/.git' -prune -o -type f -print | xargs shasum -a 256 | awk '{print $1}' | sort > after_sha.txt )
echo "    SHA-256 of the $(sort -u ../before_sha.txt | wc -l | tr -d ' ') distinct files the mover handled (untracked, ignored, run outputs, the nested clone's files): not found after: $(comm -23 <(sort -u ../before_sha.txt) ../after_sha.txt | wc -l | tr -d ' ')"
echo "    tracked files of the worktree are in git: engine-gen3 $(git rev-parse --short engine-gen3), $(git ls-tree -r engine-gen3 --name-only | wc -l | tr -d ' ') files"
