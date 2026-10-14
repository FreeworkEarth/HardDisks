#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 51 item 2; 261012 sec. 4.7.114): THE DRY RUN OF THE CONSOLIDATION main <- engine-gen3, in a
# scratch clone, never in the working trees. git clone --local (the objects are hardlinks of the repository's, no copy) without a
# checkout; a sparse checkout of two files only (the merge itself runs on the whole index; conflicted files are written regardless);
# git merge --no-ff --no-commit; then the checks of Report D51: the conflicts; for each conflicted file whether main's changes since the
# merge base are already in engine-gen3's version (every added line present, every removed line absent); the merged index against
# engine-gen3 and against main; the merge against the main checkout's local changes, untracked and ignored files (a merge refuses the
# first two and overwrites the third); the stage 1 manifests from the merged index; engine-divider-resched against engine-gen3; the
# runsheet's merge lines rehearsed on a branch of the clone; the clone's size. Deletes nothing; the clone stays.
#   usage: bash d51_dryrun.sh <new scratch folder>      (bash, not zsh: "${b}:${f}" revision syntax throughout)
set -uo pipefail
C=${1:?a new scratch folder}; [ -e "$C" ] && { echo "REFUSED: $C exists"; exit 1; }
M=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
W=$C.work; mkdir -p "$W"
echo "# d51_dryrun.sh, $(date '+%Y-%m-%d %H:%M:%S %Z'), $(git --version)"
echo; echo "## the heads (repository refs)"
for b in main engine-gen3 engine-divider-resched; do
  echo "- $b $(git -C $M rev-parse --short=8 $b), origin/$b $(git -C $M rev-parse --short=8 origin/$b), ahead/behind $(git -C $M rev-list --left-right --count $b...origin/$b | tr '\t' '/')"
done
MB=$(git -C $M merge-base main engine-gen3)
echo "- merge base main / engine-gen3: $(git -C $M log -1 --format='%h %ad %s' --date=short $MB | cut -c1-110)"
echo "- commits on main only: $(git -C $M rev-list --count engine-gen3..main); on engine-gen3 only: $(git -C $M rev-list --count main..engine-gen3)"
if git -C $M merge-base --is-ancestor engine-divider-resched engine-gen3; then echo "- engine-divider-resched IS an ancestor of engine-gen3"
else
  ME=$(git -C $M merge-base engine-divider-resched engine-gen3)
  echo "- engine-divider-resched is NOT an ancestor of engine-gen3: $(git -C $M rev-list --count engine-gen3..engine-divider-resched) commits not in it; base $(git -C $M log -1 --format='%h %ad' --date=short $ME)"
  n=0; d=""
  for f in $(git -C $M diff --name-only $ME engine-divider-resched); do
    n=$((n + 1)); [ "$(git -C $M rev-parse -q --verify "engine-divider-resched:${f}")" = "$(git -C $M rev-parse -q --verify "engine-gen3:${f}")" ] || d="$d $f"
  done
  echo "  files it changed since that base: $n; of them not identical in engine-gen3:${d:- none}"
fi
echo; echo "## the clone"
git clone -q --local --no-checkout "$M" "$C" && cd "$C" || exit 1
git sparse-checkout set --no-cone '/0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md' '/hspist3/00ALLINONE.c'
git checkout -q main
echo "- $C: main $(git rev-parse --short=8 HEAD), sparse (two files written)"
echo; echo "## git merge --no-ff --no-commit engine-gen3"
git merge --no-ff --no-commit origin/engine-gen3 > "$W/merge.txt" 2>&1; echo "- exit $?"; sed 's/^/    /' "$W/merge.txt"
git diff --name-only --diff-filter=U > "$W/conflicts.txt"
echo "- conflicts: $(wc -l < "$W/conflicts.txt" | tr -d ' ')"
while read -r f; do
  git -C $M diff "${MB}" main -- "$f" | grep '^+' | grep -v '^+++' | sed 's/^+//' > "$W/add.txt"
  git -C $M diff "${MB}" main -- "$f" | grep '^-' | grep -v '^---' | sed 's/^-//' > "$W/rem.txt"
  git show "origin/engine-gen3:${f}" > "$W/theirs.txt"
  na=$(wc -l < "$W/add.txt" | tr -d ' '); ma=0; pr=0
  while IFS= read -r l; do grep -qxF -- "$l" "$W/theirs.txt" || ma=$((ma + 1)); done < "$W/add.txt"
  while IFS= read -r l; do grep -qxF -- "$l" "$W/theirs.txt" && pr=$((pr + 1)); done < "$W/rem.txt"
  echo "  - $f: main's commits on it since the base: $(git -C $M log --format=%h ${MB}..main -- "$f" | tr '\n' ' ')"
  echo "    main added $na lines, $ma of them missing in engine-gen3's version; main removed $(wc -l < "$W/rem.txt" | tr -d ' ') lines, $pr of them still in it: $([ $ma -eq 0 ] && [ $pr -eq 0 ] && echo 'CONTAINED (engine-gen3 wins, nothing of main lost)' || echo 'NOT CONTAINED: STOP')"
  git checkout --theirs -- "$f" && git add "$f"
done < "$W/conflicts.txt"
echo; echo "## the merged index (conflicts taken from engine-gen3)"
echo "- against engine-gen3: $(git diff --cached --name-status origin/engine-gen3 | cut -f1 | sort | uniq -c | tr -s ' ' | tr '\n' ';')"
git diff --cached --name-only --diff-filter=A origin/engine-gen3 | awk -F/ '{ if ($1=="hspist3") print $1"/"$2; else print $1"/"$2"/"$3"/"$4 }' | sort | uniq -c | sed 's/^/    /'
echo "- files of engine-gen3 changed or removed by the merge: $(git diff --cached --name-only --diff-filter=MD origin/engine-gen3 | wc -l | tr -d ' ')"
echo "- against main: $(git diff --cached --name-status main | cut -f1 | sort | uniq -c | tr -s ' ' | tr '\n' ';') $(git diff --cached --shortstat main)"
git diff --cached --name-only --diff-filter=M main | sed 's/^/    M /'
echo; echo "## the merge against the main checkout"
git diff --cached --name-only main | sort > "$W/changed.txt"
git -C $M status --porcelain --untracked-files=all > "$W/mstat.txt"
grep -v '^??' "$W/mstat.txt" | cut -c4- | sed 's/^"//; s/"$//' | sort > "$W/dirty.txt"
grep '^??' "$W/mstat.txt" | cut -c4- | sed 's/^"//; s/"$//' | sort > "$W/untr.txt"
git -C $M status --porcelain --ignored=matching --untracked-files=all | grep '^!!' | cut -c4- | sed 's/^"//; s/"$//' | sort > "$W/ign.txt"
echo "- main checkout: $(wc -l < "$W/dirty.txt" | tr -d ' ') local changes, $(wc -l < "$W/untr.txt" | tr -d ' ') untracked, $(wc -l < "$W/ign.txt" | tr -d ' ') ignored entries"
echo "- in the merge's paths: local changes $(comm -12 "$W/dirty.txt" "$W/changed.txt" | wc -l | tr -d ' '), untracked $(comm -12 "$W/untr.txt" "$W/changed.txt" | wc -l | tr -d ' '), ignored $(comm -12 "$W/ign.txt" "$W/changed.txt" | wc -l | tr -d ' ')"
comm -12 "$W/dirty.txt" "$W/changed.txt" | sed 's/^/    local change: /'
echo; echo "## the stage 1 manifests from the merged index"
for m in s1_inputs.sha256 s1_merge_inputs.sha256; do
  git show ":hspist3/cluster/gen3_stage1_261012/$m" | while read -r h f; do
    a=$(git show ":hspist3/${f}" | shasum -a 256 | cut -d' ' -f1); [ "$a" = "$h" ] && echo OK || echo "FAIL $f"; done | sort | uniq -c | sed "s/^/- $m: /"
done
echo "- hspist3/00ALLINONE.c: merged $(git rev-parse :hspist3/00ALLINONE.c | cut -c1-12), engine-gen3 $(git rev-parse origin/engine-gen3:hspist3/00ALLINONE.c | cut -c1-12)"
echo; echo "## the runsheet's merge lines, rehearsed on a branch of the clone"
git merge --abort
git checkout -q -b d51_runsheet_rehearsal main
git merge --no-ff -m "consolidate: engine-gen3 into main after stage 1 (decision 51)" origin/engine-gen3 > "$W/merge2.txt" 2>&1; echo "- git merge: exit $?, $(grep -c '^CONFLICT' "$W/merge2.txt") CONFLICT line(s): $(grep '^CONFLICT' "$W/merge2.txt" | tr '\n' ' ')"
while read -r f; do git checkout --theirs -- "$f" && git add "$f"; done < "$W/conflicts.txt"
git commit -q -m "consolidate: engine-gen3 into main after stage 1 (decision 51)"
echo "- the commit: $(git log -1 --format='%h, parents %p, "%s"')"
echo "- its tree against engine-gen3: $(git diff --name-status origin/engine-gen3 HEAD | cut -f1 | sort | uniq -c | tr -s ' ' | tr '\n' ';')"
echo; echo "## the size"
echo "- the clone without .git: $(du -sh -I .git "$C" | cut -f1); .git: $(du -sh "$C/.git" | cut -f1), of it new (not hardlinked) $(find "$C/.git" -type f -links 1 -print0 | xargs -0 stat -f %z | awk '{s+=$1} END {printf "%.1f MB", s/1e6}')"
echo "- the work files: $(du -sh "$W" | cut -f1) in $W"
