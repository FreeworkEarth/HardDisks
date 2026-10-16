#!/bin/sh
# tests/run_tests.sh -- Demo 0a tests of brief 01_DEMO0_BRIEF.md sec. 6 (items 1-4), plus the rules of sec. 2 that a test can
# see. Run through `make test` from demo0/. Headless throughout; the GUI = headless test uses SDL's dummy video driver, which
# opens no window on screen. Outputs go to a new temporary folder that is left in place (the OS cleans its temp dir);
# nothing is written into the repository and nothing is deleted.
set -u
BIN=./bin/pid_demo0
CHK=./bin/pid_check
T=$(mktemp -d "${TMPDIR:-/tmp}/pid_demo0_test.XXXXXX") || exit 1
PASS=0; FAIL=0
ok()  { PASS=$((PASS + 1)); echo "PASS $*"; }
bad() { FAIL=$((FAIL + 1)); echo "FAIL $*"; }
echo "test outputs: $T"

# ---- rules: one generator, no libc rand(); -ffp-contract=off in the build
if nm -u $BIN $CHK | grep -qE '^_?s?rand$'; then bad "rule 1: a binary imports libc rand()/srand()"; else ok "rule 1: neither binary imports libc rand() or srand() (nm -u)"; fi
if $BIN --version | grep -q -- '-ffp-contract=off'; then ok "rule 2: build line carries -ffp-contract=off"; else bad "rule 2: -ffp-contract=off missing from the build line"; fi

# ---- 1. replay: same seed and flags -> every state file identical; another seed -> final.txt differs
A="--seed=9300 --steps=150 --write-every=10 --headless"
$BIN $A --out=$T/rep_a >/dev/null && $BIN $A --out=$T/rep_b >/dev/null && $BIN --seed=9301 --steps=150 --write-every=10 --headless --out=$T/rep_c >/dev/null
n=0; d=0
for f in $T/rep_a/state_*.txt $T/rep_a/final.txt; do n=$((n + 1)); cmp -s "$f" "$T/rep_b/$(basename "$f")" || d=$((d + 1)); done
if [ "$n" -ge 16 ] && [ "$d" -eq 0 ]; then ok "1 replay: seed 9300 twice, $n state files cmp-identical"; else bad "1 replay: $d of $n files differ"; fi
if cmp -s $T/rep_a/final.txt $T/rep_c/final.txt; then bad "1 replay: seed 9301 gives the same final.txt as 9300"; else ok "1 replay: seed 9301 gives a different final.txt"; fi

# ---- 1. GUI = headless: a window run with --steps=N writes the same files as a headless run
SDL_VIDEODRIVER=dummy $BIN --seed=9300 --steps=150 --write-every=10 --tpf=7 --tick-ms=0 --exit-at-end --out=$T/gui >$T/gui.log 2>&1
rc=$?; n=0; d=0
for f in $T/rep_a/state_*.txt $T/rep_a/final.txt; do n=$((n + 1)); cmp -s "$f" "$T/gui/$(basename "$f")" || d=$((d + 1)); done
if [ "$rc" -eq 0 ] && [ "$d" -eq 0 ] && grep -q '^mode window' $T/gui/manifest.txt; then ok "1 GUI = headless: window run (dummy driver, 7 ticks per frame) wrote the same $n state files"; else bad "1 GUI = headless: rc $rc, $d of $n files differ (log $T/gui.log)"; fi

# ---- 0a.1 regression: the 0a parameters reproduce the committed 0a result (tests/ref; sha256 84a6a12f..., the final.txt of
# the 0a screenshot run, demo0/shots/demo0a_seed9300_tick200_manifest.txt) in every line except params and build
REF=tests/ref/demo0a_seed9300_tick200_final.txt
$BIN --seed=9300 --steps=200 --headless --start-rings=0 --pbranch-exp=0 --pbranch=0.08 --ntry=8 --budget=3000 --hull=off --out=$T/reg0a >/dev/null
if [ "$(shasum -a 256 $REF | cut -c1-64)" = 84a6a12fa3f500fb88ec44e2f3b532b1e482fdfd730fa4246300fd2ab4f9d851 ] && \
   grep -v -e '^params ' -e '^build ' $REF > $T/reg_ref.txt && grep -v -e '^params ' -e '^build ' $T/reg0a/final.txt > $T/reg_new.txt && \
   cmp -s $T/reg_ref.txt $T/reg_new.txt; then
  ok "0a.1 regression: --start-rings=0 --pbranch-exp=0 with the 0a values reproduces the 0a final.txt (84a6a12f) except params and build"
else bad "0a.1 regression: the 0a result is not reproduced (diff $T/reg_ref.txt $T/reg_new.txt)"; fi

# ---- 0a.2 presets: --preset=mat is exactly its flags; an explicit flag overrides the preset in either order
$BIN --seed=9300 --steps=60 --headless --preset=mat --out=$T/pre_mat >/dev/null
$BIN --seed=9300 --steps=60 --headless --pbranch=1 --pbranch-exp=0 --delta=1.57 --ntry=48 --budget=9000 --out=$T/pre_flags >/dev/null
$BIN --seed=9300 --steps=60 --headless --budget=4500 --preset=mat --out=$T/pre_over >/dev/null
if cmp -s $T/pre_mat/final.txt $T/pre_flags/final.txt && grep -q '^budget 4500$' $T/pre_over/final.txt && grep -q '^preset mat$' $T/pre_mat/manifest.txt && \
   ! $BIN --seed=1 --steps=1 --headless --preset=sheet --out=$T/pre_bad >/dev/null 2>&1 && [ ! -e $T/pre_bad ]; then
  ok "0a.2 presets: --preset=mat equals its flags, an explicit --budget overrides it, an unknown preset is refused"
else bad "0a.2 presets"; fi
$BIN --seed=9300 --steps=60 --write-every=1 --headless --preset=mat --out=$T/inv_mat >/dev/null
out=$($CHK verify $T/inv_mat/state_*.txt $T/inv_mat/final.txt); rc=$?
if [ "$rc" -eq 0 ]; then ok "2+3 invariants and ledger, --preset=mat, seed 9300, every tick: $out"; else bad "2+3 --preset=mat: $out"; fi

# ---- 0b (brief sec. 3 as amended by sec. 8): --preset=mat45, patches A under the seed and B, adaptation, waiting tips,
# random tip order. 5: guided with beta=0 equals unbiased in every line except rule and params; a fed path by tick 400.
G="--preset=mat45 --steps=400 --write-every=50 --headless"
$BIN --seed=9300 $G --rule=unbiased --out=$T/b_unb >/dev/null
$BIN --seed=9300 $G --rule=guided --beta=0 --out=$T/b_g0 >/dev/null
n=0; d=0
for f in $T/b_unb/state_*.txt $T/b_unb/final.txt; do
  n=$((n + 1)); grep -v -e '^rule ' -e '^params ' "$f" > $T/b1.txt; grep -v -e '^rule ' -e '^params ' "$T/b_g0/$(basename "$f")" > $T/b2.txt; cmp -s $T/b1.txt $T/b2.txt || d=$((d + 1))
done
if [ "$d" -eq 0 ] && ! cmp -s $T/b_unb/final.txt $T/b_g0/final.txt; then ok "5 guided --beta=0 = unbiased (mat45, seed 9300): $n state files identical except rule and params"; else bad "5 guided beta=0: $d of $n files differ"; fi
$BIN --seed=9300 $G --rule=guided --out=$T/b_g >/dev/null
info=$($CHK info $T/b_g/state_000400.txt)
if echo "$info" | awk '{for (i=1;i<=NF;i++) if ($i ~ /^touch/) t += $(i+1)} END {exit !(t > 0)}'; then ok "5 guided (default beta): patch B reached by tick 400 on seed 9300 ($info)"; else bad "5 guided: no fed path by tick 400 ($info)"; fi
SDL_VIDEODRIVER=dummy $BIN --seed=9300 --preset=mat45 --rule=guided --steps=120 --write-every=10 --tpf=9 --tick-ms=0 --exit-at-end --out=$T/b_gui >$T/b_gui.log 2>&1
$BIN --seed=9300 --preset=mat45 --rule=guided --steps=120 --write-every=10 --headless --out=$T/b_hl >/dev/null
n=0; d=0; for f in $T/b_hl/state_*.txt $T/b_hl/final.txt; do n=$((n + 1)); cmp -s "$f" "$T/b_gui/$(basename "$f")" || d=$((d + 1)); done
if [ "$d" -eq 0 ]; then ok "1 GUI = headless under 0b (mat45, guided, random tip order): $n state files identical"; else bad "1 GUI = headless 0b: $d of $n differ"; fi
for v in "b400:--preset=mat45 --rule=guided --steps=400" "bfast:--preset=mat45 --rule=guided --steps=300 --deltathin=0.2 --gamma=0.1"; do
  n=${v%%:*}; a=${v#*:}
  $BIN --seed=9301 --write-every=1 --headless --out=$T/inv_$n $a >/dev/null
  out=$($CHK verify $T/inv_$n/state_*.txt $T/inv_$n/final.txt); rc=$?
  last=$($CHK info $T/inv_$n/final.txt)
  if [ "$rc" -eq 0 ]; then ok "2+3 0b invariants and ledger every tick, seed 9301, $a: $out; end: $last"; else bad "2+3 0b $a: $out"; fi
done
d=0; for f in $T/b_g/state_000050.txt $T/b_g/final.txt $T/inv_bfast/final.txt; do $CHK roundtrip "$f" $T/rt_b.txt && cmp -s "$f" $T/rt_b.txt || d=$((d + 1)); done
if [ "$d" -eq 0 ]; then ok "4 round trip with [patches] and removed disks: 3 files byte-identical"; else bad "4 round trip 0b: $d files differ"; fi

# ---- 0c (brief sec. 3 "0c", sec. 6 test 6, sec. 9): the hull in final.txt, checked by pid_check hull: (a) every loop
# closed to 1e-9, (b) every disk centre inside, (c) raster area at 2/4/8/16 px per sigma (finest within 0.1 %, falling
# log-log trend), (e) reversed edge order gives the same primitives and area, and the stored hull equals a recomputation
H6=""
for v in "0a:--preset=0a --seed=9300 --steps=200" "mat:--preset=mat --seed=9300 --steps=60" "mat45_400:--preset=mat45 --rule=guided --seed=9300 --steps=400" \
         "mat45_1000:--preset=mat45 --rule=guided --seed=9300 --steps=1000" "thin05:--preset=mat45 --rule=guided --seed=9301 --steps=400 --deltathin=0.05" \
         "rings3:--start-rings=3 --budget=6000 --seed=9302 --steps=200" "k1:--preset=0a --seed=7 --steps=1 --K=1"; do
  n=${v%%:*}; a=${v#*:}
  $BIN --headless --write-every=0 --out=$T/h_$n $a >/dev/null
  out=$($CHK hull $T/h_$n/final.txt); rc=$?
  short=$(echo "$out" | sed 's#^[^:]*: ##')
  if [ "$rc" -eq 0 ]; then ok "6 hull ($a): $short"; else bad "6 hull ($a): $out"; fi
  H6="$H6 $T/h_$n/final.txt"
done
$BIN --headless --preset=0a --seed=9300 --steps=40 --write-every=10 --hull=all --out=$T/h_all >/dev/null
out=$($CHK hull $T/h_all/state_*.txt $T/h_all/final.txt); rc=$?
if [ "$rc" -eq 0 ] && [ "$(grep -l '^\[hull\]' $T/h_all/state_*.txt | wc -l | tr -d ' ')" -eq 5 ]; then ok "6 hull --hull=all: a hull in every state file, all 6 pass"; else bad "6 hull --hull=all: $out"; fi
d=0; for f in $H6; do $CHK roundtrip "$f" $T/rt_h.txt && cmp -s "$f" $T/rt_h.txt || d=$((d + 1)); done
if [ "$d" -eq 0 ]; then ok "4 round trip with [hull] and [areas]: 7 files byte-identical"; else bad "4 round trip with hull: $d files differ"; fi
out=$(python3 tests/hull_reference.py $T/h_mat45_400/final.txt 2>&1); rc=$?
case "$out" in SKIP*) echo "$out";; *) if [ "$rc" -eq 0 ]; then ok "$out"; else bad "$out"; fi;; esac
SDL_VIDEODRIVER=dummy $BIN --preset=mat45 --rule=guided --seed=9300 --steps=100 --write-every=0 --tpf=9 --tick-ms=0 --view=network --exit-at-end --out=$T/h_gui >$T/h_gui.log 2>&1
$BIN --preset=mat45 --rule=guided --seed=9300 --steps=100 --write-every=0 --headless --out=$T/h_hl >/dev/null
if cmp -s $T/h_gui/final.txt $T/h_hl/final.txt; then ok "1 GUI = headless with --view=network (hull drawn every frame): final.txt with its hull identical"; else bad "1 GUI = headless, view network"; fi

# ---- 2 and 3. invariants and ledger on every written state (every tick), and the tick-to-tick sequence
$BIN --seed=9300 --steps=200 --write-every=1 --headless --out=$T/inv >/dev/null
out=$($CHK verify $T/inv/state_*.txt $T/inv/final.txt); rc=$?
if [ "$rc" -eq 0 ]; then ok "2+3 invariants and ledger, seed 9300, every tick: $out"; else bad "2+3 seed 9300: $out"; fi
for s in 9301 9302 9303 9304 9305 1 2 3; do
  $BIN --seed=$s --steps=200 --write-every=5 --headless --out=$T/inv_$s >/dev/null
  out=$($CHK verify $T/inv_$s/state_*.txt $T/inv_$s/final.txt); rc=$?
  if [ "$rc" -eq 0 ]; then ok "2+3 invariants and ledger, seed $s, every 5th tick: $out"; else bad "2+3 seed $s: $out"; fi
done
for v in "0a:--start-rings=0 --pbranch-exp=0 --pbranch=0.08 --ntry=8 --budget=3000" "rings0:--start-rings=0" "rings2:--start-rings=2" "rings3:--start-rings=3 --budget=6000"; do
  n=${v%%:*}; a=${v#*:}
  $BIN --seed=9302 --steps=200 --write-every=1 --headless --out=$T/var_$n $a >/dev/null
  out=$($CHK verify $T/var_$n/state_*.txt $T/var_$n/final.txt); rc=$?
  if [ "$rc" -eq 0 ]; then ok "2+3 invariants and ledger, seed 9302, every tick, $a: $out"; else bad "2+3 $a: $out"; fi
done
# the checker must catch a broken state: move one disk onto its non-adjacent neighbour, break the ledger
awk 'BEGIN{done=0} /^\[disks\]/{s=1} /^\[tips\]/{s=0} { if (s && !done && $1 == "2") { $2 = "100"; $3 = "60"; done = 1 } print }' $T/inv/final.txt > $T/broken_overlap.txt
sed 's/^m_used .*/m_used 1234.5/' $T/inv/final.txt > $T/broken_ledger.txt
if ! $CHK verify $T/broken_overlap.txt >/dev/null && ! $CHK verify $T/broken_ledger.txt >/dev/null; then ok "2+3 the checker rejects a moved disk and a broken ledger"; else bad "2+3 the checker accepted a broken state"; fi

# ---- 4. state round trip: read -> write is byte-identical
d=0
for f in $T/inv/state_000000.txt $T/inv/state_000050.txt $T/inv/final.txt $T/inv_9303/final.txt; do
  $CHK roundtrip "$f" $T/roundtrip.txt && cmp -s "$f" $T/roundtrip.txt || { d=$((d + 1)); echo "  differs: $f"; }
done
if [ "$d" -eq 0 ]; then ok "4 round trip: 4 state files read and written back byte-identical"; else bad "4 round trip: $d files differ"; fi

# ---- manifest: the sha256 of every state file
(cd $T/inv && grep -v '^#' manifest.txt | grep -E '^[0-9a-f]{64}  ' | shasum -a 256 -c --quiet) && ok "manifest: every listed sha256 matches its state file ($(grep -cE '^[0-9a-f]{64}  ' $T/inv/manifest.txt) files)" || bad "manifest: sha256 mismatch"

# ---- rule 4 and the command line: nothing in the cwd, no reuse of a full run folder, exact flags
mkdir $T/cwd && (cd $T/cwd && "$OLDPWD/$BIN" --seed=5 --steps=20 --headless --out=$T/cwd_run >/dev/null)
if [ -z "$(ls -A $T/cwd)" ]; then ok "rule 4: a run writes nothing into its cwd"; else bad "rule 4: files appeared in the cwd: $(ls -A $T/cwd)"; fi
before=$(shasum -a 256 $T/rep_a/final.txt)
$BIN --seed=1 --steps=5 --headless --out=$T/rep_a >/dev/null 2>&1; rc=$?
if [ "$rc" -eq 2 ] && [ "$before" = "$(shasum -a 256 $T/rep_a/final.txt)" ]; then ok "rule 4: a non-empty --out is refused and left unchanged"; else bad "rule 4: non-empty --out not refused (rc $rc)"; fi
r1=0; for a in "--see=1" "--seed" "--steps=5x" "--sed=1" "--seed=-1" "--R0=0" "--rule=fancy" "--start-rings=-1" "--start-rings=30" "--pbranch-exp=-1" "--patch=1,2,3" "--tip-order=up" "--hull=some" "--alpha=2" "--view=tree"; do
  $BIN --seed=1 --steps=5 --headless --out=$T/cli_bad $a >/dev/null 2>&1 && r1=1; done
if [ "$r1" -eq 0 ] && [ ! -e $T/cli_bad ]; then ok "command line: 15 bad or unknown flags refused before any folder is made"; else bad "command line: a bad flag was accepted"; fi

echo "tests: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
