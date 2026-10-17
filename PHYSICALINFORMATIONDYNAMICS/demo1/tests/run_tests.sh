#!/bin/sh
# tests/run_tests.sh -- Demo 1 delivery 1a: tests 1-2 of brief 02_DEMO1_BRIEF.md sec. 6, with the corrections of
# 02_DEMO1_BRIEF_PRECHECK.md part A. Run through `make test` from demo1/; the Makefile builds bin/ first.
#   Test 2  each engine file is its original plus the exact 4-line header (whose commits and blob git confirms), byte for byte;
#           engine/ holds exactly the three files and README_FORK.md.
#   Test 1  the fork and the unmodified original gen3, each linked into its own build of the unchanged M1/M2 harness sources
#           (gen2 linked into both as one object) and of tests/box_final_state.c, give byte-identical stdout for M1 audit,
#           M1 audit --quick, M2 audit --quick and the final-state driver, each a separate process, with empty stderr; each
#           output also passes its own verdict checks (cell counts, verdict lines, health counters, exact line count, end),
#           because two equally broken outputs would compare equal (pre-check A5). The fork's object links into a program
#           with nothing but libc and libm.
# The eight processes run in parallel, each capped at 900 s (perl alarm: macOS has no `timeout`); Ctrl-C stops them all.
# Outputs: a new folder runs/box_reference_<cc>_<UTC time>/ (git-ignored) with manifest.txt. Nothing else is written and
# nothing is deleted. One `make test` at a time per checkout: bin/ is shared.
set -u
: "${CC:?run through make test}" "${IS_CLANG:?run through make test}" "${CSTD:?run through make test}"
HSP=../../hspist3/edmd_core
D=$(pwd)
if [ "$IS_CLANG" = 1 ]; then CCID=clang; else CCID=gcc; fi
OUT=runs/box_reference_${CCID}_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p runs && mkdir "$OUT" || { echo "cannot create a new $OUT (another run started this second?)"; exit 1; }
PASS=0; FAIL=0; LOG=$OUT/results.txt
ok()   { PASS=$((PASS + 1)); echo "PASS $*" | tee -a "$LOG"; }
bad()  { FAIL=$((FAIL + 1)); echo "FAIL $*" | tee -a "$LOG"; }
note() { echo "INFO $*" | tee -a "$LOG"; }
now()  { perl -MTime::HiRes=time -e 'printf "%.2f\n", time'; }
sha()  { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -c1-64; else sha256sum "$1" | cut -c1-64; fi; }
nlines() { wc -l < "$1" | tr -d ' '; }
echo "test outputs: $OUT"

# ------------------------------------------------------------------ what is about to be tested (written before the runs)
{
  echo "date_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host $(uname -srm), $(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu) cpus"
  echo "repo_head $(git rev-parse HEAD 2>/dev/null); demo1/ entries not committed: $(git status --porcelain --untracked-files=all -- . 2>/dev/null | wc -l | tr -d ' ')"
  echo "compiler CC=$CC is_clang=$IS_CLANG: $($CC --version 2>/dev/null | head -1)"
  echo "flags $CSTD"
  echo "## sources (sha256)"
  for f in Makefile tests/run_tests.sh tests/box_final_state.c engine/README_FORK.md engine/edmd_gen3.c engine/edmd_gen3.h engine/edmd.h; do echo "$(sha "$f")  demo1/$f"; done
  echo "## originals (git blob)"
  for f in edmd_gen3.c edmd_gen3.h edmd.h edmd.c tests/gen3_m1_harness.c tests/gen3_m2_harness.c tests/gen3_state_io.h; do echo "$(git hash-object "$HSP/$f")  hspist3/edmd_core/$f"; done
  echo "## the tested objects and programs (sha256)"
  for f in bin/gen2.o bin/orig/edmd_gen3.o bin/fork/edmd_gen3.o bin/orig/gen3_m1_harness bin/fork/gen3_m1_harness bin/orig/gen3_m2_harness \
           bin/fork/gen3_m2_harness bin/orig/box_final_state bin/fork/box_final_state; do echo "$(sha "$f")  demo1/$f"; done
  echo "## build: every command make runs for 'all' (make -n -B), in order"
  cat bin/build_lines.txt 2>/dev/null || echo "(bin/build_lines.txt missing: run through make test)"
} > "$OUT/snapshot.txt"

# ------------------------------------------------------------------ test 2: the fork's sources
t2=0
for f in edmd_gen3.c edmd_gen3.h edmd.h; do
  fk=engine/$f; or=$HSP/$f
  cur=$(git hash-object "$or")
  l2=$(sed -n 2p "$fk")
  cpy=$(echo "$l2" | sed -n 's/^   copied at origin\/main \([0-9a-f]\{40\}\); last change to it \([0-9a-f]\{40\}\);$/\1/p')
  lst=$(echo "$l2" | sed -n 's/^   copied at origin\/main \([0-9a-f]\{40\}\); last change to it \([0-9a-f]\{40\}\);$/\2/p')
  {
    echo "/* PID fork, PHYSICALINFORMATIONDYNAMICS/demo1/engine (Demo 1, delivery 1a, 2026-10-16): hspist3/edmd_core/$f"
    echo "   copied at origin/main $cpy; last change to it $lst;"
    echo "   source blob $cur. Lines 5 on are that file unchanged: its line n is line n + 4 here."
    echo '   Changes from delivery 1b on are marked at their site by a comment starting "PID:" (README_FORK.md). */'
  } > "$OUT/header_$f.txt"
  if [ -z "$cpy" ]; then
    bad "2 source: $fk line 2 is not '   copied at origin/main <commit>; last change to it <commit>;'"; t2=1
  elif ! head -n 4 "$fk" | cmp -s - "$OUT/header_$f.txt"; then
    bad "2 source: $fk lines 1-4 are not the header template with the original's current blob $cur (diff $OUT/header_$f.txt)"; t2=1
  elif ! tail -n +5 "$fk" | cmp -s - "$or"; then
    bad "2 source: $fk from line 5 differs from $or"; t2=1
  elif [ "$(git rev-parse "$cpy:hspist3/edmd_core/$f" 2>/dev/null)" != "$cur" ]; then
    bad "2 source: at the recorded copy commit $cpy, hspist3/edmd_core/$f is not blob $cur"; t2=1
  elif [ "$(git log -1 --format=%H "$cpy" -- "$or" 2>/dev/null)" != "$lst" ]; then
    bad "2 source: the last change to hspist3/edmd_core/$f before $cpy is not the recorded $lst"; t2=1
  fi
done
extra=$(ls -A engine | grep -v -x -e edmd_gen3.c -e edmd_gen3.h -e edmd.h -e README_FORK.md)
if [ -n "$extra" ]; then bad "2 source: engine/ holds other entries: $(echo $extra)"; t2=1; fi
if [ ! -f engine/README_FORK.md ]; then bad "2 source: engine/README_FORK.md is missing"; t2=1; fi
if [ "$t2" -eq 0 ]; then ok "2 source: edmd_gen3.c, edmd_gen3.h, edmd.h = the exact 4-line header (commits and blobs confirmed by git) plus their originals, byte for byte"; fi

# ------------------------------------------------------------------ test 1a: the fork alone
nmok=1
nm -u bin/fork/edmd_gen3.o > "$OUT/nm_raw_fork.txt" 2>&1 || nmok=0
nm -u bin/orig/edmd_gen3.o > "$OUT/nm_raw_orig.txt" 2>&1 || nmok=0
awk '{print $NF}' "$OUT/nm_raw_fork.txt" | sort > "$OUT/nm_u_fork.txt"
awk '{print $NF}' "$OUT/nm_raw_orig.txt" | sort > "$OUT/nm_u_orig.txt"
printf 'int main(void){ return 0; }\n' > "$OUT/stub_main.c"
$CC $CSTD -o "$OUT/stub_link" "$OUT/stub_main.c" bin/fork/edmd_gen3.o -lm > "$OUT/stub_link.txt" 2>&1; lrc=$?
if [ "$lrc" -eq 0 ] && "$OUT/stub_link" && [ "$nmok" -eq 1 ] && [ -s "$OUT/nm_u_fork.txt" ] && cmp -s "$OUT/nm_u_fork.txt" "$OUT/nm_u_orig.txt"; then
  ok "1 fork alone: the fork's object links into a stub program with libc and -lm only; its $(nlines "$OUT/nm_u_fork.txt") undefined symbols equal the original object's"
else bad "1 fork alone: stub link rc $lrc ($OUT/stub_link.txt), nm ok $nmok, or nm -u empty or different from the original object's"; fi
if cmp -s bin/fork/edmd_gen3.o bin/orig/edmd_gen3.o; then objcmp="identical"; else objcmp="different"; fi
note "1 objects: bin/fork/edmd_gen3.o and bin/orig/edmd_gen3.o are byte-$objcmp (information; a compiler may embed the source path)"

# ------------------------------------------------------------------ test 1b: the eight processes, in parallel
trap 'trap - INT TERM HUP; echo "interrupted: stopping the test processes"; kill 0' INT TERM HUP
launch() {   # launch <tag> <binary> [args...]: own work dir as cwd, capped at 900 s; writes <tag>.out.txt, .err.txt, .rc.txt
  tag=$1; shift
  (
    W="$D/$OUT/work_$tag"; mkdir -p "$W" && cd "$W" || exit 1
    t0=$(now)
    perl -e 'alarm 900; exec { $ARGV[0] } @ARGV or die "exec failed: $!\n"' -- "$@" > "$D/$OUT/$tag.out.txt" 2> "$D/$OUT/$tag.err.txt"
    rc=$?
    t1=$(now)
    echo "$rc $(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b - a}')" > "$D/$OUT/$tag.rc.txt"
  ) &
}
T0=$(now)
for side in orig fork; do
  launch m1f_$side "$D/bin/$side/gen3_m1_harness" audit
  launch m1q_$side "$D/bin/$side/gen3_m1_harness" audit --quick
  launch m2q_$side "$D/bin/$side/gen3_m2_harness" audit --quick
  launch fs_$side  "$D/bin/$side/box_final_state"
done
echo "eight processes started (the longest, M1 audit, takes about 5 min on 4 cores)"
wait
trap - INT TERM HUP
note "the eight processes took $(awk -v a="$T0" -v b="$(now)" 'BEGIN{printf "%.0f", b - a}') s together"

# The verdicts of one output: prints "VERDICTS OK", or the problems found. A harness exits 0 even after RUN FAILED or FATAL,
# so these checks are what catches a failure that both sides share. The expected counts are those of the harness sources
# whose blobs the snapshot records (M1: 6 cells of 3 audit classes; M2: 13 cells of 5); a harness change shows up here.
verdicts() {   # verdicts <m1|m2|fs> <file>
  awk -v kind="$1" '
    function trim(s) { gsub(/^ +| +$/, "", s); return s }
    function zero(s,  n, k, a) { n = split(s, a, /[^0-9.eE+-]+/); for (k = 1; k <= n; ++k) if (a[k] != "" && a[k] + 0 != 0) return 0; return 1 }
    function bad(m) { ++nbad; if (nbad <= 5) print "line " NR ": " m }
    function fail(m) { ++nbad; print m }
    /RUN FAILED|FATAL/ { bad($0) }
    { last = $0 }
    kind != "fs" && /^### Cell / { ++cells }
    kind != "fs" && /^audits do not steer/ {
      want = (kind == "m1") ? "audits do not steer (A = B): YES; deterministic (B = C): YES" : "audits do not steer (A = B): YES; same-seed bit identity (B = C): YES"
      if ($0 == want) ++steer; else bad($0) }
    kind == "m2" && /^run flag \(edmd3_health_clean\)/ { if ($0 == "run flag (edmd3_health_clean): A 1, B 1, C 1") ++flag; else bad($0) }
    kind != "fs" && /^\| [BC] gen3, audits off/ { n = split($0, f, "|"); if (trim(f[n - 1]) == "yes") ++same; else bad($0) }
    kind != "fs" && /must be 0/ { s = $0; gsub(/\((all )?must be 0\)/, "", s); if (!zero(s)) bad($0) }
    kind != "fs" && /^\| (pairs|outer walls|crossings|divider faces|pistons) \|/ {
      n = split($0, f, "|"); ++rows
      if (trim(f[4]) != "0" || trim(f[5]) != "0" || trim(f[7]) != "0") bad("missing, extra or |dt| > 1e-10 x horizon: " $0)
      if (kind == "m2" && trim(f[9]) != "0" && trim(f[9]) != "-") bad("deferred earlier than eligible: " $0) }
    kind != "fs" && /^\[EDMD3-HEALTH\]/ {
      ++health
      for (k = 2; k <= NF; ++k) { if (split($k, kv, "=") != 2) continue
        if (kv[1] ~ /^(overlap_repair|wall_overdue|obj_overlap_repair|past_event|clamp_repair|cell_repair|grid_escape|stagnation|local_findings|full_findings|body_findings)$/ && kv[2] != "0") bad($k) } }
    kind == "m1" && /^max \|dt\| over matched events/ { ++dts; if ($0 !~ /; duplicate disagreements 0; disks outside their cell 0$/) bad($0) }
    kind == "m1" && /^Z gen3 - gen2 = / { ++zl }
    kind == "m1" && /^\| [0-9.]+ \| (yes|no) \|/ { n = split($0, f, "|"); ++wrow; w = trim(f[2]); a = trim(f[3])
      if ((w == "24.0") != (a == "yes")) bad("cell-width start-up check: " $0) }
    kind == "m2" && /^\| gen2 \| / { ++g2 }
    kind == "fs" && /^## cell / { ++cells; idx = 0 }
    kind == "fs" && /^run_flag / { if ($2 == "1") ++flag; else bad($0) }
    kind == "fs" && /^schedule_audit / { ++sa; for (k = 4; k <= NF; ++k) if ($k ~ /^[0-9]+$/ && $k != "0") bad("schedule audit: " $0) }
    kind == "fs" && /^restart_from_checkpoint / { if ($0 ~ /equal hash 1 health 1 ledger 1 disks 1$/) ++rs; else bad($0) }
    kind == "fs" && /^tail / { bad("a record longer than the original: " $0) }
    kind == "fs" && /^disk / { if (NF != 7 || $2 != idx "") bad("disk line " idx ": " $0); ++idx; ++disks }
    END {
      if (kind == "m1") { wc = 6; wr = 3 } else if (kind == "m2") { wc = 13; wr = 5 } else wc = 2
      if (cells != wc) fail("cells " cells ", want " wc)
      if (kind != "fs") {
        if (steer != cells) fail("audits-do-not-steer YES lines " steer " for " cells " cells")
        if (same != 2 * cells) fail("B and C equal to A: " same " of " 2 * cells " rows")
        if (rows != wr * cells) fail("schedule-audit class rows " rows ", want " wr * cells)
        if (health != cells) fail("[EDMD3-HEALTH] lines " health ", want " cells) }
      if (kind == "m1" && (dts != cells || zl != cells || wrow != 4)) fail("m1 per-cell lines: max-dt " dts ", Z " zl ", start-up rows " wrow " (want " cells ", " cells ", 4)")
      if (kind == "m2" && (flag != cells || g2 != 2 * cells)) fail("m2 per-cell lines: run flag " flag ", gen2 rows " g2 " (want " cells ", " 2 * cells ")")
      if (kind == "fs" && (flag != 2 || sa != 2 || rs != 2 || disks != 800 || last != "# end")) fail("fs: run flags " flag ", audits " sa ", restarts " rs ", disks " disks ", last line \"" last "\" (want 2, 2, 2, 800, \"# end\")")
      if (nbad > 5) print "(" nbad - 5 " more problems)"
      if (nbad == 0) print "VERDICTS OK"
    }' "$2"
}
want_lines() { case $1 in m1f|m1q) echo 224;; m2q) echo 749;; fs) echo 833;; esac; }
nchecked=0
for c in m1f m1q m2q fs; do
  case $c in m1f) what="M1 audit"; kind=m1;; m1q) what="M1 audit --quick"; kind=m1;; m2q) what="M2 audit --quick"; kind=m2;; fs) what="final-state driver"; kind=fs;; esac
  set -- $(cat "$OUT/${c}_orig.rc.txt" 2>/dev/null) x x; rco=$1; so=$2
  set -- $(cat "$OUT/${c}_fork.rc.txt" 2>/dev/null) x x; rcf=$1; sf=$2
  if [ "$rco" != 0 ] || [ "$rcf" != 0 ]; then
    bad "1 $what: return codes original $rco, fork $rcf (142 = stopped at 900 s; see $OUT/${c}_orig.err.txt, $OUT/${c}_fork.err.txt)"; continue; fi
  probs=""
  cmp -s "$OUT/${c}_orig.out.txt" "$OUT/${c}_fork.out.txt" || probs="$probs stdout differs (diff $OUT/${c}_orig.out.txt $OUT/${c}_fork.out.txt);"
  for side in orig fork; do
    [ -s "$OUT/${c}_$side.err.txt" ] && probs="$probs $side stderr not empty ($OUT/${c}_$side.err.txt);"
    n=$(nlines "$OUT/${c}_$side.out.txt"); w=$(want_lines $c)
    [ "$n" = "$w" ] || probs="$probs $side has $n lines, want $w;"
    v=$(verdicts $kind "$OUT/${c}_$side.out.txt" 2>&1)
    nchecked=$((nchecked + 1))
    [ "$v" = "VERDICTS OK" ] || probs="$probs $side verdicts: $(echo "$v" | head -4 | tr '\n' ';')"
  done
  if [ -z "$probs" ]; then
    ok "1 $what: fork = original, stdout byte-identical ($(nlines "$OUT/${c}_fork.out.txt") lines, sha256 $(sha "$OUT/${c}_fork.out.txt" | cut -c1-16)), stderr empty, verdicts good on both sides; ${so} s and ${sf} s"
  else bad "1 $what:$probs"; fi
done
note "verdict checks ran on $nchecked of 8 outputs (each must print exactly VERDICTS OK)"

# ------------------------------------------------------------------ the manifest (pre-check A10)
M=$OUT/manifest.txt
{
  echo "# PID Demo 1, delivery 1a: the box reference (brief 02_DEMO1_BRIEF.md sec. 6 tests 1-2; pre-check part A)"
  cat "$OUT/snapshot.txt"
  echo "## the fork alone"
  echo "objects bin/fork/edmd_gen3.o vs bin/orig/edmd_gen3.o: byte-$objcmp"
  echo "nm -u of the fork object ($(nlines "$OUT/nm_u_fork.txt") symbols; equal to the original object's: $(cmp -s "$OUT/nm_u_fork.txt" "$OUT/nm_u_orig.txt" && echo yes || echo no)): $(tr '\n' ' ' < "$OUT/nm_u_fork.txt")"
  echo "stub link with -lm only: rc $lrc"
  echo "## runs: tag | binary and arguments | rc | wall s | stdout lines | stdout sha256 | stderr bytes"
  for c in m1f m1q m2q fs; do for side in orig fork; do
    case $c in m1f) a="gen3_m1_harness audit";; m1q) a="gen3_m1_harness audit --quick";; m2q) a="gen3_m2_harness audit --quick";; fs) a="box_final_state";; esac
    set -- $(cat "$OUT/${c}_$side.rc.txt" 2>/dev/null) x x
    echo "${c}_$side | bin/$side/$a | $1 | $2 | $(nlines "$OUT/${c}_$side.out.txt") | $(sha "$OUT/${c}_$side.out.txt") | $(wc -c < "$OUT/${c}_$side.err.txt" | tr -d ' ')"
  done; done
  echo "## what is compared"
  echo "Per pair of runs: stdout byte for byte (cmp), with no line excluded (both sides use one compiler and one flag set, so"
  echo "even the build line is equal), and stderr must be empty. The harnesses print event hashes, event counts, health and"
  echo "audit counters, contact-audit maxima, KE drift and Z, but no final positions. box_final_state prints every disk's final"
  echo "position, velocity and collision count, the bodies, the impulse, virial and tie accessors, the collision-hook digest,"
  echo "the event-log digest, a checkpoint-restart verdict, and the health, audit, ledger and tolerance records (in %a or raw"
  echo "bytes). Each output also passes its own verdicts (exact cell, row and line counts; every verdict line; the gen3 health"
  echo "counters 0; the last line) so that a failure shared by both sides is not read as agreement."
  echo "Cloud (Linux) and Mac hashes differ for the same sources (pre-check A7): compare only within one machine and compiler."
  echo "## results"
  cat "$LOG"
  echo "summary: $PASS passed, $FAIL failed"
} > "$M"
echo "manifest: $M"
echo "$PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
