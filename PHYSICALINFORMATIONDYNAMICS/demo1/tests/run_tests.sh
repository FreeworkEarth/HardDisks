#!/bin/sh
# tests/run_tests.sh -- Demo 1, deliveries 1a and 1b: the tests of brief 02_DEMO1_BRIEF.md sec. 6 as amended by sec. 8, with
# the corrections of 02_DEMO1_BRIEF_PRECHECK.md. Run through `make test` from demo1/; the Makefile rebuilds bin/ first.
#   Test 2 (the 1b form, pre-check B21) each engine file starts with the exact 4-line header (its copy commit and blob confirmed
#           by git; the "last change" commit only reported, sec. 8); edmd.h is its original byte for byte below the header; every
#           hunk in which edmd_gen3.c or edmd_gen3.h differs from its original adds a line carrying "/* PID: ". The hunks are
#           git's (git diff --no-index, Myers): the system diff cuts them differently (GNU diff on Linux, BSD diff on macOS
#           13 and later), so a hunk that is marked on one machine can be split into an unmarked piece on another.
#   Test 3 (= test 1, with the 1b fork: its walls unused) the fork and the unmodified original gen3, each linked into its own
#           build of the unchanged M1/M2 harness sources (gen2 linked into both as one object) and of tests/box_final_state.c,
#           give byte-identical stdout for M1 audit, M1 audit --quick, M2 audit --quick and the final-state driver (whose
#           fork side prints the bytes 1b appended to the public records on "tail" lines: stripped before the comparison and
#           required to be zero), each a separate process, with empty stderr; each output also passes its own verdict checks
#           (cell counts, verdict lines, health counters, exact line count, end), because two equally broken outputs would
#           compare equal (pre-check A5). The fork's object links into a program with nothing but libc and libm, and an empty
#           wall list changes nothing (prim_tests empty).
#   Tests 4-5 prim_tests pair 0.05 (4000 sigma-time) and pair 0.30 (400): a closed 100 x 40 sigma corridor, straight and
#           turned by 30 degrees; test 6: prim_tests y 0.05, y 0.30, yaim 0, yaim 1 (a Y of segments and arcs; disks aimed at its
#           reflex vertices and at an arc-segment joint, with and without a 1.7e-12 px vertex split); test 10: prim_controls
#           (four planted errors, each reported by its own detector; the fourth, a wall event moved in time, has a grazing
#           disk as its baseline: the same contact to rounding, which the audit must not count as an error). Each of these prints "check ... OK|FAIL" lines; all must
#           be OK, their number complete, the last line "failed checks: 0".
#   Test 7  Demo 0 built headless (pre-check B12) grows --preset=mat45 --rule=guided --seed=9300 --steps=400; its final.txt is
#           identified by structure (sec. 8: 416 alive disks, 415 edges, 10 loops, 833 segments + 330 arcs, lumen_area
#           6659.169369 to 1e-6 relative); pid_demo1 runs N = 2000 for 200 sigma-time in its lumen four times: A with gen3's
#           audits, B and C without (the same seed), D with another seed. Each: rc 0, empty stderr, its checks OK (clean = 1,
#           |KE drift| <= 1e-12, no escape at any of 21 stops; in A also nothing missing or extra and the wall gaps <= tol_face).
#           pid_demo1 refuses a Demo 0 file without [hull].
#   Test 8  B = C and A = B byte for byte in every pstate_*.txt, density_*.txt and measure.csv; D differs from B.
#   Test 9  every density file of B: the particle counts sum to N, the material column to the sum of pi r^2 over the alive
#           disks of final.txt to 1e-9 relative (B14; M_used and the lens total printed beside it), the kinetic-energy column
#           (1c, brief sec. 9) to the KE of measure.csv at the same t to 1e-12 relative.
#   Test 11 (1c) the window, bin/pid_view1, runs test 7's B under SDL's dummy video driver, advancing 0.37 sigma-time per frame
#           (frames between the stops) with a screenshot at t = 100: rc 0, its checks OK, and every pstate_*.txt, density_*.txt
#           and measure.csv byte-identical to the headless B. Without SDL2 (pkg-config sdl2) pid_view1 is not built and the
#           test is reported as not run.
#   Test 12 (studio) bin/pid_studio, the window with panels, under the dummy driver: (a) every Demo 0 preset (0a, 0a1, mat, mat45,
#           mat45r; rule guided, seed 9300, 150 ticks) grown in the studio (--auto --grow-only) writes the state file Demo 0
#           headless writes, byte for byte except its build line (the presets are copied from demo0/pid_main.c); (b) test 7's
#           configuration grown, applied and run in the studio (--auto) gives test 7's final.txt (except the build line) and
#           all 43 data files of B byte for byte. Without SDL2_ttf pid_studio is not built and the test is reported as not run.
# Three batches of processes run in parallel (8, 8 and 6), each process capped at 900 s (perl alarm: macOS has no `timeout`); Ctrl-C
# stops them all. Outputs: a new folder runs/box_reference_<cc>_<UTC time>/ (git-ignored) with manifest.txt. Nothing else is
# written and nothing is deleted. One `make test` at a time per checkout: bin/ is shared.
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
  for f in Makefile tests/run_tests.sh tests/box_final_state.c tests/prim_tests.c tests/prim_controls.c tests/demo0_view_stub.c pid_demo1.h \
           pid_demo1.c pid_draw1.h pid_draw1.c pid_view1.c pid_studio.c tests/bmp2png.py tests/make_shots.sh engine/README_FORK.md \
           engine/edmd_gen3.c engine/edmd_gen3.h engine/edmd.h; do echo "$(sha "$f")  demo1/$f"; done
  echo "## Demo 0 sources compiled unchanged from ../demo0 (git blob)"
  for f in pid.h pid_main.c pid_grow.c pid_state.c pid_rng.c pid_sha256.c pid_hull.c; do echo "$(git hash-object "../demo0/$f")  demo0/$f"; done
  echo "## originals (git blob)"
  for f in edmd_gen3.c edmd_gen3.h edmd.h edmd.c tests/gen3_m1_harness.c tests/gen3_m2_harness.c tests/gen3_state_io.h; do echo "$(git hash-object "$HSP/$f")  hspist3/edmd_core/$f"; done
  echo "## the tested objects and programs (sha256)"
  for f in bin/gen2.o bin/orig/edmd_gen3.o bin/fork/edmd_gen3.o bin/orig/gen3_m1_harness bin/fork/gen3_m1_harness bin/orig/gen3_m2_harness \
           bin/fork/gen3_m2_harness bin/orig/box_final_state bin/fork/box_final_state bin/demo0/pid_rng.o bin/fork/prim_tests bin/fork/prim_controls \
           bin/pid_demo1 bin/demo0/pid_demo0_headless bin/pid_view1 bin/pid_studio; do
    if [ -f "$f" ]; then echo "$(sha "$f")  demo1/$f"; else echo "(not built)  demo1/$f"; fi; done
  echo "sdl2 $(pkg-config --modversion sdl2 2>/dev/null || echo 'not found by pkg-config'); SDL2_ttf $(pkg-config --modversion SDL2_ttf 2>/dev/null || echo 'not found by pkg-config')"
  echo "## build: every command make runs for 'all' (make -n -B), in order"
  cat bin/build_lines.txt 2>/dev/null || echo "(bin/build_lines.txt missing: run through make test)"
} > "$OUT/snapshot.txt"

# ------------------------------------------------------------------ test 2 (1b form): the fork's sources
t2=0
for f in edmd_gen3.c edmd_gen3.h edmd.h; do
  fk=engine/$f; or=$HSP/$f
  l2=$(sed -n 2p "$fk"); l3=$(sed -n 3p "$fk")
  cpy=$(echo "$l2" | sed -n 's/^   copied at origin\/main \([0-9a-f]\{40\}\); last change to it \([0-9a-f]\{40\}\);$/\1/p')
  lst=$(echo "$l2" | sed -n 's/^   copied at origin\/main \([0-9a-f]\{40\}\); last change to it \([0-9a-f]\{40\}\);$/\2/p')
  blob=$(echo "$l3" | sed -n 's/^   source blob \([0-9a-f]\{40\}\)\. Lines 5 on are that file unchanged: its line n is line n + 4 here\.$/\1/p')
  {
    echo "/* PID fork, PHYSICALINFORMATIONDYNAMICS/demo1/engine (Demo 1, delivery 1a, 2026-10-16): hspist3/edmd_core/$f"
    echo "   copied at origin/main $cpy; last change to it $lst;"
    echo "   source blob $blob. Lines 5 on are that file unchanged: its line n is line n + 4 here."
    echo '   Changes from delivery 1b on are marked at their site by a comment starting "PID:" (README_FORK.md). */'
  } > "$OUT/header_$f.txt"
  if [ -z "$cpy" ] || [ -z "$blob" ]; then
    bad "2 source: $fk lines 2-3 do not record the copy commit, the last change and the source blob"; t2=1; continue
  elif ! head -n 4 "$fk" | cmp -s - "$OUT/header_$f.txt"; then
    bad "2 source: $fk lines 1-4 are not the header template (diff $OUT/header_$f.txt)"; t2=1; continue
  elif [ "$(git rev-parse "$cpy:hspist3/edmd_core/$f" 2>/dev/null)" != "$blob" ]; then
    bad "2 source: at the recorded copy commit $cpy, hspist3/edmd_core/$f is not the recorded blob $blob"; t2=1; continue
  elif [ "$(git hash-object "$or")" != "$blob" ]; then
    bad "2 source: $or is now blob $(git hash-object "$or"), the fork records $blob (the original changed after the fork)"; t2=1; continue
  fi
  lgit=$(git log -1 --format=%H "$cpy" -- "$or" 2>/dev/null)
  [ "$lgit" = "$lst" ] || note "2 source: git names $lgit as the last change to hspist3/edmd_core/$f before the copy, the header $lst (clone depth and merge history decide it: provenance, not identity; sec. 8)"
  tail -n +5 "$fk" > "$OUT/body_$f"
  if [ "$f" = edmd.h ]; then
    cmp -s "$OUT/body_$f" "$or" || { bad "2 source: engine/edmd.h below its header differs from $or (it must stay byte-identical)"; t2=1; }
    continue
  fi
  git diff --no-index --no-color --no-ext-diff --no-textconv -U0 --diff-algorithm=myers --indent-heuristic -- "$or" "$OUT/body_$f" \
    > "$OUT/diff_$f.txt"
  hk=$(awk '/^@@/{ if (n) { if (!m) bad++; } n++; m = 0; next } /^\+\+\+|^---/ { next } /^\+/ { if (index($0, "/* PID: ")) m = 1 }
           END { if (n) { if (!m) bad++ }; print n + 0, bad + 0 }' "$OUT/diff_$f.txt")
  set -- $hk
  if [ "$2" -ne 0 ]; then bad "2 source: $2 of the $1 changed hunks of engine/$f add no line with \"/* PID: \" ($OUT/diff_$f.txt)"; t2=1
  else note "2 source: engine/$f differs from its original in $1 hunks, each adding a \"/* PID: \" line"; fi
done
extra=$(ls -A engine | grep -v -x -e edmd_gen3.c -e edmd_gen3.h -e edmd.h -e README_FORK.md)
if [ -n "$extra" ]; then bad "2 source: engine/ holds other entries: $(echo $extra)"; t2=1; fi
if [ ! -f engine/README_FORK.md ]; then bad "2 source: engine/README_FORK.md is missing"; t2=1; fi
if [ "$t2" -eq 0 ]; then ok "2 source: the exact 4-line headers (copy commits and blobs confirmed by git), edmd.h unchanged below it, every changed hunk of edmd_gen3.c and .h marked /* PID: */"; fi

# ------------------------------------------------------------------ test 1a: the fork alone
nmok=1
nm -u bin/fork/edmd_gen3.o > "$OUT/nm_raw_fork.txt" 2>&1 || nmok=0
nm -u bin/orig/edmd_gen3.o > "$OUT/nm_raw_orig.txt" 2>&1 || nmok=0
awk '{print $NF}' "$OUT/nm_raw_fork.txt" | sort > "$OUT/nm_u_fork.txt"
awk '{print $NF}' "$OUT/nm_raw_orig.txt" | sort > "$OUT/nm_u_orig.txt"
printf 'int main(void){ return 0; }\n' > "$OUT/stub_main.c"
$CC $CSTD -o "$OUT/stub_link" "$OUT/stub_main.c" bin/fork/edmd_gen3.o -lm > "$OUT/stub_link.txt" 2>&1; lrc=$?
if [ "$lrc" -eq 0 ] && "$OUT/stub_link" && [ "$nmok" -eq 1 ] && [ -s "$OUT/nm_u_fork.txt" ]; then
  ok "3 fork alone: the fork's object links into a stub program with libc and -lm only ($(nlines "$OUT/nm_u_fork.txt") undefined symbols)"
else bad "3 fork alone: stub link rc $lrc ($OUT/stub_link.txt), nm ok $nmok, or nm -u empty"; fi
cmp -s "$OUT/nm_u_fork.txt" "$OUT/nm_u_orig.txt" && nmeq=yes || nmeq="no ($(comm -3 "$OUT/nm_u_orig.txt" "$OUT/nm_u_fork.txt" | tr -d '\t' | tr '\n' ' '))"
note "3 the fork's undefined symbols equal the original object's: $nmeq (information: 1b may call more of libm)"
if cmp -s bin/fork/edmd_gen3.o bin/orig/edmd_gen3.o; then objcmp="identical"; else objcmp="different"; fi
note "3 objects: bin/fork/edmd_gen3.o and bin/orig/edmd_gen3.o are byte-$objcmp (information; 1b changes the fork)"

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
echo "batch 1: eight processes started (the longest, M1 audit, takes about 5 min on 4 cores)"
wait
note "batch 1 (test 3) took $(awk -v a="$T0" -v b="$(now)" 'BEGIN{printf "%.0f", b - a}') s"
T1=$(now)
P="$D/bin/fork/prim_tests"
launch t45a "$P" pair 0.05 4000 0x5EED30
launch t45b "$P" pair 0.30 400 0x5EED31
launch t6a  "$P" y 0.05 400 0x5EED32
launch t6b  "$P" y 0.30 400 0x5EED33
launch t6c  "$P" yaim 0 50 0x5EED34
launch t6d  "$P" yaim 1 50 0x5EED34
launch t3e  "$P" empty 100 0x5EED20
launch t10  "$D/bin/fork/prim_controls"
echo "batch 2: eight processes started (tests 4-6, 10 and the empty wall list; about 2 min)"
wait
trap - INT TERM HUP
note "batch 2 (tests 4-6, 10) took $(awk -v a="$T1" -v b="$(now)" 'BEGIN{printf "%.0f", b - a}') s"

# ------------------------------------------------------------------ tests 7-9: a real hull
T2=$(now)
H7="$D/$OUT/t7_demo0"
"$D/bin/demo0/pid_demo0_headless" --preset=mat45 --rule=guided --seed=9300 --steps=400 --headless --out="$H7" > "$OUT/t7_demo0.out.txt" 2> "$OUT/t7_demo0.err.txt"
d0rc=$?
F7="$H7/final.txt"
trap 'trap - INT TERM HUP; echo "interrupted: stopping the test processes"; kill 0' INT TERM HUP
P1="$D/bin/pid_demo1"
launch t7a "$P1" --state="$F7" --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --audit-every=20000 --out="$D/$OUT/t7a"
launch t7b "$P1" --state="$F7" --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --out="$D/$OUT/t7b"
launch t7c "$P1" --state="$F7" --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --out="$D/$OUT/t7c"
launch t7d "$P1" --state="$F7" --N=2000 --T=200 --seed=0x5EED71 --write-every=10 --out="$D/$OUT/t7d"
launch t7n "$P1" --state="$H7/state_000350.txt" --N=10 --T=1 --seed=1 --out="$D/$OUT/t7n"
if [ -x "$D/bin/pid_view1" ]; then
  launch t11 env SDL_VIDEODRIVER=dummy "$D/bin/pid_view1" --state="$F7" --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --exit-at-end \
    --speed=0.37 --shot-at=100 --out="$D/$OUT/t11"
fi
if [ -x "$D/bin/pid_studio" ]; then
  launch t12r env SDL_VIDEODRIVER=dummy "$D/bin/pid_studio" --preset=mat45 --rule=guided --grow-seed=9300 --steps=400 --tpf=4 --tick=400 \
    --N=2000 --seed=0x5EED70 --T=200 --write-every=10 --speed=0.37 --auto --exit-at-end --shot-tick=200 --shot-at=100 --out="$D/$OUT/t12r"
fi
echo "batch 3: Demo 0 headless (rc $d0rc), then the pid_demo1 processes (and the window) started (tests 7-9, 11; under a minute)"
wait
trap - INT TERM HUP
note "batch 3 (tests 7-9) took $(awk -v a="$T2" -v b="$(now)" 'BEGIN{printf "%.0f", b - a}') s"

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
    kind == "fs" && /^tail / { if ($NF ~ /[1-9a-f]/) bad("bytes appended to a record are not zero: " $0) }
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
    bad "3 $what: return codes original $rco, fork $rcf (142 = stopped at 900 s; see $OUT/${c}_orig.err.txt, $OUT/${c}_fork.err.txt)"; continue; fi
  probs=""
  for side in orig fork; do grep -v '^tail ' "$OUT/${c}_$side.out.txt" > "$OUT/${c}_$side.cmp.txt"; done   # 1b appends to the records
  cmp -s "$OUT/${c}_orig.cmp.txt" "$OUT/${c}_fork.cmp.txt" || probs="$probs stdout differs (diff $OUT/${c}_orig.cmp.txt $OUT/${c}_fork.cmp.txt);"
  [ "$c" = fs ] && [ "$(grep -c '^tail ' "$OUT/${c}_orig.out.txt")" -ne 0 ] && probs="$probs the original prints tail lines;"
  for side in orig fork; do
    [ -s "$OUT/${c}_$side.err.txt" ] && probs="$probs $side stderr not empty ($OUT/${c}_$side.err.txt);"
    n=$(nlines "$OUT/${c}_$side.cmp.txt"); w=$(want_lines $c)
    [ "$n" = "$w" ] || probs="$probs $side has $n lines, want $w;"
    v=$(verdicts $kind "$OUT/${c}_$side.out.txt" 2>&1)
    nchecked=$((nchecked + 1))
    [ "$v" = "VERDICTS OK" ] || probs="$probs $side verdicts: $(echo "$v" | head -4 | tr '\n' ';')"
  done
  if [ -z "$probs" ]; then
    tl=$(grep -c '^tail ' "$OUT/${c}_fork.out.txt")
    ok "3 $what: fork = original, stdout byte-identical ($(nlines "$OUT/${c}_fork.cmp.txt") lines, sha256 $(sha "$OUT/${c}_fork.cmp.txt" | cut -c1-16))$([ "$tl" -gt 0 ] && echo ", the fork's $tl tail lines zero"), stderr empty, verdicts good on both sides; ${so} s and ${sf} s"
  else bad "3 $what:$probs"; fi
  set -- $(cat "$OUT/${c}_orig.rc.txt" 2>/dev/null) x x
done
note "verdict checks ran on $nchecked of 8 outputs (each must print exactly VERDICTS OK)"

# tests 3 (empty wall list), 4-6 and 10: rc 0, empty stderr, every check OK and their number complete, the last line
for c in t3e t45a t45b t6a t6b t6c t6d t10; do
  case $c in t3e) what="3 empty wall list"; want=1;; t45a) what="4-5 corridor, straight and turned 30 degrees, eta 0.05"; want=19;;
             t45b) what="4-5 corridor, straight and turned 30 degrees, eta 0.30"; want=15;; t6a) what="6 Y, eta 0.05"; want=7;;
             t6b) what="6 Y, eta 0.30"; want=7;; t6c) what="6(d) Y, aimed disks"; want=8;; t6d) what="6(d) Y, aimed disks, vertex split 1.7e-12 px"; want=8;;
             t10) what="10 positive controls"; want=8;; esac
  f="$OUT/$c.out.txt"; set -- $(cat "$OUT/$c.rc.txt" 2>/dev/null) x x; rc=$1; tw=$2
  nok=$(grep -c '^check .* OK$' "$f"); nfl=$(grep -c '^check .* FAIL$' "$f"); last=$(tail -n 1 "$f")
  if [ "$rc" = 0 ] && [ ! -s "$OUT/$c.err.txt" ] && [ "$nfl" -eq 0 ] && [ "$nok" -eq "$want" ] && [ "$last" = "failed checks: 0" ]; then
    ok "$what: $nok of $nok checks OK; ${tw} s"
  else bad "$what: rc $rc, $nok OK of $want, $nfl FAIL, last line '$last' ($f, $OUT/$c.err.txt)"; fi
done

# test 7: the hull, by structure (sec. 8, B13), from final.txt itself
st=$(awk '/^\[/ { sec = $0; next } /^#/ { next }
  sec == "[disks]" { alive[$1] = $8; par[$1] = $5; n++ }
  sec == "[hull]" && /^loop / { nl++ } sec == "[hull]" && /^seg / { ns++ } sec == "[hull]" && /^arc / { na++ }
  sec == "[areas]" && $1 == "lumen_area" { la = $2 }
  END { for (i in alive) if (alive[i] == 1) { a++; if (par[i] >= 0 && alive[par[i]] == 1) e++ }
        printf "%d %d %d %d %d %.9f %d", a, e, nl, ns, na, la, (la > 0 && (la - 6659.169369 < 0 ? 6659.169369 - la : la - 6659.169369) <= 1e-6 * 6659.169369) }' "$F7" 2>/dev/null)
set -- $st x x x x x x x
if [ "$d0rc" = 0 ] && [ ! -s "$OUT/t7_demo0.err.txt" ] && [ "$1 $2 $3 $4 $5 $7" = "416 415 10 833 330 1" ]; then
  ok "7 hull: Demo 0 headless (mat45, guided, seed 9300, 400 ticks) gives the sec. 8 hull: $1 alive disks, $2 edges, $3 loops, $4 segments + $5 arcs, lumen_area $6"
else bad "7 hull: Demo 0 rc $d0rc ($OUT/t7_demo0.err.txt); structure '$st', want 416 415 10 833 330 6659.169369 1"; fi
for c in t7a t7b t7c t7d; do
  case $c in t7a) want=5; what="A, gen3's audits on";; t7b) want=3; what="B";; t7c) want=3; what="C (= B)";; t7d) want=3; what="D, another seed";; esac
  f="$OUT/$c.out.txt"; set -- $(cat "$OUT/$c.rc.txt" 2>/dev/null) x x; rc=$1; tw=$2
  nok=$(grep -c '^check .* OK$' "$f"); nfl=$(grep -c '^check .* FAIL$' "$f"); last=$(tail -n 1 "$f")
  if [ "$rc" = 0 ] && [ ! -s "$OUT/$c.err.txt" ] && [ "$nfl" -eq 0 ] && [ "$nok" -eq "$want" ] && [ "$last" = "failed checks: 0" ] \
     && [ "$(ls "$OUT/$c" | grep -c '^pstate_')" -eq 21 ]; then
    ok "7 real hull, N 2000, 200 sigma-time, $what: $nok of $nok checks OK (clean 1, KE drift, no escape at 21 stops$([ $c = t7a ] && echo ', audits')); ${tw} s"
  else bad "7 real hull $what: rc $rc, $nok OK of $want, $nfl FAIL, last line '$last', $(ls "$OUT/$c" 2>/dev/null | grep -c '^pstate_') pstate files ($f)"; fi
done
set -- $(cat "$OUT/t7n.rc.txt" 2>/dev/null) x x
if [ "$1" = 2 ] && grep -q 'has no \[hull\]' "$OUT/t7n.err.txt"; then ok "7 refusal: pid_demo1 refuses a Demo 0 file without [hull] (state_000350.txt; rc 2)"
else bad "7 refusal: pid_demo1 on a file without [hull]: rc $1 ($OUT/t7n.err.txt)"; fi
# test 8: determinism, and the audits do not steer
same_run() {   # same_run <dir1> <dir2>: the same data files, byte for byte (manifest.txt aside: it records argv)
  [ "$(ls "$1" | grep -v -x manifest.txt)" = "$(ls "$2" | grep -v -x manifest.txt)" ] || return 1
  for g in $(ls "$1" | grep -v -x manifest.txt); do cmp -s "$1/$g" "$2/$g" || return 1; done
  return 0
}
nd=$(ls "$OUT/t7b" 2>/dev/null | grep -v -x manifest.txt | wc -l | tr -d ' ')
if same_run "$OUT/t7b" "$OUT/t7c" && same_run "$OUT/t7a" "$OUT/t7b" && [ "$nd" -eq 43 ] \
   && ! cmp -s "$OUT/t7b/pstate_000200.txt" "$OUT/t7d/pstate_000200.txt" && ! cmp -s "$OUT/t7b/measure.csv" "$OUT/t7d/measure.csv"; then
  ok "8 determinism: B = C and A = B in all $nd data files (21 pstate, 21 density, measure.csv; sha256 of B's measure.csv $(sha "$OUT/t7b/measure.csv" | cut -c1-16)); D (another seed) differs"
else bad "8 determinism: B = C, A = B or D != B fails ($OUT/t7a, t7b, t7c, t7d; $nd data files in B, want 43)"; fi
# test 9: the density files of B
t9=$(awk 'BEGIN { pi = atan2(0, -1) }
  FNR == 1 { fi++ }
  fi == 1 && /^\[/ { sec = $0; next } fi == 1 && /^#/ { next }
  fi == 1 && sec == "[header]" && $1 == "m_used" { mu = $2 }
  fi == 1 && sec == "[disks]" && $8 == 1 { s += pi * $4 * $4 }
  fi == 2 { if (FNR > 1) { split($0, q, ","); ke[q[1]] = q[2] + 0 }; next }
  fi > 2 && /^#/ { next }
  fi > 2 { if (NF != 5) badl++; cnt[FILENAME] += $3; mat[FILENAME] += $4; kc[FILENAME] += $5; files[FILENAME] = 1 }
  END { nf = 0; worst = 0; badc = 0; wk = 0; nk = 0
        for (f in files) { nf++; if (cnt[f] != 2000) badc++; r = (mat[f] - s) / s; if (r < 0) r = -r; if (r > worst) worst = r
          t = f; sub(/.*density_/, "", t); sub(/[.]txt$/, "", t)
          if (!(t in ke) || !(ke[t] > 0)) nk++; else { rk = (kc[f] - ke[t]) / ke[t]; if (rk < 0) rk = -rk; if (rk > wk) wk = rk } }
        printf "%d %d %d %.3e %.9f %.9f %.9f %.3e %d", nf, badc, badl + 0, worst, s, mu, s - mu, wk, nk }' "$F7" "$OUT/t7b/measure.csv" "$OUT"/t7b/density_*.txt 2>/dev/null)
set -- $t9 x x x x x x x x x
if [ "$1" = 21 ] && [ "$2" = 0 ] && [ "$3" = 0 ] && [ "$9" = 0 ] && awk -v w="$4" -v k="$8" 'BEGIN { exit !(w + 0 <= 1e-9 && k + 0 <= 1e-12) }'; then
  ok "9 density: $1 files, particles sum to N = 2000 in each; material sums to sum pi r^2 = $5 sigma^2 over the alive disks (worst relative difference $4); M_used $6, lens total $7 sigma^2; the ke column sums to the KE of measure.csv (worst relative difference $8)"
else bad "9 density: '$t9' (files, count mismatches, bad lines, worst relative material difference, sum pi r^2, M_used, lens, worst relative KE difference, stops without a KE)"; fi
# test 11 (1c): the window writes what the headless run writes
if [ ! -x bin/pid_view1 ]; then note "11 window = headless: not run (bin/pid_view1 is not built: SDL2 not found by pkg-config)"
else
  f="$OUT/t11.out.txt"; set -- $(cat "$OUT/t11.rc.txt" 2>/dev/null) x x; rc=$1; tw=$2
  nok=$(grep -c '^check .* OK$' "$f"); nfl=$(grep -c '^check .* FAIL$' "$f"); last=$(tail -n 1 "$f")
  nd=0; nsame=0
  for g in $(ls "$OUT/t7b" 2>/dev/null | grep -v -x manifest.txt); do nd=$((nd + 1)); cmp -s "$OUT/t7b/$g" "$OUT/t11/$g" && nsame=$((nsame + 1)); done
  if [ "$rc" = 0 ] && [ ! -s "$OUT/t11.err.txt" ] && [ "$nfl" -eq 0 ] && [ "$nok" -eq 3 ] && [ "$last" = "failed checks: 0" ] && [ "$nd" -eq 43 ] \
     && [ "$nsame" -eq 43 ] && [ -s "$OUT/t11/shot_000100.bmp" ]; then
    ok "11 window = headless: bin/pid_view1 under SDL_VIDEODRIVER=dummy, 0.37 sigma-time per frame, writes all 43 data files of B byte for byte (3 of 3 checks OK, a screenshot at t = 100); ${tw} s"
  else bad "11 window = headless: rc $rc, $nok OK of 3, $nfl FAIL, last line '$last', $nsame of $nd data files equal to B's, shot $(ls "$OUT/t11" 2>/dev/null | grep -c '^shot_') ($f, $OUT/t11.err.txt)"; fi
fi

# test 12 (studio): growth = Demo 0 for every preset; grow, apply, run = test 7
if [ ! -x bin/pid_studio ]; then note "12 studio: not run (bin/pid_studio is not built: SDL2_ttf not found by pkg-config)"
else
  n12=0; e12=0; l12=""
  for p in 0a 0a1 mat mat45 mat45r; do
    n12=$((n12 + 1))
    bin/demo0/pid_demo0_headless --preset=$p --rule=guided --seed=9300 --steps=150 --headless --out="$OUT/t12_d0_$p" > "$OUT/t12_d0_$p.out.txt" 2>&1
    SDL_VIDEODRIVER=dummy perl -e 'alarm 900; exec { $ARGV[0] } @ARGV' -- bin/pid_studio --preset=$p --rule=guided --grow-seed=9300 --steps=150 --tpf=7 \
      --auto --grow-only --exit-at-end --out="$OUT/t12_st_$p" > "$OUT/t12_st_$p.out.txt" 2>&1
    a="$OUT/t12_st_$p/net_1_000150/final.txt"; b="$OUT/t12_d0_$p/final.txt"
    if [ -s "$a" ] && [ -s "$b" ] && grep -v '^build ' "$a" > "$OUT/t12_a.txt" && grep -v '^build ' "$b" > "$OUT/t12_b.txt" && cmp -s "$OUT/t12_a.txt" "$OUT/t12_b.txt"; then
      e12=$((e12 + 1)); l12="$l12 $p $(grep -c '^seg ' "$a") seg"
    else l12="$l12 $p DIFFERENT"; fi
  done
  f="$OUT/t12r.out.txt"; set -- $(cat "$OUT/t12r.rc.txt" 2>/dev/null) x x; rc=$1; tw=$2
  nok=$(grep -c '^check .* OK$' "$f"); nfl=$(grep -c '^check .* FAIL$' "$f")
  nd=0; nsame=0
  for g in $(ls "$OUT/t7b" 2>/dev/null | grep -v -x manifest.txt); do nd=$((nd + 1)); cmp -s "$OUT/t7b/$g" "$OUT/t12r/run_1/$g" && nsame=$((nsame + 1)); done
  grep -v '^build ' "$OUT/t12r/net_1_000400/final.txt" > "$OUT/t12_a.txt" 2>/dev/null; grep -v '^build ' "$F7" > "$OUT/t12_b.txt"
  if cmp -s "$OUT/t12_a.txt" "$OUT/t12_b.txt"; then neq=1; else neq=0; fi
  if [ "$e12" -eq 5 ] && [ "$rc" = 0 ] && [ ! -s "$OUT/t12r.err.txt" ] && [ "$nfl" -eq 0 ] && [ "$nok" -eq 3 ] && [ "$neq" = 1 ] && [ "$nd" -eq 43 ] \
     && [ "$nsame" -eq 43 ] && [ -s "$OUT/t12r/shot_grow_tick000200.bmp" ] && [ -s "$OUT/t12r/shot_run1_t000100.bmp" ]; then
    ok "12 studio: growth = Demo 0 headless for all 5 presets (except the build line:$l12); test 7 grown, applied and run in the studio gives its final.txt and all 43 data files of B; ${tw} s"
  else bad "12 studio: presets equal $e12 of 5 ($l12); run: rc $rc, $nok OK of 3, $nfl FAIL, network equal $neq, $nsame of $nd data files equal to B's ($f, $OUT/t12r.err.txt)"; fi
fi

# ------------------------------------------------------------------ the manifest (pre-check A10)
M=$OUT/manifest.txt
{
  echo "# PID Demo 1, deliveries 1a-1b: the box reference and the static walls (brief 02_DEMO1_BRIEF.md sec. 6 and 8; pre-check)"
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
  for c in t3e t45a t45b t6a t6b t6c t6d t10; do
    set -- $(cat "$OUT/$c.rc.txt" 2>/dev/null) x x
    echo "$c | bin/fork/$([ $c = t10 ] && echo prim_controls || echo prim_tests) | $1 | $2 | $(nlines "$OUT/$c.out.txt") | $(sha "$OUT/$c.out.txt") | $(wc -c < "$OUT/$c.err.txt" | tr -d ' ')"
  done
  for c in t7a t7b t7c t7d t7n $([ -x bin/pid_view1 ] && echo t11) $([ -x bin/pid_studio ] && echo t12r); do
    set -- $(cat "$OUT/$c.rc.txt" 2>/dev/null) x x
    case $c in t11) prog="bin/pid_view1, SDL_VIDEODRIVER=dummy (argv in t11/manifest.txt)";; t12r) prog="bin/pid_studio --auto, SDL_VIDEODRIVER=dummy (run argv in t12r/run_1/manifest.txt)";; t7n) prog="bin/pid_demo1 --state=t7_demo0/state_000350.txt --N=10 --T=1 --seed=1 (refused: no manifest)";;
               *) prog="bin/pid_demo1 (argv in $c/manifest.txt)";; esac
    echo "$c | $prog | $1 | $2 | $(nlines "$OUT/$c.out.txt") | $(sha "$OUT/$c.out.txt") | $(wc -c < "$OUT/$c.err.txt" | tr -d ' ')"
  done
  echo "t7_demo0 | bin/demo0/pid_demo0_headless --preset=mat45 --rule=guided --seed=9300 --steps=400 --headless | $d0rc | - | - | final.txt $(sha "$F7" 2>/dev/null) | $(wc -c < "$OUT/t7_demo0.err.txt" | tr -d ' ')"
  echo "## tests 4-6 and 10: the measured values (the full outputs are in this run folder)"
  grep -h -e 'pressure on the walls' -e '^pair eta' -e 'contact audit' -e 'KE drift' -e 'schedule audit (A)' -e 'smallest signed' -e '^control' "$OUT"/t45a.out.txt \
    "$OUT"/t45b.out.txt "$OUT"/t6?.out.txt "$OUT"/t10.out.txt 2>/dev/null
  echo "## test 7: the measured values of A and B"
  grep -h -e '^demo0 ' -e '^seeding' -e '^events' -e '^health' -e '^relative KE' -e '^schedule audit' -e '^contact audit' "$OUT"/t7a.out.txt "$OUT"/t7b.out.txt 2>/dev/null
  echo "## what is compared"
  echo "Per pair of runs: stdout byte for byte (cmp), with no line excluded (both sides use one compiler and one flag set, so"
  echo "even the build line is equal), and stderr must be empty. The harnesses print event hashes, event counts, health and"
  echo "audit counters, contact-audit maxima, KE drift and Z, but no final positions. box_final_state prints every disk's final"
  echo "position, velocity and collision count, the bodies, the impulse, virial and tie accessors, the collision-hook digest,"
  echo "the event-log digest, a checkpoint-restart verdict, and the health, audit, ledger and tolerance records (in %a or raw"
  echo "bytes). Each output also passes its own verdicts (exact cell, row and line counts; every verdict line; the gen3 health"
  echo "counters 0; the last line) so that a failure shared by both sides is not read as agreement."
  echo "Tests 4-6 and 10 (prim_tests, prim_controls) and 7 (pid_demo1) judge their own runs and print one \"check\" line per criterion;"
  echo "the outputs are kept here in full. Test 8 compares the data files of runs A-D of test 7 byte for byte; test 9 sums the"
  echo "density files of B against final.txt independently (awk)."
  echo "Cloud (Linux) and Mac hashes differ for the same sources (pre-check A7): compare only within one machine and compiler."
  echo "## results"
  cat "$LOG"
  echo "summary: $PASS passed, $FAIL failed"
} > "$M"
echo "manifest: $M"
echo "$PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
