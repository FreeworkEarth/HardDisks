#!/bin/sh
# tests/make_shots.sh -- Demo 1c: the two PNGs of brief 02_DEMO1_BRIEF.md sec. 9, test 7's configuration at t = 100 sigma-time
# (Demo 0 headless: --preset=mat45 --rule=guided --seed=9300 --steps=400; pid_view1: N 2000, T 200, seed 0x5EED70), drawn by
# bin/pid_view1 under SDL's dummy video driver (no display needed) and converted by tests/bmp2png.py (Python standard library):
#   shots/demo1c_mat45_seed9300_t100_network_particles.png   --view=network: the hull, the graph, the food patches, the particles
#   shots/demo1c_mat45_seed9300_t100_ke_overlay.png          --overlay=ke --show=hull,food: the particles' kinetic energy per 2 sigma
#                                                             cell (top colour 3 x the mean over the non-empty cells)
#   shots/demo1c_studio_mat45r_tick600.png                   bin/pid_studio (where SDL2_ttf is installed): the panels and a mat45r
#                                                             network growing, at tick 600 (hull and disks)
#   shots/demo1c_shots_manifest.txt                          the commands, the hashes, the run checks
# Run through `make shots` from demo1/. The runs go to runs/shots_<cc>_<UTC time>/ (git-ignored); only shots/ is written besides.
set -u
: "${CC:?run through make shots}" "${IS_CLANG:?run through make shots}"
if [ "$IS_CLANG" = 1 ]; then CCID=clang; else CCID=gcc; fi
[ -x bin/pid_view1 ] || { echo "bin/pid_view1 is not built (SDL2 not found by pkg-config)"; exit 1; }
sha() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -c1-64; else sha256sum "$1" | cut -c1-64; fi; }
OUT=runs/shots_${CCID}_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p runs && mkdir "$OUT" || exit 1
mkdir -p shots
bin/demo0/pid_demo0_headless --preset=mat45 --rule=guided --seed=9300 --steps=400 --headless --out="$OUT/demo0" > "$OUT/demo0.out.txt" 2>&1 \
  || { echo "Demo 0 failed ($OUT/demo0.out.txt)"; exit 1; }
F7="$OUT/demo0/final.txt"
COMMON="--state=$F7 --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --exit-at-end --speed=1 --shot-at=100"
fail=0
for v in network_particles ke_overlay; do
  case $v in network_particles) X="--view=network";; ke_overlay) X="--overlay=ke --show=hull,food";; esac
  # shellcheck disable=SC2086
  SDL_VIDEODRIVER=dummy bin/pid_view1 $COMMON $X --out="$OUT/$v" > "$OUT/$v.out.txt" 2>&1
  rc=$?
  if [ "$rc" != 0 ] || [ "$(tail -n 1 "$OUT/$v.out.txt")" != "failed checks: 0" ] || [ ! -s "$OUT/$v/shot_000100.bmp" ]; then
    echo "pid_view1 $v: rc $rc ($OUT/$v.out.txt)"; fail=1; continue; fi
  python3 -I tests/bmp2png.py "$OUT/$v/shot_000100.bmp" "shots/demo1c_mat45_seed9300_t100_$v.png" || fail=1
done
if [ -x bin/pid_studio ]; then
  SDL_VIDEODRIVER=dummy bin/pid_studio --preset=mat45r --rule=guided --grow-seed=9300 --steps=1000 --tpf=20 --show=disks,hull,food,particles \
    --alpha=0.35 --shot-tick=600 --auto --grow-only --exit-at-end --out="$OUT/studio" > "$OUT/studio.out.txt" 2>&1 \
    && python3 -I tests/bmp2png.py "$OUT/studio/shot_grow_tick000600.bmp" shots/demo1c_studio_mat45r_tick600.png || fail=1
fi
[ "$fail" = 0 ] || exit 1
M=shots/demo1c_shots_manifest.txt
{
  echo "# PID Demo 1c: the two screenshots of 02_DEMO1_BRIEF.md sec. 9 (make shots; tests/make_shots.sh)"
  echo "date_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host $(uname -srm); compiler CC=$CC: $($CC --version 2>/dev/null | head -1); sdl2 $(pkg-config --modversion sdl2 2>/dev/null)"
  echo "repo_head $(git rev-parse HEAD 2>/dev/null); demo1/ entries not committed: $(git status --porcelain --untracked-files=all -- . 2>/dev/null | wc -l | tr -d ' ')"
  echo "demo0: bin/demo0/pid_demo0_headless --preset=mat45 --rule=guided --seed=9300 --steps=400 --headless; final.txt sha256 $(sha "$F7")"
  echo "     ($(grep -c '^seg ' "$F7") segments, $(grep -c '^arc ' "$F7") arcs, $(grep -c '^loop ' "$F7") loops; $(grep '^lumen_area' "$F7"))"
  echo "window: SDL_VIDEODRIVER=dummy bin/pid_view1 --state=<final.txt> --N=2000 --T=200 --seed=0x5EED70 --write-every=10 --exit-at-end --speed=1 --shot-at=100"
  echo "        plus --view=network (network_particles) or --overlay=ke --show=hull,food (ke_overlay); BMP -> PNG by tests/bmp2png.py"
  echo "bin/pid_view1 sha256 $(sha bin/pid_view1)"
  for v in network_particles ke_overlay; do
    echo "## $v"
    grep -h -e '^shot ' -e '^check ' -e '^failed checks' "$OUT/$v.out.txt" | sed "s|$OUT/||"
    echo "pstate_000100.txt sha256 $(sha "$OUT/$v/pstate_000100.txt"); shot_000100.bmp sha256 $(sha "$OUT/$v/shot_000100.bmp")"
    echo "$(sha "shots/demo1c_mat45_seed9300_t100_$v.png")  demo1c_mat45_seed9300_t100_$v.png ($(wc -c < "shots/demo1c_mat45_seed9300_t100_$v.png" | tr -d ' ') bytes)"
  done
  if [ -x bin/pid_studio ]; then
    echo "## studio: SDL_VIDEODRIVER=dummy bin/pid_studio --preset=mat45r --rule=guided --grow-seed=9300 --steps=1000 --tpf=20"
    echo "   --show=disks,hull,food,particles --alpha=0.35 --shot-tick=600 --auto --grow-only --exit-at-end (font: see its terminal output)"
    echo "bin/pid_studio sha256 $(sha bin/pid_studio)"
    echo "$(sha shots/demo1c_studio_mat45r_tick600.png)  demo1c_studio_mat45r_tick600.png ($(wc -c < shots/demo1c_studio_mat45r_tick600.png | tr -d ' ') bytes)"
  fi
  if cmp -s "$OUT/network_particles/pstate_000100.txt" "$OUT/ke_overlay/pstate_000100.txt"; then echo "both shots show the same state (pstate_000100.txt equal)"
  else echo "the two runs' pstate_000100.txt differ"; fi
  echo "Rendering (SDL's software renderer, the dummy driver) is the cloud's; a window on the Mac draws the same picture by another"
  echo "renderer, so the PNG bytes are not compared across machines; the state files are (make test, test 11)."
} > "$M"
echo "shots: $(ls shots | tr '\n' ' ')"
