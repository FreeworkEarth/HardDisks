#!/bin/bash
# ##CHRIS 2026-10-13 (plan-author decision 48, item 3): one test scene of the viewer, in the window's scene loop but with the SDL dummy
# video driver (no window on the screen), from its own new scratch folder (the scene writes its CSVs into the working directory; the
# font is linked in), with a hard kill after LIMIT seconds. stdout and stderr go to the evidence folder.
# usage: run_scene.sh NAME BINARY LIMIT_S -- <the viewer's flags>
set -uo pipefail
NAME=$1; BIN=$2; LIMIT=$3; shift 4
G=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
RUN="$SP/d48/runs/$NAME"; EV="$G/experiments_gen3_d48_261013/tests"
[ -e "$RUN" ] && { echo "REFUSED: $RUN exists (a new name, nothing is deleted)"; exit 1; }
[ -e "$EV/$NAME.out" ] && { echo "REFUSED: $EV/$NAME.out exists"; exit 1; }
mkdir -p "$RUN" "$EV"
ln -s "$G/Roboto-Regular.ttf" "$RUN/Roboto-Regular.ttf"
{ echo "# $NAME: $(date '+%Y-%m-%d %H:%M:%S %Z')"; echo "# binary: $BIN ($(shasum -a 256 "$BIN" | cut -c1-16))"; printf '# flags:'; printf ' %q' "$@"; echo; } > "$EV/$NAME.cmd"
cd "$RUN"
t0=$(date +%s)
SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy perl -e 'alarm shift; exec @ARGV' "$LIMIT" "$BIN" "$@" > "$EV/$NAME.out" 2> "$EV/$NAME.err"
rc=$?
echo "# exit code $rc after $(( $(date +%s) - t0 )) s" >> "$EV/$NAME.cmd"
echo "$NAME: exit $rc in $(( $(date +%s) - t0 )) s"
grep -h "STOPPED\|STOP (\|\[EDMD3\] built\|\[GUI-KEYS\]\|\[GUI-STATS\]\|APPLY refused\|Failed" "$EV/$NAME.out" "$EV/$NAME.err" | cut -c1-260
grep -h "\[EDMD3-HEALTH\]" "$EV/$NAME.out" | sed -E 's/ (overlap_repair|contact_now|local_checks|ev_pair|validator_every)/\n   \1/g' | cut -c1-260
