#!/bin/bash
# ##CHRIS 2026-09-16: bring KOA results back, into their own campaign root so cluster and laptop
# trajectories can never be mixed inside one cell (see README.md).
#
# usage: bash cluster/fetch_results.sh <ssh-host> <remote-runs-dir> [campaign]
set -euo pipefail
HOST="${1:?usage: fetch_results.sh HOST REMOTE_RUNS_DIR [campaign]}"
REMOTE="${2:?remote runs dir, e.g. ~/harddisks/runs}"
CAMPAIGN="${3:-}"
cd "$(dirname "$0")/.."
LOCAL="experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/KOA_$(date +%Y%m%d)"
mkdir -p "$LOCAL"

rsync -avh --partial --info=progress2 \
  --include '*/' --include 'wall_x_positions_*.csv' --include 'stdout.log' \
  --include 'command.txt' --include 'speed_of_sound_psi6.csv' --include 'run_records_*.tsv' \
  --exclude '*' \
  "$HOST:$REMOTE/${CAMPAIGN:+$CAMPAIGN/}" "$LOCAL/${CAMPAIGN:+$CAMPAIGN/}"

echo
echo "fetched into $LOCAL"
echo "traces:        $(find "$LOCAL" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')"
echo "runs recorded: $(cat "$LOCAL"/run_records_*.tsv 2>/dev/null | wc -l | tr -d ' ')"
echo "failed runs:   $(awk -F'\t' '$5 != 1' "$LOCAL"/run_records_*.tsv 2>/dev/null | wc -l | tr -d ' ')"
echo "health lines:  $(awk -F'\t' '{s += $6} END {print s + 0}' "$LOCAL"/run_records_*.tsv 2>/dev/null)"
echo
echo "Remember: these are a separate campaign. Do not merge them into famB/A1v2 cells."
