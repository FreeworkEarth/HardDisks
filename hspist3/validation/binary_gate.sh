#!/usr/bin/env bash
# ##CHRIS 2026-09-20: the byte-identity gate every binary change must pass before promotion.
#
# Runs the same three seeds through two binaries on two experiments and compares the outputs
# BYTE FOR BYTE -- traces and per-event piston logs both, because a change can leave the sampled
# trace alone and still move the event stream. Pass/fail, no tolerance, no "close enough".
#
#   ./binary_gate.sh OLD NEW [workdir]
#   ./binary_gate.sh ./00ALLINONE ./00ALLINONE_sp
#
# Exit 0 = every file identical = safe to promote. Anything else and the promotion stops.
#
# Note the gate compares the two binaries on commands that DO NOT use the new flags: the point is
# that existing physics is untouched. A new flag that changes behaviour when you pass it is fine;
# a binary that changes behaviour when you don't is not.
set -uo pipefail
HS=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
OLD=${1:?usage: binary_gate.sh OLD NEW [workdir]}
NEW=${2:?usage: binary_gate.sh OLD NEW [workdir]}
W=${3:-$(mktemp -d)}
cd "$HS" || exit 1
mkdir -p "$W"

echo "gate: $(basename "$OLD") [$(shasum -a 1 "$OLD" | cut -c1-8)]  vs  $(basename "$NEW") [$(shasum -a 1 "$NEW" | cut -c1-8)]"
echo "work: $W"

BAD=0
run_pair(){                       # $1 = label, then the flags, with @DIR@/@SEED@ placeholders
  local label=$1; shift
  local sd side d f
  for side in old new; do
    local bin=$OLD; [ "$side" = new ] && bin=$NEW
    for sd in 4101 4102 4103; do
      d=$W/$label/$side/$sd; mkdir -p "$d"
      local args=()
      for f in "$@"; do args+=("${f//@DIR@/$d}"); done
      args+=("--seed=$sd")
      HD_PISTON_EVENTS="$d/events.csv" "$bin" "${args[@]}" > "$d/run.log" 2>&1
    done
  done
  # DATA files must be byte-identical. run.log and the 00_COMMAND / 01_PLOT provenance files are
  # excluded by construction: they record the binary's own path and a wall-clock timestamp, so they
  # can never match and say nothing about physics. summary.csv is excluded here and checked
  # column-wise below, because a deliberate schema addition is not a physics change.
  local n=0 bad=0
  for sd in 4101 4102 4103; do
    for f in $(cd "$W/$label/old/$sd" && ls | grep -vE '^(run\.log|summary\.csv|0[01]_.*\.md)$'); do
      n=$((n + 1))
      if ! cmp -s "$W/$label/old/$sd/$f" "$W/$label/new/$sd/$f"; then
        echo "    DIFFER  seed $sd  $f"; bad=$((bad + 1))
      fi
    done
  done
  if [ "$n" = 0 ]; then echo "  $label: NO DATA FILES PRODUCED -- gate cannot pass"; BAD=$((BAD + 1))
  elif [ "$bad" = 0 ]; then echo "  $label: $n/$n data files byte-identical across 3 seeds"
  else echo "  $label: $bad of $n data files DIFFER"; BAD=$((BAD + 1)); fi

  # summary.csv column by column. A column that changed and is NOT named in EXPECT_DIFF fails.
  if [ -f "$W/$label/old/4101/summary.csv" ]; then
    python3 - "$W/$label" "${EXPECT_DIFF:-}" <<'PY' || BAD=$((BAD + 1))
import csv, sys
base, expect = sys.argv[1], set(filter(None, sys.argv[2].split(",")))
# these three record when and by what the run was made, never what it did
IGNORE = {"timestamp", "command", "trace_path"}
bad = []
for sd in (4101, 4102, 4103):
    a = list(csv.reader(open(f"{base}/old/{sd}/summary.csv")))
    b = list(csv.reader(open(f"{base}/new/{sd}/summary.csv")))
    da, db = dict(zip(a[0], a[1])), dict(zip(b[0], b[1]))
    added = [c for c in b[0] if c not in a[0]]
    removed = [c for c in a[0] if c not in b[0]]
    if sd == 4101:
        if added:   print(f"    summary.csv: columns ADDED   {added}")
        if removed: print(f"    summary.csv: columns REMOVED {removed}")
    for k in a[0]:
        if k in IGNORE or k not in db or da[k] == db[k]:
            continue
        if k in expect:
            if sd == 4101: print(f"    summary.csv: {k} {da[k]!r} -> {db[k]!r}  (declared in EXPECT_DIFF)")
        else:
            bad.append((sd, k, da[k], db[k]))
for sd, k, x, y in bad:
    print(f"    UNDECLARED summary change  seed {sd}  {k}: {x!r} -> {y!r}")
print("    summary.csv: every other column identical" if not bad else "    summary.csv: FAIL")
sys.exit(1 if bad else 0)
PY
  fi
}

# 1. energy transfer, geometry A in the master box -- Paper 2's workhorse.
run_pair energy_transfer --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
  --seed-drift-order=drift-first --energy-transfer-summary=@DIR@/summary.csv \
  --energy-transfer-trace=@DIR@/trace.csv --particles=100 --particles-boxes=0,50,50 \
  --particle-radius=0.5 --l0=44.75 --height=10 --num-walls=2 --wall-positions=10.5,50.25 \
  --wall-mass-factors=1000000000,1000000000 --piston-right-protocol-mode=step \
  --velocity-right-piston-step=0.05 --max-right-piston-travel=3.93 --auto-piston-step \
  --wall-hold-steps=3000 --steps=6000 --fixed-dt=0.4 --energy-measurement --eff-output=wall-ke --kbt1

# 2. speed of sound -- Paper 1's workhorse, a completely different driver through the same core.
run_pair speed_of_sound --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
  --seed-drift-order=drift-first --edmd-acc=0 --particles=100 --particles-boxes=50,50 \
  --height=10 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
  --lengths=7.5 --wall-masses=200 --repeats=1 --wall-hold-steps=4000 --fixed-dt=0.4 \
  --target-oscillations=20 --oscillation-safety=1.0 --oscillation-min-steps=4000 \
  --oscillation-max-steps=40000 --speed-sound-log-stride=4 --speed-sound-run-dir=@DIR@

echo
if [ "$BAD" = 0 ]; then echo "GATE PASS -- the two binaries are byte-identical on both experiments"; exit 0
else echo "GATE FAIL -- $BAD experiment(s) differ; do NOT promote"; exit 1; fi
