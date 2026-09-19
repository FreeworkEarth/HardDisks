#!/usr/bin/env bash
# ##CHRIS 2026-09-19: Paper 2 -- run the whole ladder by hand, one step at a time.
#
# This is the script to TEST the pipeline, not the script that produced the published numbers.
# It runs the same commands with a small seed count (SEEDS, default 5) so every step finishes in
# seconds to a couple of minutes. The campaign scripts next to this file are the real thing and
# carry the published seed counts; this one exists so the apparatus can be exercised end to end.
#
#   ./paper2_run_all.sh                 list the steps and stop
#   ./paper2_run_all.sh 3               run step 3 only
#   ./paper2_run_all.sh 2 3 4           run those steps in order
#   ./paper2_run_all.sh all             run every step, pausing before each one
#   ./paper2_run_all.sh all --yes       run every step without pausing
#   SEEDS=25 ./paper2_run_all.sh all    heavier, closer to the campaign
#   OUT=/tmp/p2test ./paper2_run_all.sh all      write somewhere else
#
# Every step prints a PASS/FAIL line. FAIL is decided by the health contract, not by eye:
# a step fails if the binary aborted, if any EDMD health event fired, or if the step's own
# numeric criterion is missed. The script keeps going after a FAIL so you see the whole picture.
#
# ---------------------------------------------------------------------------------------------
# THE MASTER BOX -- one apparatus, four geometries. All of Paper 2 runs in this box.
#
#   |<-- 10 sigma empty -->|  wall_S  |<-- gas 1, 38.75 x 10 -->| divider |<-- gas 2, 38.75 x 10 -->| piston
#   0                    10.0      11.0                      49.75     50.75                     89.50
#
#   --l0=44.75 (the binary doubles it: box 0 .. 89.50), --height=10, wall thickness 1.0 (default),
#   wall_S centre 10.5, divider centre 50.25, 50 disks per gas -> eta = 0.1013.
#   The spring lives in the left compartment, which is ALWAYS empty: no gas ever touches it.
#
#   A  wall_S held, divider held           <- what Levels 0-2 measured
#   B  wall_S held, divider free (M = 1000)
#   C  wall_S free on the spring, divider REMOVED: one gas of 100 disks over 78.5 sigma
#   D  both free                           <- the chain
#
#   Compression is 10 % of the gas length in every geometry: dx = 3.93 sigma for a 38.75 sigma
#   compartment (eta 0.101 -> 0.113), dx = 7.96 sigma for geometry C's 78.5 sigma gas.
#
# GRID TRAP: every wall position and the box length must be an integer number of 1/24 sigma, or
# the binary aborts with [initial_wall_position_mismatch]. 10.5, 50.25, 44.75 all are. The piston
# travel is NOT subject to this -- 3.93 is fine.
# ---------------------------------------------------------------------------------------------
set -uo pipefail

HS=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
# Geometries C and D need --spring-wall, which the installed binary only has once the spring-hygiene
# build is promoted. Until then fall back to the staged build and say so, loudly, every run.
BIN=${BIN:-}
if [ -z "$BIN" ]; then
  BIN=$HS/00ALLINONE
  if ! "$BIN" --help 2>&1 | grep -q -- '--spring-wall' && \
     [ -x "$HS/00ALLINONE_sp" ] && "$HS/00ALLINONE_sp" --help 2>&1 | grep -q -- '--spring-wall'; then
    BIN=$HS/00ALLINONE_sp
    echo "NOTE: 00ALLINONE has no --spring-wall yet, using the staged build 00ALLINONE_sp."
    echo "      Promote it (build + 3-seed byte-identity gate) and this note goes away."
  fi
fi
OUT=${OUT:-$HS/experiments_energy_transfer/paper2_pipeline_test}
SEEDS=${SEEDS:-5}
JOBS=${JOBS:-8}
PAUSE=1
FAILED=()

cd "$HS" || exit 1

# ---- the master box, as shell arrays so no step can drift from another -----------------------
BOX=(--particle-radius=0.5 --l0=44.75 --height=10)
GAS2=(--particles=100 --particles-boxes=0,50,50 --num-walls=2 --wall-positions=10.5,50.25)
GAS1=(--particles=100 --particles-boxes=0,100  --num-walls=1 --wall-positions=10.5)
HELD=1000000000
RUN=(--mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0
     --seed-drift-order=drift-first --fixed-dt=0.4 --kbt1 --energy-measurement)
DX=3.93          # 10 % of a 38.75 sigma compartment
DXC=7.96         # the same 10.14 % of geometry C's 78.5 sigma gas

c()  { printf '\033[%sm%s\033[0m' "$1" "$2"; }
hdr(){ echo; echo "$(c '1;36' "=== step $1 -- $2")"; echo "    $3"; }
pass(){ echo "    $(c '1;32' PASS)  $1"; }
fail(){ echo "    $(c '1;31' FAIL)  $1"; FAILED+=("step $CUR: $1"); }
ask(){ [ "$PAUSE" = 1 ] || return 0; read -r -p "    press ENTER to run, s to skip: " a; [ "$a" != s ]; }

# health contract: nothing may abort, no health event, and the four counters stay zero.
health(){                     # $1 = log file, $2 = expected number of traces, $3 = trace glob
  local log=$1 want=$2 glob=$3 ok=1 ab he got
  # grep -c already prints 0 when it finds nothing; do NOT add '|| echo 0' -- grep exits 1 on
  # no-match and the substitution then contains "0\n0", which reads as a failure that is not one.
  ab=$(grep -c ABORTING "$log" 2>/dev/null); ab=${ab:-0}
  he=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$log" 2>/dev/null); he=${he:-0}
  got=$(ls $glob 2>/dev/null | wc -l | tr -d ' ')
  grep -q "Unknown option" "$log" 2>/dev/null && \
    fail "the binary rejected a flag: $(grep -m1 "Unknown option" "$log")"
  [ "$ab" = 0 ] || { fail "$ab runs aborted -- see $log"; ok=0; }
  [ "$he" = 0 ] || { fail "$he health events -- see $log"; ok=0; }
  [ "$got" = "$want" ] || { fail "$got/$want trajectories produced"; ok=0; }
  [ $ok = 1 ] && pass "$got/$want runs, 0 aborts, 0 health events"
  return $((1 - ok))
}

# one cell of N seeds at speed $u; extra flags after the first four arguments.
cell(){
  local dir=$1 u=$2 steps=$3 n=$4; shift 4
  mkdir -p "$dir"
  local k sd
  for k in $(seq 0 $((n - 1))); do
    sd=$((9200 + k)); [ -s "$dir/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$dir/ev_${sd}.csv" nice -n 5 "$BIN" "${RUN[@]}" "${BOX[@]}" "$@" \
        --energy-transfer-summary="$dir/summary.csv" --energy-transfer-trace="$dir/tr_${sd}.csv" \
        --piston-right-protocol-mode=step --velocity-right-piston-step=$u --auto-piston-step \
        --wall-hold-steps=12000 --steps=$steps --seed=$sd >> "$dir/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 0.3; done
  done
  wait
}

# mean +- standard error of one summary column, printed as "mean sem n"
col(){ python3 - "$1" "$2" <<'PY'
import sys, math, pandas as pd, numpy as np
d = pd.read_csv(sys.argv[1], low_memory=False)
d = d[d['timestamp'] != 'timestamp']                 # summary.csv re-emits its header per run
a = d[sys.argv[2]].astype(float).to_numpy()
print(f"{a.mean():.6f} {a.std(ddof=1)/math.sqrt(len(a)):.6f} {len(a)}")
PY
}

# =============================================================================================
step1(){ hdr 1 "the binary and the box" \
  "checks the binary exists, is the scientific backend, and that the master box starts up clean"
  ask || return 0
  [ -x "$BIN" ] || { fail "no binary at $BIN -- run 'make' in $HS"; return; }
  echo "    binary   $BIN  ($(shasum -a 256 "$BIN" | cut -c1-12), $(date -r "$BIN" '+%Y-%m-%d %H:%M'))"
  local d=$OUT/step1; rm -rf "$d"; mkdir -p "$d"
  cell "$d" 0.05 17000 1 "${GAS2[@]}" --wall-mass-factors=$HELD,$HELD --eff-output=wall-ke \
       --max-right-piston-travel=$DX
  # the binary prints the geometry it actually built; that is what we check, not what we asked for
  grep -iE 'L_eff|compartment|wall .* at|thickness' "$d/run.log" | head -6 | sed 's/^/    | /'
  health "$d/run.log" 1 "$d/tr_*.csv"
}

step2(){ hdr 2 "geometries A B C D start up" \
  "one seed in each of the four geometries -- catches a bad wall position or a gas in the spring compartment"
  ask || return 0
  local g ok=1
  for g in A B C D; do
    local d=$OUT/step2/$g; rm -rf "$d"
    case $g in
      A) cell "$d" 0.05 17000 1 "${GAS2[@]}" --wall-mass-factors=$HELD,$HELD --eff-output=wall-ke --max-right-piston-travel=$DX ;;
      B) cell "$d" 0.05 17000 1 "${GAS2[@]}" --wall-mass-factors=$HELD,1000 --eff-output=wall-ke --max-right-piston-travel=$DX ;;
      C) cell "$d" 0.05 20000 1 "${GAS1[@]}" --wall-mass-factors=200 --spring-k=5 --spring-wall=0 \
              --spring-eq=10.5 --eff-output=spring --max-right-piston-travel=$DXC ;;
      D) cell "$d" 0.05 20000 1 "${GAS2[@]}" --wall-mass-factors=200,1000 --spring-k=5 --spring-wall=0 \
              --spring-eq=10.5 --eff-output=spring --max-right-piston-travel=$DX ;;
    esac
    CUR="2$g"; printf "  %s: " "$g"; health "$d/run.log" 1 "$d/tr_*.csv" || ok=0
  done
  [ $ok = 1 ] && pass "all four geometries built and ran"
}

step3(){ hdr 3 "Level 0 -- the ledger closes" \
  "$SEEDS seeds in geometry A; the piston work must equal the gas kinetic-energy change to 1e-4"
  ask || return 0
  local d=$OUT/step3; mkdir -p "$d"
  cell "$d" 0.05 17000 "$SEEDS" "${GAS2[@]}" --wall-mass-factors=$HELD,$HELD --eff-output=wall-ke \
       --max-right-piston-travel=$DX
  health "$d/run.log" "$SEEDS" "$d/tr_*.csv" || return
  python3 - "$d" <<'PY'
import glob, sys, pandas as pd
worst = 0.0
for tr in sorted(glob.glob(sys.argv[1] + "/tr_*.csv")):
    t = pd.read_csv(tr, low_memory=False)
    W = float(t["PistonWork"].iloc[-1])
    g = float(t["KE_gas_total"].iloc[-1] - t["KE_gas_total"].iloc[0])
    worst = max(worst, abs(W - g) / abs(W))
print(f"    worst |W - dKE_gas| / W = {worst:.3e}   "
      + ("PASS (<= 1e-4)" if worst <= 1e-4 else "FAIL (> 1e-4)"))
sys.exit(0 if worst <= 1e-4 else 1)
PY
  [ $? = 0 ] && pass "ledger closes" || fail "ledger residual above 1e-4"
}

step4(){ hdr 4 "Level 1 -- the quasi-static baseline" \
  "extrapolate <W> to u -> 0 over three slow speeds and compare with W_qs = 7.4685 +- 0.0131 kT"
  ask || return 0
  local u
  for u in 0.005 0.02 0.05; do
    local st=17000; [ "$u" = 0.005 ] && st=60000; [ "$u" = 0.02 ] && st=30000
    cell "$OUT/step4/u$u" "$u" "$st" "$SEEDS" "${GAS2[@]}" --wall-mass-factors=$HELD,$HELD \
         --eff-output=wall-ke --max-right-piston-travel=$DX
  done
  health "$(cat $OUT/step4/u*/run.log > $OUT/step4/all.log; echo $OUT/step4/all.log)" \
         $((3 * SEEDS)) "$OUT/step4/u*/tr_*.csv" || return
  python3 - "$OUT/step4" <<'PY'
import glob, math, sys, numpy as np, pandas as pd
WQS, dWQS = 7.4685, 0.0131                      # 260919 report, finite-box path integral
rows = []
for d in sorted(glob.glob(sys.argv[1] + "/u*")):
    u = float(d.rsplit("u", 1)[1])
    s = pd.read_csv(d + "/summary.csv", low_memory=False)
    s = s[s["timestamp"] != "timestamp"]
    w = s["W_in_max"].astype(float).to_numpy()
    rows.append((u, w.mean(), w.std(ddof=1) / math.sqrt(len(w))))
    print(f"    u = {u:<6g} <W_in> = {rows[-1][1]:.4f} +- {rows[-1][2]:.4f} kT   ({len(w)} seeds)")
u = np.array([r[0] for r in rows]); w = np.array([r[1] for r in rows]); e = np.array([r[2] for r in rows])
# the excess is quadratic in u (Level 2), so extrapolate on u^2
A = np.vstack([np.ones_like(u), u ** 2]).T
C = np.diag(e ** 2)
iC = np.linalg.inv(C)
cov = np.linalg.inv(A.T @ iC @ A)
p = cov @ (A.T @ iC @ w)
w0, dw0 = p[0], math.sqrt(cov[0, 0])
gap = w0 - WQS; dgap = math.hypot(dw0, dWQS)
n = abs(gap) / dgap
print(f"    W(u -> 0) = {w0:.4f} +- {dw0:.4f}   against W_qs = {WQS} +- {dWQS}")
print(f"    gap = {gap:+.4f} +- {dgap:.4f} kT  ->  {n:.2f} sigma  " + ("PASS (<= 2)" if n <= 2 else "FAIL (> 2)"))
sys.exit(0 if n <= 2 else 1)
PY
  [ $? = 0 ] && pass "quasi-static baseline recovered" || fail "extrapolated W disagrees with W_qs"
}

step5(){ hdr 5 "Level 2 -- the finite-rate excess is quadratic" \
  "five speeds; <W> - W_qs must grow as A u^2 with A of order 26 for a stepped piston.
    At a handful of seeds the per-speed A values scatter hard, some negative -- that is the
    thermal spread of W, not a bug. Only the fitted A over all five speeds is the check."
  ask || return 0
  local u
  for u in 0.02 0.05 0.10 0.15 0.20; do
    local st=17000; [ "$u" = 0.02 ] && st=30000
    cell "$OUT/step5/u$u" "$u" "$st" "$SEEDS" "${GAS2[@]}" --wall-mass-factors=$HELD,$HELD \
         --eff-output=wall-ke --max-right-piston-travel=$DX
  done
  health "$(cat $OUT/step5/u*/run.log > $OUT/step5/all.log; echo $OUT/step5/all.log)" \
         $((5 * SEEDS)) "$OUT/step5/u*/tr_*.csv" || return
  python3 - "$OUT/step5" <<'PY'
import glob, math, sys, numpy as np, pandas as pd
WQS = 7.4685
u, w, e = [], [], []
for d in sorted(glob.glob(sys.argv[1] + "/u*")):
    s = pd.read_csv(d + "/summary.csv", low_memory=False)
    s = s[s["timestamp"] != "timestamp"]
    a = s["W_in_max"].astype(float).to_numpy()
    u.append(float(d.rsplit("u", 1)[1])); w.append(a.mean()); e.append(a.std(ddof=1) / math.sqrt(len(a)))
u, w, e = map(np.array, (u, w, e))
for i in range(len(u)):
    print(f"    u = {u[i]:<6g} <W> - W_qs = {w[i]-WQS:+.4f} +- {e[i]:.4f}   A = {(w[i]-WQS)/u[i]**2:7.2f}")
# one-parameter fit through the origin in u^2
A = float(np.sum((w - WQS) * u**2 / e**2) / np.sum(u**4 / e**2))
dA = float(1.0 / math.sqrt(np.sum(u**4 / e**2)))
print(f"    A_step = {A:.2f} +- {dA:.2f}   (campaign value 25.9 +- 4.3)")
n = abs(A - 25.9) / math.hypot(dA, 4.3)
print("    " + (f"PASS ({n:.2f} sigma from the campaign value)" if n <= 3 else f"FAIL ({n:.2f} sigma)"))
sys.exit(0 if n <= 3 else 1)
PY
  [ $? = 0 ] && pass "quadratic excess reproduced" || fail "A_step disagrees with the campaign"
}

step6(){ hdr 6 "Level 3 -- the spring captures some of it" \
  "geometry C at two speeds; the spring must gain energy and the ledger must still close"
  ask || return 0
  local u
  for u in 0.02 0.20; do
    local st=20000; [ "$u" = 0.02 ] && st=40000
    cell "$OUT/step6/u$u" "$u" "$st" "$SEEDS" "${GAS1[@]}" --wall-mass-factors=200 --spring-k=5 \
         --spring-wall=0 --spring-eq=10.5 --eff-output=spring --max-right-piston-travel=$DXC
  done
  health "$(cat $OUT/step6/u*/run.log > $OUT/step6/all.log; echo $OUT/step6/all.log)" \
         $((2 * SEEDS)) "$OUT/step6/u*/tr_*.csv" || return
  python3 - "$OUT/step6" <<'PY'
import glob, math, sys, numpy as np, pandas as pd
ok = True
for d in sorted(glob.glob(sys.argv[1] + "/u*")):
    u = float(d.rsplit("u", 1)[1]); eps = []
    for tr in sorted(glob.glob(d + "/tr_*.csv")):
        t = pd.read_csv(tr, low_memory=False)
        v = np.abs(t["PistonR_v"].to_numpy(float)); mv = np.nonzero(v > 1e-12)[0]
        i0, i1 = mv[0], mv[-1]
        W = float(t["PistonWork"].iloc[i1] - t["PistonWork"].iloc[i0])
        E = t["SpringE"].to_numpy(float)
        base = float(np.mean(E[max(0, i0 - 200):i0 + 1]))
        eps.append(float(np.max(E[i1:] - base)) / W)
    m = float(np.mean(eps)); s = float(np.std(eps, ddof=1) / math.sqrt(len(eps)))
    good = m > 0
    ok &= good
    print(f"    u = {u:<6g} eps = E_spring,max / W_in = {m:.5f} +- {s:.5f}  "
          + ("(positive)" if good else "(NOT POSITIVE -- the spring is not being driven)"))
print("    NOTE the Level 3 pilot found this observable thermal-noise dominated: the 1/2 kT pedestal")
print("         is far above the quasi-static capture E_qs. A positive eps here is a plumbing check,")
print("         not a physics result. See 260920_paper2_level3_REPORT.md.")
sys.exit(0 if ok else 1)
PY
  [ $? = 0 ] && pass "spring is driven and the trace carries SpringE" || fail "spring did not gain energy"
}

step7(){ hdr 7 "the pictures (needs a screen)" \
  "opens the GUI in each geometry and saves one screenshot per render mode -- skip this over ssh"
  ask || return 0
  local sh=$HS/experiments_energy_transfer/_run_scripts/paper2_geometry_pictures.sh
  if [ -x "$sh" ]; then bash "$sh" "$OUT/step7" && pass "8 pictures in $OUT/step7" || fail "picture run failed"
  else echo "    $(c '1;33' SKIP)  $sh not present"; fi
}

step8(){ hdr 8 "the analysis scripts" \
  "re-runs the published analyses against the published run directories -- reads only, changes nothing"
  ask || return 0
  local s ok=1
  for s in paper2_level2_20260917.py paper2_geometry_fix_20260918.py paper2_level3_20260919.py; do
    printf "    %-34s " "$s"
    if python3 "$HS/validation/$s" > "$OUT/${s%.py}.out" 2>&1; then echo "$(c '1;32' ok)   -> $OUT/${s%.py}.out"
    else echo "$(c '1;31' failed)  -> $OUT/${s%.py}.out"; ok=0; fi
  done
  [ $ok = 1 ] && pass "all analyses ran" || fail "an analysis script failed"
}

# =============================================================================================
STEPS=(step1 step2 step3 step4 step5 step6 step7 step8)
NAMES=("the binary and the box" "geometries A B C D start up" "Level 0 -- the ledger closes" \
       "Level 1 -- the quasi-static baseline" "Level 2 -- the finite-rate excess" \
       "Level 3 -- the spring" "the pictures (needs a screen)" "the analysis scripts")

usage(){
  echo "Paper 2 -- run the ladder step by step.   output: $OUT   seeds per cell: $SEEDS"
  echo
  for i in "${!STEPS[@]}"; do printf "  %d  %s\n" $((i + 1)) "${NAMES[$i]}"; done
  echo
  echo "  ./paper2_run_all.sh 3          one step"
  echo "  ./paper2_run_all.sh all        all of them, pausing before each"
  echo "  ./paper2_run_all.sh all --yes  all of them, no pauses"
  echo "  SEEDS=25 JOBS=10 ./paper2_run_all.sh all       heavier"
}

ARGS=()
for a in "$@"; do case $a in --yes|-y) PAUSE=0 ;; *) ARGS+=("$a") ;; esac; done
[ ${#ARGS[@]} -eq 0 ] && { usage; exit 0; }
[ "${ARGS[0]}" = all ] && ARGS=(1 2 3 4 5 6 7 8)

mkdir -p "$OUT"
echo "output -> $OUT    seeds per cell: $SEEDS    parallel jobs: $JOBS"
for n in "${ARGS[@]}"; do
  CUR=$n
  case $n in [1-8]) "${STEPS[$((n - 1))]}" ;; *) echo "no step '$n'"; usage; exit 1 ;; esac
done

echo
if [ ${#FAILED[@]} -eq 0 ]; then echo "$(c '1;32' 'ALL STEPS PASSED')"
else echo "$(c '1;31' "${#FAILED[@]} FAILURES")"; printf '  %s\n' "${FAILED[@]}"; fi
echo "everything written under $OUT"
