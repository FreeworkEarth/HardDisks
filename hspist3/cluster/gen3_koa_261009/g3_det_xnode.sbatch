#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. Cross-node determinism of gen3 (part (b) of the KOA gate): the
# same binary and the same seed on TWO DIFFERENT nodes (--nodes=2, one task per node) must give byte-identical outputs: the divider
# trace, the psi6(t) file, speed_of_sound_psi6.csv, and the run record's event hash. Two cells: N = 400 at eta 0.70 (M = 500, M4
# seeding, exact box, 2000 sigma-time held, 200 periods released) and the T-prime cell (N = 100, pi/8, M = 300, the T-prime protocol).
# Run from ~/harddisks_gen3/hspist3 after build_gen3_koa.sh: sbatch cluster/gen3_koa_261009/g3_det_xnode.sbatch
#SBATCH --job-name=gen3-xnode
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: run from ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
head=$(git rev-parse --short HEAD)
sha256sum --status -c logs/BUILD_KOA_LAST.hash && ./00ALLINONE --version | head -1 | grep -q -- "git $head  target koa" \
  || { echo "STOP: ./00ALLINONE is not the recorded clean koa build of HEAD $head"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/gen3_xnode_${SLURM_JOB_ID}"; mkdir -p "$OUT"
echo "nodes: $SLURM_JOB_NODELIST | out $OUT"; ./00ALLINONE --version | head -1
runpair() {   # runpair <cell name> <arguments...>: the same command on node 0 (A) and node 1 (B), each in its own folder
  local cell=$1; shift
  for tag in A B; do mkdir -p "$OUT/$cell/$tag"; done
  srun --nodes=1 --ntasks=1 --relative=0 --exact bash -c "cd '$OUT/$cell/A' && hostname > node.txt && HD_CONTACT_AUDIT=1 '$PWD/00ALLINONE' $* --speed-sound-run-dir='$OUT/$cell/A' > run.log 2>&1" &
  srun --nodes=1 --ntasks=1 --relative=1 --exact bash -c "cd '$OUT/$cell/B' && hostname > node.txt && HD_CONTACT_AUDIT=1 '$PWD/00ALLINONE' $* --speed-sound-run-dir='$OUT/$cell/B' > run.log 2>&1" &
  wait
  echo "== $cell: node A $(cat "$OUT/$cell/A/node.txt"), node B $(cat "$OUT/$cell/B/node.txt")"
  for f in $(cd "$OUT/$cell/A" && ls wall_x_positions_*.csv psi6_t_*.csv speed_of_sound_psi6.csv 2>/dev/null); do
    cmp -s "$OUT/$cell/A/$f" "$OUT/$cell/B/$f" && echo "   $f: IDENTICAL" || echo "   $f: **DIFFERENT**"
  done
  ha=$(grep -o "hash=[0-9a-f]*" "$OUT/$cell/A/run.log"); hb=$(grep -o "hash=[0-9a-f]*" "$OUT/$cell/B/run.log")
  [ -n "$ha" ] && [ "$ha" = "$hb" ] && echo "   event $ha on both: IDENTICAL" || echo "   event hash A '$ha' B '$hb': **DIFFERENT or missing**"
  grep -h "^\[EDMD3-HEALTH\]" "$OUT/$cell/A/run.log" "$OUT/$cell/B/run.log" | grep -o "clean=[01]" | tr '\n' ' '; echo
}
L0=$(python3 -c "import math; print('%.17g' % (400 * math.pi / (8 * 20 * 0.70)))")
runpair N400_eta070_M500 --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=400 --particles-boxes=200,200 --height=20 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
  --lengths=$L0 --wall-masses=500 --repeats=1 --seed=20261219 --wall-hold-steps=120000 --fixed-dt=0.4 --target-oscillations=200 \
  --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=60 \
  --speed-sound-exact-seed=4242 --engine=gen3 --gen3-seeding=lattice --gen3-exact-box
runpair Tprime_cell_M300 --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=100 --particles-boxes=50,50 --height=10.000000 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
  --lengths=10.000000 --wall-masses=300 --repeats=1 --seed=20261219 --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=200 \
  --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=51 \
  --speed-sound-exact-seed=4243 --engine=gen3
echo "cross-node determinism (gen3) done; paste from 'nodes:' to here"
