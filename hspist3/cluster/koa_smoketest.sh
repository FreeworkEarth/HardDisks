#!/bin/bash
# ##CHRIS 2026-10-13: KOA smoke test for the Paper 1 confinement campaign (261012 sec. 1.10). Run it BEFORE any array.
#   cd ~/harddisks/hspist3 && sbatch cluster/koa_smoketest.sh        (or: bash cluster/koa_smoketest.sh on sandbox)
#
# STEPS: (1) make -B koa   (2) ./00ALLINONE --version   (3) determinism self-test: the same seed twice, cmp must say
# identical   (4) the pi/8 pilot -- the pre-registered anchor cell (eta = pi/8, H = L_0 = 10, N_s = 50, t = 0.05),
# method B, the nine A1v2 masses x ONE seed, 200 oscillations; seeds from tests_20260913.run_seed(20261013, 0, m, 0)
# -- then compare with the Mac below. One script (cluster/confinement_pilot.py) runs and analyses on both machines.
#
# MAC TARGET (same seeds; 2026-10-13; binary "00ALLINONE git 05215ea-dirty target release",
#             CFLAGS -O3 -march=native -ffp-contract=off; "-dirty" = unrelated files, core sources = 05215ea):
#   eta (trace)     = 0.392699          L_0 (trace) = 10.0      H = 10.0      N_s = 50      r = 0.5
#   t               = 0.05 (set by --wall-thickness; speed-of-sound mode does not write it, so it is checked
#                     by the mode-equivalence gate on the Mac, 261012 sec. 1.10, not here)
#   L_eff           = L_0 - 2r - t/2 = 8.975000
#   T_total         = 69944.5 sigma-time (sum of the nine planned durations); health lines = 0
#   per-mass nu     = M50 0.07519225, M100 0.05870708, M200 0.04395559, M300 0.03664223, M500 0.02912381,
#                     M750 0.02392290, M1000 0.02092928, M1500 0.01742544, M2000 0.01506197
#   c_s             = 3.85886 +- 0.05150   (through-origin slope; +- = 1-sigma mass scatter of the implied c_s)
#
# GATES (pre-registered, 261012 sec. 1.10):
#   determinism     KOA run twice with the same seed -> cmp IDENTICAL (same binary, same node).
#   geometry        eta, L_0, H, L_eff equal to the Mac to 1e-6; t equal by the Mac mode-equivalence gate.
#   statistics      |c_s(KOA) - c_s(Mac)| <= 0.05150, the Mac pilot's 1-sigma mass scatter, fixed now.
#   NOT required    byte-identity Mac vs KOA. The Mac is arm64 (clang, Apple libm), KOA x86-64-v2 (gcc, glibc libm).
#                   Even with -ffp-contract=off on both, transcendental functions (log, cos, exp in the velocity
#                   draw) are not correctly rounded and differ in the last bit between the two libraries; one ulp
#                   in an initial velocity is amplified by the chaotic dynamics within a few hundred collisions.
#                   The two runs are therefore independent realisations, and the gate is statistical.
#
# AMENDED 2026-10-14 (261012 sec. 1.10.1; printed by cluster/smoketest_gate_width_20261014.py) -- the statistics gate:
#   0.05150 above is np.std(imp, ddof=1) (cluster/confinement_pilot.py:64): the SCATTER (SD) of the nine per-mass implied
#   c_s, NOT a standard error. The pilot c_s is the through-origin slope, a weighted mean with w = x^2 (effective n = 4.45),
#   so SE = s sqrt(sum w^2)/sum w = 0.02441 and two independent pilots differ with sigma_diff = sqrt(2) SE = 0.03452.
#   The gate is now 2 sigma_diff:  |c_s(KOA) - c_s(Mac)| <= 0.06903  (false-fail under the null 4.55 %). The old 0.05150
#   sat at 1.49 sigma_diff and would have failed a correct KOA build 13.6 % of the time.
# FILLED 2026-10-02 (Task K1/K2; KOA facts from Chris's terminal, 2026-10-03 UTC; the placeholders are gone):
#   partition sandbox (4 h; this job needs < 1 h), account uh, scratch /mnt/lustre/koa/scratch/charing, environment
#   cluster/koa_env.sh (module compiler/GCC/14.3.0 + the conda env ~/envs/hd for SDL2/SDL2_ttf/GLEW and python).
#   /home is noexec on the login node, so the build runs HERE, inside the job (cluster/build_koa.sh); the binary stays in
#   the checkout and runs on the compute node. build_koa.sh stops unless build_git == the checkout's clean HEAD.
#   Determinism: the two runs are two separate srun steps. Inside this one-node job both land on the SAME node by
#   construction; that the bytes are also identical ACROSS nodes (the point of -march=x86-64-v2) is tested by
#   cluster/koa_crossnode_det.sh, which places the two steps on two different nodes. The gates below are unchanged.
#   Submit from ~/harddisks/hspist3 after `mkdir -p logs` (Slurm does not create the log directory):
#       sbatch cluster/koa_smoketest.sh            # a re-run needs a fresh name: SMOKE_TAG=koa_pi8_try2 sbatch ...
#SBATCH --job-name=conf-smoke
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks/hspist3}"
SCRATCH="/mnt/lustre/koa/scratch/charing"
REL=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013
# ##CHRIS 2026-10-05 (engine gate, 261012 sec. 4.4): the data root follows the clone this job was submitted from
# ($SCRATCH/harddisks for ~/harddisks, $SCRATCH/harddisks_resched for the engine-branch clone), so a second build
# generation never writes into the 279282b tree. For ~/harddisks this is the old path, unchanged.
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
OUT="$SCRATCH/$ROOT_NAME/hspist3/$REL/${SMOKE_TAG:-koa_pi8_H10_L10}"
[ -e "$OUT" ] && { echo "$OUT exists -- not overwriting; choose a fresh name with SMOKE_TAG=..."; exit 2; }

echo "== (0) node";    echo "$(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | job $SLURM_JOB_ID"
source cluster/koa_env.sh || exit 1
# ##CHRIS 2026-10-02: preflight -- the analysis modules must import BEFORE anything is built or run (seconds, not minutes)
echo "== (0b) analysis modules"
python3 -c "import sys; sys.path[:0] = ['validation', '.']; import tests_20260913, plot_speed_of_sound_edmd, paper1_populate_cs_err_20261002" \
  || { echo "STOP: an analysis module does not import -- plot_speed_of_sound_edmd.py must be tracked, pushed and pulled (runsheet step 1)"; exit 1; }
echo "== (1) build";   bash cluster/build_koa.sh || exit 1
echo "== (2) version"; gcc --version | head -1; ./00ALLINONE --version | head -2
./00ALLINONE --version | grep -q -- "-ffp-contract=off" || { echo "STOP: binary lacks -ffp-contract=off"; exit 1; }
echo "== (3) determinism (two srun steps, same node)"
for tag in A B; do
  srun --ntasks=1 --cpus-per-task=1 --exact python3 cluster/confinement_pilot.py det1 --tag $tag --bin ./00ALLINONE --out "$OUT/_determinism" || exit 1
done
python3 cluster/confinement_pilot.py detcmp --out "$OUT/_determinism" || exit 1
echo "== (4) pi/8 pilot"; python3 cluster/confinement_pilot.py run --bin ./00ALLINONE --out "$OUT" --jobs "${SLURM_CPUS_PER_TASK:-9}" || exit 1
python3 cluster/confinement_pilot.py analyse --out "$OUT" | tee "$OUT/pilot_analysis.txt"
echo "== gates"
python3 - "$OUT/pilot_analysis.txt" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
eta = float(re.search(r"eta \(trace\) = \[([\d.]+)\]", t).group(1)); L0 = float(re.search(r"L_0 \(trace\) = \[([\d.]+)\]", t).group(1))
Le = float(re.search(r"L_eff = L_0 - 2r - t/2 = ([\d.]+)", t).group(1)); cs = float(re.search(r"c_s = ([\d.]+) \+-", t).group(1))
hl = int(re.search(r"health lines = (\d+)", t).group(1))
g = [("eta", abs(eta - 0.392699) <= 1e-6), ("L_0", abs(L0 - 10.0) <= 1e-6), ("L_eff", abs(Le - 8.975) <= 1e-6),
     ("health", hl == 0), ("c_s within 0.06903 (2 sigma_diff) of 3.85886", abs(cs - 3.85886) <= 0.06903)]
for n, ok in g: print(f"  {n}: {'PASS' if ok else 'FAIL'}")
print(f"  c_s(KOA) - c_s(Mac) = {cs - 3.85886:+.5f}")
print("SMOKE TEST", "PASSED -- the method-A pilot may be submitted (after the go)" if all(ok for _, ok in g) else "FAILED -- STOP")
PY
