#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G; corrected by plan-author decision 13, item 2, sec. 4.7.29): PREPARED, NOT RUN. Build
# the gen3 binaries on KOA from the clean clone ~/harddisks_gen3 (engine-gen3 at the commit named in the runsheet): the default
# `make koa` (-> ./00ALLINONE, recorded in logs/BUILD_KOA_LAST.hash as cluster/build_koa.sh does), the long-double option `make koa-ld`
# (-> bin/00ALLINONE_ld, logs/BUILD_KOA_LD.hash), and the M1/M2 harnesses with the koa target's flags (-> bin/gen3_m1, bin/gen3_m2,
# logs/BUILD_KOA_HARNESS.hash; -std=gnu11 because edmd.c uses M_PI, which strict -std=c11 does not define with glibc). Both drivers must
# carry "git <HEAD>" with no -dirty.
# Runs ONLY inside a Slurm job (home is noexec on the login node): from hspist3/ of ~/harddisks_gen3, at the compute-node prompt of
#   srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
#   bash cluster/gen3_koa_261009/build_gen3_koa.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 2
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"; exit 2; }
ROOT_NAME=$(basename "$(dirname "$PWD")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: build from the gen3 clone ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "STOP: the clone has modified tracked files"; git status --short | head; exit 2; }
[ "$(git -C kissfft rev-parse HEAD 2>/dev/null)" = febd4caeed32e33ad8b2e0bb5ea77542c40f18ec ] || { echo "STOP: kissfft is not the pinned commit febd4ca"; exit 2; }
source cluster/koa_env.sh || exit 2
head=$(git rev-parse --short HEAD); mkdir -p logs bin
echo "== HEAD $head $(git rev-parse HEAD) ($(git log -1 --format=%s | cut -c1-80))"
echo "== long-double build (koa-ld)"
make -B koa-ld || { echo "STOP: make koa-ld failed"; exit 4; }
v=$(./00ALLINONE --version | head -1); echo "$v"
echo "$v" | grep -q -- "git $head  target koa-ld" || { echo "STOP: the long-double build is not the clean HEAD: $v"; exit 6; }
mv ./00ALLINONE bin/00ALLINONE_ld && (cd bin && sha256sum 00ALLINONE_ld) > logs/BUILD_KOA_LD.hash && cat logs/BUILD_KOA_LD.hash
echo "== default build (koa)"
make -B koa || { echo "STOP: make koa failed"; exit 4; }
v=$(./00ALLINONE --version | head -1); echo "$v"
echo "$v" | grep -q -- "git $head  target koa" || { echo "STOP: the default build is not the clean HEAD: $v"; exit 6; }
missing=$(ldd ./00ALLINONE | grep -c "not found") || true
[ "${missing:-0}" -eq 0 ] || { echo "STOP: unresolved shared libraries"; exit 5; }
sha256sum 00ALLINONE > logs/BUILD_KOA_LAST.hash && cat logs/BUILD_KOA_LAST.hash
echo "== the M1/M2 harnesses (the koa target's flags)"
HF="-std=gnu11 -O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off"
gcc $HF -o bin/gen3_m1 edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm || { echo "STOP: the M1 harness build failed"; exit 4; }
gcc $HF -o bin/gen3_m2 edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm || { echo "STOP: the M2 harness build failed"; exit 4; }
(cd bin && sha256sum gen3_m1 gen3_m2) > logs/BUILD_KOA_HARNESS.hash && cat logs/BUILD_KOA_HARNESS.hash
echo "== long double on this node: $(printf '#include <float.h>\n#include <stdio.h>\nint main(void){printf("LDBL_MANT_DIG %%d, sizeof %%zu\\n", LDBL_MANT_DIG, sizeof(long double));return 0;}\n' > /tmp/ld_$$.c && gcc -o /tmp/ld_$$ /tmp/ld_$$.c && /tmp/ld_$$)"
echo "BUILD OK: ./00ALLINONE (koa), bin/00ALLINONE_ld (koa-ld), bin/gen3_m1 and bin/gen3_m2 of $head"
