#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G of the plan-author programme of sec. 4.7.12): PREPARED, NOT RUN. Build the gen3
# binaries on KOA from the clean clone ~/harddisks_gen3 (branch engine-gen3 at the commit named in the runsheet): the default
# `make koa` (-> ./00ALLINONE, recorded in logs/BUILD_KOA_LAST.hash as cluster/build_koa.sh does) and the long-double option
# `make koa-ld` (-> bin/00ALLINONE_ld, recorded in logs/BUILD_KOA_LD.hash). Both must carry "git <HEAD>" with no -dirty.
# Runs ONLY inside a Slurm job (home is noexec on the login node): from hspist3/ of ~/harddisks_gen3, inside
#   srun -p sandbox -A uh -c 2 --mem=4G -t 00:30:00 --pty /bin/bash
#   bash cluster/gen3_koa_261009/build_gen3_koa.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 2
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"; exit 2; }
ROOT_NAME=$(basename "$(dirname "$PWD")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: build from the gen3 clone ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "STOP: the clone has modified tracked files"; git status --short | head; exit 2; }
source cluster/koa_env.sh || exit 2
head=$(git rev-parse --short HEAD); mkdir -p logs bin
echo "== HEAD $head ($(git log -1 --format=%s | cut -c1-80)); branch $(git rev-parse --abbrev-ref HEAD)"
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
echo "== long double on this node: $(printf '#include <float.h>\n#include <stdio.h>\nint main(void){printf("LDBL_MANT_DIG %%d, sizeof %%zu\\n", LDBL_MANT_DIG, sizeof(long double));return 0;}\n' > /tmp/ld_$$.c && gcc -o /tmp/ld_$$ /tmp/ld_$$.c && /tmp/ld_$$)"
echo "BUILD OK: ./00ALLINONE (koa) and bin/00ALLINONE_ld (koa-ld) of $head"
