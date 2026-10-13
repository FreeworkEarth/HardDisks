#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 a, option (i); registration 261012 sec. 4.8.5): THE STAGE 1 CAMPAIGN BUILD, made ONCE
# in the NEW clone ~/harddisks_gen3_s1 (engine-gen3 f72d764, fetched by hash). This is f72d764's cluster/gen3_koa_261009/build_gen3_koa.sh
# (byte-identical to a643bc7's, which built the gate binary) with THREE changes and nothing else: (1) it builds in ~/harddisks_gen3_s1 --
# f72d764's script refuses every clone not named harddisks_gen3, and ~/harddisks_gen3 (a643bc7) is not touched; (2) it is run from the
# scripts clone ~/harddisks_gen3_k2 and goes to the campaign clone itself; (3) it refuses when the clone's HEAD is not f72d764, or when a
# binary or the build-time record already exists (the campaign build is never rebuilt). The make targets (make -B koa-ld, make -B koa),
# the harness flags, the hash files (logs/BUILD_KOA_LD.hash, BUILD_KOA_LAST.hash, BUILD_KOA_HARNESS.hash), the environment (the clone's
# own cluster/koa_env.sh, identical between a643bc7 and f72d764) and the checks are build_gen3_koa.sh's lines, unchanged.
# Then record_s1_build.sh writes the build-time record. Runs ONLY inside a Slurm job (home is noexec on the login node), at the
# compute-node prompt of
#   srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
#   bash ~/harddisks_gen3_k2/hspist3/cluster/gen3_s1_gate_261012/build_s1_koa.sh 2>&1 | tail -20
set -uo pipefail
cd "$HOME/harddisks_gen3_s1/hspist3" || { echo "STOP: no campaign clone ~/harddisks_gen3_s1"; exit 2; }
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"; exit 2; }
ROOT_NAME=$(basename "$(dirname "$PWD")")
[ "$ROOT_NAME" = harddisks_gen3_s1 ] || { echo "STOP: build from the campaign clone ~/harddisks_gen3_s1, not $ROOT_NAME"; exit 2; }
[ "$(git rev-parse HEAD)" = f72d764a8a7e837644c26bd8abc175d9e2cac268 ] || { echo "STOP: the campaign clone's HEAD is $(git rev-parse HEAD), not f72d764a8a7e837644c26bd8abc175d9e2cac268"; exit 2; }
for f in logs/BUILD_RECORD_S1.txt 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2; do [ ! -e "$f" ] || { echo "STOP: $PWD/$f exists -- the campaign build is made once, never rebuilt"; exit 2; }; done
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
