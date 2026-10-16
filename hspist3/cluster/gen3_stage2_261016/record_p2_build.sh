#!/bin/bash
# ##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE BUILD-TIME RECORD of the P2 build, as record_s1_build.sh
# wrote stage 1's (cluster/gen3_s1_gate_261012/record_s1_build.sh), for ONE fresh clone ~/harddisks_gen3_p2 that is built at 8aaf5fef
# by 8aaf5fef's own cluster/build_koa.sh (make -B koa, cluster/koa_env.sh: the stage 1 build's target, flags and environment) and then
# checked out at the package commit (runsheet block a). Run ONCE in the sandbox shell, from ~/harddisks_gen3_p2, right after that checkout:
#   bash hspist3/cluster/gen3_stage2_261016/record_p2_build.sh
# It checks, then writes ~/harddisks_gen3_p2/hspist3/logs/BUILD_RECORD_P2.txt (read-only), which every P2 job checks WITHOUT git
# (p2_koa.sbatch). Git is used here, at build time, not at run time. The checks:
#   - the record does not exist yet (the build is never rebuilt);
#   - the binary 00ALLINONE exists, equals build_koa.sh's hash file logs/BUILD_KOA_LAST.hash, and prints "00ALLINONE  git 8aaf5fe  target koa"
#     (the build commit, clean: no -dirty; git may abbreviate the hash to more than 7 digits);
#   - every build input on disk -- the tracked sources, headers and Makefile, cluster/build_koa.sh and cluster/koa_env.sh, the kissfft
#     sources -- has the SHA-256 that p2_build_sources.sha256 lists for 8aaf5fef (made on the Mac with git show; kissfft at febd4ca), so the
#     inputs the binary was built from are 8aaf5fef's and the package checkout changed none of them;
#   - kissfft is at febd4caeed32e33ad8b2e0bb5ea77542c40f18ec; no tracked file of the clone has changed;
#   - build_koa.sh's own record logs/BUILD_KOA_<this job>.txt exists (its compiler, CPU, node and time go into this record).
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
C="$HOME/harddisks_gen3_p2/hspist3"; REC="$C/logs/BUILD_RECORD_P2.txt"; WANT=8aaf5fef27c92391fc4fadd5133324740c29ef06
[ -n "${SLURM_JOB_ID:-}" ] || stop "run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"
[ ! -e "$REC" ] || stop "$REC exists -- the P2 build is never rebuilt; nothing written"
cd "$C" || stop "no folder $C"
source cluster/koa_env.sh || stop "the clone's cluster/koa_env.sh failed"
head=$(git -C "$C" rev-parse HEAD) || stop "git rev-parse failed in $C"
[ -z "$(git -C "$C" status --porcelain --untracked-files=no)" ] || stop "tracked files of the clone have changed: $(git -C "$C" status --porcelain --untracked-files=no | head -3 | tr '\n' ' ')"
[ "$(git -C kissfft rev-parse HEAD 2>/dev/null)" = febd4caeed32e33ad8b2e0bb5ea77542c40f18ec ] || stop "kissfft is not the pinned commit febd4ca"
[ -x 00ALLINONE ] || stop "$C/00ALLINONE is missing (run cluster/build_koa.sh first)"
sha256sum --status -c logs/BUILD_KOA_LAST.hash || stop "00ALLINONE does not match logs/BUILD_KOA_LAST.hash"
V=$(./00ALLINONE --version | head -1)
[[ "$V" =~ ^00ALLINONE\ \ git\ 8aaf5fe[0-9a-f]*\ \ target\ koa$ ]] || stop "the driver's version line is '$V', not '00ALLINONE  git 8aaf5fe  target koa'"
n=$(grep -c . cluster/gen3_stage2_261016/p2_build_sources.sha256)
ok=$(sha256sum -c cluster/gen3_stage2_261016/p2_build_sources.sha256 2>/dev/null | grep -c ': OK$')
[ "$ok" -eq "$n" ] || { sha256sum -c cluster/gen3_stage2_261016/p2_build_sources.sha256 2>&1 | grep -v ': OK$' | head -5; stop "build inputs: $ok of $n equal 8aaf5fef's"; }
BL="logs/BUILD_KOA_${SLURM_JOB_ID}.txt"
[ -r "$BL" ] || stop "no $BL: build_koa.sh did not run in this job"
{ echo "# the build-time record of the P2 build (decision 57): written once by record_p2_build.sh after cluster/build_koa.sh at 8aaf5fef"
  echo "commit_build $WANT"; echo "commit_scripts $head"; echo "version $V"
  echo "sha256 $(sha256sum 00ALLINONE | cut -d' ' -f1) 00ALLINONE"
  echo "build_inputs $ok of $n equal to 8aaf5fef's (cluster/gen3_stage2_261016/p2_build_sources.sha256)"
  grep -E '^(built_utc|host|cpu|compiler|make compiler|libs) ' "$BL" | sed 's/^/build_koa /'; } > "$REC.tmp$$" && mv "$REC.tmp$$" "$REC" || stop "cannot write $REC"
chmod a-w "$REC"
cat "$REC"
echo "BUILD RECORD OK: $REC"
