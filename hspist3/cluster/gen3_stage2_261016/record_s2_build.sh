#!/bin/bash
# ##CHRIS 2026-10-16 (plan-author decision 58; notes 261012 sec. 4.7.142): THE BUILD-TIME RECORD of the S2 build. It is record_p2_build.sh
# with THESE CHANGES, all of them: the clone ~/harddisks_gen3_s2, the record logs/BUILD_RECORD_S2.txt, the build commit 8308fb00 (the package
# commit's parent; its version line "00ALLINONE  git 8308fb0  target koa"), the manifest s2_build_sources.sha256. Run ONCE inside the sandbox job
# of paste A0+A, from ~/harddisks_gen3_s2, right after the checkout of the package commit:
#   bash hspist3/cluster/gen3_stage2_261016/record_s2_build.sh
# The checks (P2's): the record does not exist yet; the binary exists, equals logs/BUILD_KOA_LAST.hash and prints the build commit's version
# line (clean: no -dirty; git may abbreviate the hash to more than 7 digits); every build input on disk has the SHA-256 that
# s2_build_sources.sha256 lists for 8308fb00 (made on the Mac with git show; kissfft at febd4ca); kissfft at febd4ca; no tracked file changed;
# build_koa.sh's own record of this job exists.
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
C="$HOME/harddisks_gen3_s2/hspist3"; REC="$C/logs/BUILD_RECORD_S2.txt"; WANT=8308fb004fdcb7e9c82197451cc04ab49ecd2f13
[ -n "${SLURM_JOB_ID:-}" ] || stop "run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"
[ ! -e "$REC" ] || stop "$REC exists -- the S2 build is never rebuilt; nothing written"
cd "$C" || stop "no folder $C"
source cluster/koa_env.sh || stop "the clone's cluster/koa_env.sh failed"
head=$(git -C "$C" rev-parse HEAD) || stop "git rev-parse failed in $C"
[ -z "$(git -C "$C" status --porcelain --untracked-files=no)" ] || stop "tracked files of the clone have changed: $(git -C "$C" status --porcelain --untracked-files=no | head -3 | tr '\n' ' ')"
[ "$(git -C kissfft rev-parse HEAD 2>/dev/null)" = febd4caeed32e33ad8b2e0bb5ea77542c40f18ec ] || stop "kissfft is not the pinned commit febd4ca"
[ -x 00ALLINONE ] || stop "$C/00ALLINONE is missing (run cluster/build_koa.sh first)"
sha256sum --status -c logs/BUILD_KOA_LAST.hash || stop "00ALLINONE does not match logs/BUILD_KOA_LAST.hash"
V=$(./00ALLINONE --version | head -1)
[[ "$V" =~ ^00ALLINONE\ \ git\ 8308fb0[0-9a-f]*\ \ target\ koa$ ]] || stop "the driver's version line is '$V', not '00ALLINONE  git 8308fb0  target koa'"
n=$(grep -c . cluster/gen3_stage2_261016/s2_build_sources.sha256)
ok=$(sha256sum -c cluster/gen3_stage2_261016/s2_build_sources.sha256 2>/dev/null | grep -c ': OK$')
[ "$ok" -eq "$n" ] || { sha256sum -c cluster/gen3_stage2_261016/s2_build_sources.sha256 2>&1 | grep -v ': OK$' | head -5; stop "build inputs: $ok of $n equal 8308fb00's"; }
BL="logs/BUILD_KOA_${SLURM_JOB_ID}.txt"
[ -r "$BL" ] || stop "no $BL: build_koa.sh did not run in this job"
{ echo "# the build-time record of the S2 build (decision 58): written once by record_s2_build.sh after cluster/build_koa.sh at 8308fb00"
  echo "commit_build $WANT"; echo "commit_scripts $head"; echo "version $V"
  echo "sha256 $(sha256sum 00ALLINONE | cut -d' ' -f1) 00ALLINONE"
  echo "build_inputs $ok of $n equal to 8308fb00's (cluster/gen3_stage2_261016/s2_build_sources.sha256)"
  grep -E '^(built_utc|host|cpu|compiler|make compiler|libs) ' "$BL" | sed 's/^/build_koa /'; } > "$REC.tmp$$" && mv "$REC.tmp$$" "$REC" || stop "cannot write $REC"
chmod a-w "$REC"
cat "$REC"
echo "BUILD RECORD OK: $REC"
