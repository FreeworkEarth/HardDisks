#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 a; registration 261012 sec. 4.8.5; the stage 1 gate package): THE BUILD-TIME RECORD of
# the stage 1 campaign build. Run ONCE in the sandbox shell right after cluster/gen3_koa_261009/build_gen3_koa.sh in the NEW clone
# ~/harddisks_gen3_s1 (engine-gen3 f72d764, fetched by hash; its build inputs are f72d764's by definition). It checks, then writes
# ~/harddisks_gen3_s1/hspist3/logs/BUILD_RECORD_S1.txt -- the commit and the SHA-256 of every binary the build made -- which every stage 1
# job checks WITHOUT git (s1_koa.sbatch, s1_gate_koa.sbatch, s1_xnode.sbatch). Git is used here, at build time, not at run time.
# The checks: the record does not exist yet (the campaign build is never rebuilt); the clone's HEAD is f72d764a8a7e837644c26bd8abc175d9e2cac268
# and no tracked file has changed; the four binaries exist, the two drivers print "00ALLINONE  git f72d764  target koa" / "koa-ld", and
# build_gen3_koa.sh's own hash files (logs/BUILD_KOA_LAST.hash, BUILD_KOA_LD.hash, BUILD_KOA_HARNESS.hash) hold the binaries' SHA-256.
#   bash ~/harddisks_gen3_k2/hspist3/cluster/gen3_s1_gate_261012/record_s1_build.sh
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
C="$HOME/harddisks_gen3_s1/hspist3"; REC="$C/logs/BUILD_RECORD_S1.txt"; WANT=f72d764a8a7e837644c26bd8abc175d9e2cac268
[ ! -e "$REC" ] || stop "$REC exists -- the campaign build is never rebuilt; nothing written"
cd "$C" || stop "no folder $C"
head=$(git -C "$C" rev-parse HEAD) || stop "git rev-parse failed in $C"
[ "$head" = "$WANT" ] || stop "the clone's HEAD is $head, not $WANT"
[ -z "$(git -C "$C" status --porcelain --untracked-files=no)" ] || stop "tracked files of the clone have changed: $(git -C "$C" status --porcelain --untracked-files=no | head -3 | tr '\n' ' ')"
for f in 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2; do [ -x "$f" ] || stop "$C/$f is missing"; done
V=$(./00ALLINONE --version | head -1); VL=$(bin/00ALLINONE_ld --version | head -1)
[ "$V" = "00ALLINONE  git f72d764  target koa" ] || stop "the driver's version line is '$V'"
[ "$VL" = "00ALLINONE  git f72d764  target koa-ld" ] || stop "the long-double driver's version line is '$VL'"
sha256sum --status -c logs/BUILD_KOA_LAST.hash || stop "00ALLINONE does not match logs/BUILD_KOA_LAST.hash"
(cd bin && sha256sum --status -c ../logs/BUILD_KOA_LD.hash && sha256sum --status -c ../logs/BUILD_KOA_HARNESS.hash) || stop "bin/ does not match the build's hash files"
{ echo "# the build-time record of the stage 1 campaign build (decision 39 item 5 a): written once by record_s1_build.sh after build_gen3_koa.sh"
  echo "commit $head"; echo "version $V"; echo "version_ld $VL"
  for f in 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2; do echo "sha256 $(sha256sum "$f" | cut -d' ' -f1) $f"; done
  echo "compiler $(gcc --version | head -1)"; echo "built $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"; } > "$REC.tmp$$" && mv "$REC.tmp$$" "$REC" || stop "cannot write $REC"
chmod a-w "$REC"
cat "$REC"
echo "BUILD RECORD OK: $REC"
