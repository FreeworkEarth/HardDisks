#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.29; plan-author decision 13, item 2c; repeated on the gate commit itself by decision 15, item 1a,
# sec. 4.7.33: then the commit argument IS the runsheet's, and the STEP 1 line says so): the MAC DRY RUN of the runsheet's steps 1 and 2. A fresh
# clone of the PUSHED engine-gen3 from GitHub into <scratch>/harddisks_gen3; the runsheet's step-1 commands as written, with one
# substitution: the commit is the pushed engine-gen3 head given as the argument, because the registration commit a643bc7 is pushed only
# after Report K. Then, as information, a fetch by hash of a commit that is NOT the branch tip (the head's parent): the runsheet relies
# on that once engine-gen3 has moved on. Then both driver builds with the Mac equivalents of the koa targets (make release-ld and make
# release: the koa flags without -march=x86-64-v2, which arm64 does not take), the M1/M2 harnesses with build_gen3_koa.sh's flags (the
# same without -march/-mtune), and ctrl_min with both builds through the clone's own audit_runs_261007.py, against the 7b08827
# reference outputs of ctrl_min ("plain vs ref" must be IDENTICAL). Network: read-only clones from GitHub. Nothing outside <scratch>
# is written.
# usage: bash hspist3/cluster/gen3_koa_261009/dryrun_mac.sh <empty scratch dir> <pushed engine-gen3 commit, 40 hex> <its parent, 40 hex> <ctrl_min reference dir>
set -uo pipefail
S=${1:?scratch dir}; C=${2:?pushed commit}; P=${3:?a non-tip commit}; REFD=${4:?ctrl_min reference outputs}
[ -d "$S" ] && [ -z "$(ls -A "$S")" ] || { echo "STOP: $S must be an existing empty folder"; exit 2; }
export GIT_TERMINAL_PROMPT=0
echo "== dry run $(date '+%Y-%m-%d %H:%M:%S %Z') in $S; pushed engine-gen3 per ls-remote: $(git ls-remote https://github.com/FreeworkEarth/HardDisks.git engine-gen3 | cut -f1)"
G=$(grep -o 'git fetch --depth 1 origin [0-9a-f]\{40\}' "$(dirname "$0")/runsheet_gen3_koa_261009.txt" | head -1 | awk '{print $NF}')   # the runsheet's gate commit
if [ "$C" = "$G" ]; then echo "== STEP 1, as written: the runsheet's gate commit $C itself, fetched by hash"
else echo "== STEP 1, as written, with $C standing in for the runsheet's gate commit $G (not fetchable yet)"; fi
cd "$S" || exit 2
set -x
git clone --depth 1 --branch engine-gen3 --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_gen3
cd harddisks_gen3
git fetch --depth 1 origin "$C"
git -c advice.detachedHead=false checkout "$C"
git sparse-checkout set --no-cone "/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"
git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft
git -C hspist3/kissfft -c advice.detachedHead=false checkout febd4caeed32e33ad8b2e0bb5ea77542c40f18ec
git log -1 --format='%h %s' | cut -c1-80
set +x
echo "== information: a fetch by hash of a commit that is not the branch tip ($P), in a second fresh clone"
( cd "$S" && git clone -q --depth 1 --branch engine-gen3 --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git nontip \
  && cd nontip && git fetch -q --depth 1 origin "$P" && git -c advice.detachedHead=false checkout -q "$P" && echo "   checked out $(git rev-parse HEAD)" ) \
  || echo "   **the non-tip fetch failed**"
cd "$S/harddisks_gen3/hspist3" || exit 2
echo "== the clone: HEAD $(git rev-parse HEAD); tracked files modified: $(git status --porcelain --untracked-files=no | wc -l | tr -d ' '); kissfft $(git -C kissfft rev-parse HEAD); files checked out: $(git ls-files | wc -l | tr -d ' ') ($(du -sh . | cut -f1) incl. kissfft)"
echo "== STEP 2, Mac equivalents: make -B release-ld (koa-ld), make -B release (koa)"
mkdir -p logs bin
make -B release-ld > logs/make_release_ld.txt 2>&1 || { echo "STOP: make release-ld failed"; tail -5 logs/make_release_ld.txt; exit 4; }
./00ALLINONE --version | head -1; mv ./00ALLINONE bin/00ALLINONE_ld
make -B release > logs/make_release.txt 2>&1 || { echo "STOP: make release failed"; tail -5 logs/make_release.txt; exit 4; }
./00ALLINONE --version | head -1
HF="-std=gnu11 -O2 -ffp-contract=off"
cc $HF -o bin/gen3_m1 edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm 2> logs/m1_build.txt && echo "   M1 harness built ($(grep -c 'warning:' logs/m1_build.txt) warnings)" || echo "   **the M1 harness build failed**"
cc $HF -o bin/gen3_m2 edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm 2> logs/m2_build.txt && echo "   M2 harness built ($(grep -c 'warning:' logs/m2_build.txt) warnings)" || echo "   **the M2 harness build failed**"
shasum -a 256 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2 | sed 's/^/   /'
echo "== ctrl_min with both builds (the clone's audit_runs_261007.py), against the 7b08827 reference outputs"
for tag in release release-ld; do
  b="$PWD/00ALLINONE"; [ "$tag" = release-ld ] && b="$PWD/bin/00ALLINONE_ld"
  O="$S/ctrl_min_$tag"
  python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin "$b" --out "$O" --cases ctrl_min --jobs 2 > "$S/ctrl_min_$tag.runner.txt" 2>&1 \
    || { echo "   **ctrl_min ($tag) failed**"; continue; }
  cp -R "$REFD" "$O/ctrl_min/ref"
  python3 cluster/resched_gate_261005/audit_runs_261007.py report --out "$O" > "$O/report.txt" 2>&1
  grep -E "^\| ctrl_min \| git" "$O/report.txt" | awk -F'|' -v t="$tag" '{print "   " t ": version" $3 "| audit vs plain:" $(NF-2) "| plain vs ref:" $(NF-1)}'
  grep -h "^\[EDMD-CONTACT\]" "$O/ctrl_min/plain/run.log" | sed 's/^/   /'
done
echo "== dry run done $(date '+%Y-%m-%d %H:%M:%S %Z')"
