git --version
cd ~/harddisks_gen3_k2
git status --porcelain --untracked-files=no | wc -l
git fetch --depth 1 origin c70d314e927580e35267045d58a22fdf10fb9850
git -c advice.detachedHead=false checkout c70d314e927580e35267045d58a22fdf10fb9850
git log -1 --format='%h %s' | cut -c1-80
cd hspist3 && sha256sum cluster/gen3_preview_261011/preview_koa.sbatch cluster/gen3_preview_261011/preview_merge.sbatch cluster/gen3_preview_261011/preview_stopped_chunks.sh cluster/gen3_preview_261011/test_preview_stopped_chunks.py cluster/gen3_testR_261011/test_merge_koa_plan.py cluster/gen3_testR_261011/test_chunk_logs_koa.py cluster/gen3_preview_261011/preview_inputs.sha256 cluster/gen3_preview_261011/preview_merge_inputs.sha256
sha256sum -c cluster/gen3_preview_261011/preview_inputs.sha256 cluster/gen3_preview_261011/preview_merge_inputs.sha256 | grep -c ': OK$'
sha256sum -c cluster/gen3_testR_261011/testR_inputs.sha256 cluster/gen3_testR_261011/testR_merge_inputs.sha256 | grep -c ': OK$'
(cd ~/harddisks_gen3/hspist3 && sha256sum -c ~/harddisks_gen3_k2/hspist3/cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256)
cat ~/harddisks_gen3/hspist3/logs/BUILD_KOA_LAST.hash
~/harddisks_gen3/hspist3/00ALLINONE --version | head -1
source cluster/koa_env.sh
mkdir -p logs/p0_unit_tmp && TMPDIR=$PWD/logs/p0_unit_tmp python3 cluster/gen3_testR_261011/test_merge_koa_plan.py > logs/p0_unit_merge.txt 2>&1; tail -1 logs/p0_unit_merge.txt
TMPDIR=$PWD/logs/p0_unit_tmp python3 cluster/gen3_testR_261011/test_chunk_logs_koa.py > logs/p0_unit_chunks.txt 2>&1; tail -1 logs/p0_unit_chunks.txt
TMPDIR=$PWD/logs/p0_unit_tmp python3 cluster/gen3_preview_261011/test_preview_stopped_chunks.py > logs/p0_unit_stopped.txt 2>&1; tail -1 logs/p0_unit_stopped.txt
(cd validation && python3 -c "import paper1_populate_cs_err_20261002, resched_testT_followup_261007, tests_20260913; print('the E0 and nu_d imports: OK')")
ls /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/p/rehearsal_run2/lustre/harddisks_gen3/preview_koa 2>&1 | head -2
exit
