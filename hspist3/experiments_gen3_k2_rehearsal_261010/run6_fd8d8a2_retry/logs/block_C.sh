git --version
cd ~/harddisks_gen3_k2
git status --porcelain --untracked-files=no | wc -l
git fetch --depth 1 origin fd8d8a28037678f7cb3cf7c891ca52ecc60e67a9
git -c advice.detachedHead=false checkout fd8d8a28037678f7cb3cf7c891ca52ecc60e67a9
git log -1 --format='%h %s' | cut -c1-80
cd hspist3 && sha256sum cluster/gen3_koa_261009/merge_testG_koa.py cluster/gen3_koa_261009/merge_testG_k2.sbatch cluster/gen3_koa_261009/testG_chunk_logs_koa.sh cluster/gen3_koa_261009/test_merge_testG_plan_261010.py cluster/gen3_koa_261009/test_testG_chunk_logs_261010.py cluster/confinement_20261013/reduce_B.py
source cluster/koa_env.sh
mkdir -p logs/k2c_unit_tmp && TMPDIR=$PWD/logs/k2c_unit_tmp python3 cluster/gen3_koa_261009/test_merge_testG_plan_261010.py > logs/k2c_unit_merge.txt 2>&1; tail -1 logs/k2c_unit_merge.txt
TMPDIR=$PWD/logs/k2c_unit_tmp python3 cluster/gen3_koa_261009/test_testG_chunk_logs_261010.py > logs/k2c_unit_chunks.txt 2>&1; tail -1 logs/k2c_unit_chunks.txt
bash cluster/gen3_koa_261009/testG_chunk_logs_koa.sh ~/harddisks_gen3/hspist3/logs /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d24/rehearsal_fd8d8a2_b/lustre/harddisks_gen3/testG_koa/chunks $(( ($(wc -l < ~/harddisks_gen3/hspist3/cluster/gen3_koa_261009/tasks_testG_koa_261009.txt) + 24) / 25 )) > logs/k2c_chunk_check.txt 2>&1; echo "chunk check exit $?"; tail -25 logs/k2c_chunk_check.txt
python3 cluster/gen3_koa_261009/merge_testG_koa.py /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d24/rehearsal_fd8d8a2_b/lustre/harddisks_gen3/testG_koa --plan-only > logs/k2c_plan_only.txt 2>&1; head -60 logs/k2c_plan_only.txt
exit
