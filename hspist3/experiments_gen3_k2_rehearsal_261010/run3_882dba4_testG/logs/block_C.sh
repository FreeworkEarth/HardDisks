git --version
cd ~/harddisks_gen3_k2
git status --porcelain --untracked-files=no | wc -l
git fetch --depth 1 origin 882dba47f3c2c6db28e3e067e8d113ef5205327b
git -c advice.detachedHead=false checkout 882dba47f3c2c6db28e3e067e8d113ef5205327b
git log -1 --format='%h %s' | cut -c1-80
cd hspist3 && sha256sum cluster/gen3_koa_261009/merge_testG_koa.py cluster/gen3_koa_261009/merge_testG_k2.sbatch cluster/gen3_koa_261009/test_merge_testG_plan_261010.py cluster/confinement_20261013/reduce_B.py
source cluster/koa_env.sh
mkdir -p logs/k2b_unit_tmp && TMPDIR=$PWD/logs/k2b_unit_tmp python3 cluster/gen3_koa_261009/test_merge_testG_plan_261010.py | tail -13
python3 cluster/gen3_koa_261009/merge_testG_koa.py /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d22/rehearsal_882dba4/lustre/harddisks_gen3/testG_koa --plan-only > logs/k2b_plan_only.txt 2>&1; head -60 logs/k2b_plan_only.txt
exit
