squeue -u charing
cd ~/harddisks_gen3_k2/hspist3 && sbatch --array=2 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_testR_261011/testR_koa.sbatch
