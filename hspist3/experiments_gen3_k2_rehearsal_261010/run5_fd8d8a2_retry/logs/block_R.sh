squeue -u charing
cd ~/harddisks_gen3/hspist3 && sbatch --array=2 --exclude=cn-05-08-07 cluster/gen3_koa_261009/testG_koa.sbatch
