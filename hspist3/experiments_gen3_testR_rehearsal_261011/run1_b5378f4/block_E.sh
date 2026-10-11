squeue -u charing
sinfo -p shared -s
cd ~/harddisks_gen3_k2/hspist3 && sbatch --array=2 cluster/gen3_testR_261011/testR_koa.sbatch
