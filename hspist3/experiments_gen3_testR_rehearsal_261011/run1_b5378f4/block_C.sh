               both engines; 1 core, at most 1 h on shared; about 12 min):
squeue -u charing
cd ~/harddisks_gen3_k2/hspist3 && sbatch --array=1 cluster/gen3_testR_261011/testR_koa.sbatch
