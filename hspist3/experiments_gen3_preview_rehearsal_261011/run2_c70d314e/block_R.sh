squeue -u charing
cd ~/harddisks_gen3_k2/hspist3 && sbatch --array=2 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_preview_261011/preview_koa.sbatch
