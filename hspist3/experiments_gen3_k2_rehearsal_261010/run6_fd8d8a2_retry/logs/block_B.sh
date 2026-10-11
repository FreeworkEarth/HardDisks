squeue -u charing
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
