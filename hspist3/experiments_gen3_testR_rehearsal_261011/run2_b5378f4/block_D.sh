squeue -u charing
cd ~/harddisks_gen3_k2/hspist3
cat logs/gen3-testR_*_1.out
for e in gen2 gen3; do for M in 50 100 300; do c=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d28/rehearsal_run2/lustre/harddisks_gen3/testR_koa/chunks/chunk_1/testR/$e/epi8_H_H10_L10/m_$M; echo "$e M=$M traces $(ls $c/wall_x_positions_* 2>/dev/null | wc -l) runs $(grep -c '^##RUN' $c/run.log) gen3-records-clean $(grep -c '^\[EDMD3-HEALTH\] .*: clean=1 ' $c/run.log) build '$(cat $c/.build_git)'"; done; done
