##CHRIS 2026-10-10 (plan-author decision 20, task 3; 261012 sec. 4.7.47): the text evidence of the Mac rehearsal of every KOA job script
(hspist3/cluster/gen3_koa_261009/rehearse_k2_mac.py). The data stayed in the session scratchpad (not copied here).
  run1_9d4490c/        every part at engine-gen3 9d4490c: the unit test, the gate job (a)-(d), K2-2 and its two STOP checks, the
                       cross-node array (2 pretended nodes) and compare, Test G stage 1, stage 2 and a chunk beyond the list, the
                       registered merge (it STOPPED: the fault), the fetch, the verdict (A as committed, B component check)
  run2_2988351_testG/  after the merge fix, at engine-gen3 2988351: the unit test and the whole Test G path again, with the merge from
                       the K2 clone (merge_testG_k2.sbatch); PASS
  harness_event_hash_compare_preliminary.txt  decision 20 task 6, PRELIMINARY: the 32 "A gen3, audits on" rows of KOA's harness
                       outputs (Chris's read-only grep, koa_harness_audits_on_rows_from_chris_grep_261010.txt) against the Mac's
Each run folder: rehearsal_driver_output.txt (the pass/fail table), logs/ (every job log, the make logs), shims/, the mini task list and
its SHA-256, the Mac builds' SHA file used for K2-2, the scratch tree listing, the disk use, SHA256_FILES.txt.
The job logs carry a .txt suffix after their original names (.out): the repository's .gitignore ignores *.out.
##CHRIS 2026-10-10 (plan-author decisions 21 and 22; 261012 sec. 4.7.49): the Test G path again, for the merge that plans before it moves
and runs without the array's job id (driver rehearse_k2_mac.py; the K2 clone cloned at 2988351, as on KOA, and updated in place by the
runsheet's K2-4 block C; all five K2-4 paste blocks run as written in runsheet_k2_261010.txt):
  run3_882dba4_testG/  mini list v2, REAL chunks of every block (registered chunks 1, 2, 33, 65, 97, 105, 121, 129 and the first 5 lines
                       of 130: 205 lines, 9 chunks); 14 of 15 PASS. It found the fault fixed in dc158b5: the merge summary's listing of
                       the run folders grew with the number of runs (block_C.out.txt: one line of 48 kB). The FAIL is verdict B, the
                       component check, on unpaired mini data (13 gen2 against 12 gen3 runs at M = 300 and 1500): see sec. 4.7.49.
                       merge_plan/ holds the merge's own plan files (block C's --plan-only, then the job's) and .merge_complete.
  run4_dc158b5_testG/  after the fix: the real chunks 1 and 2 (the run folders of one speed-of-sound folder per engine from two chunks);
                       PASS (14 of 14; verdict B not run: no paired runs at M = 300 and 1500 in this list).
