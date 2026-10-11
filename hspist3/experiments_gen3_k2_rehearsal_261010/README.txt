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
