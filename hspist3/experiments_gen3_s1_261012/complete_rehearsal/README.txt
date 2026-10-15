##CHRIS 2026-10-13 (plan-author decision 47 item 3): the rehearsal of the reduction completion (cluster/gen3_s1_complete_261013/), the
contingency for D's merge (job 15055588). Only text is kept here; the rehearsal's trees stay in the session scratchpad.
- rc1_f434f13/   engine-gen3 f434f13b: stopped after the build. Its clone home/harddisks_gen3_rc was a full checkout: an LFS file outside
                 hspist3 failed, and it wrote 383 MB (the data cap of the message was broken by it: Report D47). Fixed in 60a744f3: the
                 clone as RC block B's (no checkout, K2's sparse set, then the commit).
- rc2_60a744f/   engine-gen3 60a744f3: COMPLETION REHEARSAL 11 of 13. K1 and K4 (killed inside reduce_s1) fail only the check that no file
                 equal to REF before the completion is written again. In both, reduce_s1 of the cut folder had written its first tables and
                 the completion ran it again (the logs: completion_K1, completion_K4), which rewrote those tables with the same content. The
                 tree equals REF file by file. rehearse_complete_mac.py now separates these from finished files (not yet rerun: the data cap).
                 logs/: every completion's log, the killed merges' logs (each ends at "FAILED lines in the chunk logs (logs/): 0": the
                 merge script's block-buffered lines are lost at a kill), the sbatch shim's record.
- test_rc_block_c.py, test_rc_block_c_output.txt   RC block C's guarded line with shims: PASS 7 of 7.
