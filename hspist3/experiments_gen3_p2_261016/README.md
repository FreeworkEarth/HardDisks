# experiments_gen3_p2_261016: the Mac rehearsal of the P2 package (plan-author decision 57; Report D57-1)

##CHRIS 2026-10-15. The package is hspist3/cluster/gen3_stage2_261016/ at commit c69804b0 (the runsheet, rehearse_p2_mac.py and the openrsync
fix of fetch_p2_koa.sh in the commit after it). Nothing here is P2 data: the four runs are the rehearsal's, on the Mac.

| folder | what |
|---|---|
| rehearsal_261015c | THE REHEARSAL OF RECORD: 12 of 12 steps PASS (REHEARSAL.md). The runsheet's blocks A, B (twice), C, D, D2, the merge again, E, as typed, each block's transcript (block_*.txt: every line as "$ <line>" and its output); the fetched merged tables without the run folders (fetched_merged/); p2_tables_mac.py on the fetch (p2_tables/); the fetch filter's test with the Mac's openrsync (fetch_filter_test_openrsync.txt); the second test, after the gzipped traces were added to the fetch for the P2 report's trace figure (decision 57 addendum 1; fetch_filter_test2_openrsync.txt, on a copy of the merged root with two decoys). |
| rehearsal_261015b | the second attempt, kept: 9 of 11. Both failures were the rehearsal script's own: its block-A check read a transcript whose commands came after their outputs (bash -v), and the runsheet's block D held a "  ... later" note that ended the block before its log line (block D was then split into D and D2). The package ran clean there too. |
| rehearsal_261015 | the first attempt, empty: it stopped before any block, because rehearse_k2_mac.kissfft() clones kissfft from the locked worktree, whose clone the consolidation moved into the main checkout. |

The rehearsal's scratch roots (the clones, the build, the Lustre stand-in) are in CC's scratchpad, d57/rehearsal_1 .. rehearsal_3. The four task
lines rehearsed are the package's lines 221, 271, 321, 371 (eta 0.78, alpha 5, s = 0; G1 and G2, pin 0 and turned), unchanged, at full length.
