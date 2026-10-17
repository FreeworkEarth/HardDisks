##CHRIS 2026-10-16 (plan-author decision 58 item 3; Report D58-1, notes 261012 sec. 4.7.143)

# The Mac rehearsals of the S2 package (914db066, built at its parent 8308fb00)

Each folder is one run of `cluster/gen3_stage2_261016/rehearse_s2_mac.py --pkg 914db066e95205b1cb4f253d32827600a1a719b5`: the runsheet's
own pastes as Chris types them, with the shims and substitutions the script lists, on four task lines at full length (N 900 turned
0.76 k697, N 900 pin 0 0.78 k685, N 1600 0.72 k866 and 0.80 k1024; alpha 5, s 0). The KOA side of the package (the sbatch files, the
worker, the task list, the chunk table, the manifests, the build record script, the fetch) is 914db066 in all three; only the Mac-side
files changed between them (the runsheet, the rehearsal script, the tables script). Nothing was deleted (decision 40): each run went
into a new folder.

| folder | runsheet | rehearsal script | tables script | result | of record |
|---|---|---|---|---|---|
| rehearsal_261016/ | v1 (runsheet_s2_261017_as_rehearsed.txt) | v1 (rehearse_s2_mac_as_run.py) | v1 (s2_tables_mac.py at 914db066) | 12 of 12 PASS | no |
| rehearsal_261016_b/ | v2, after the first review | v2 | v2 | 16 of 16 PASS; the check of the fixes found that four steps (check 0 among them) could pass on the echoed command text | no |
| rehearsal_261016_c/ | v3, after the check of the fixes (= the committed runsheet_s2_261017.txt) | v3 (= the committed one; checks read printed lines only) | v3 (= the committed one) | see its REHEARSAL.md | YES |

The four runs are byte-identical in all three (the same events, engine hashes and uncompressed output SHA-256).
Each folder: REHEARSAL.md (the step table), rehearsal_stdout.txt (the script's output), one transcript per paste or check
(every command as `$ <line>`, then the output), fetched_merged/ (the fetched tables, without the traces and the run folders), s2_tables/
(s2_tables_mac.py on the fetched folder) and s2_tables_stdout.txt. The scratch roots ($SP/d58/rehearsal_1, _2, _3: clones, the local
Lustre tree, the build) stay in the session scratchpad.
