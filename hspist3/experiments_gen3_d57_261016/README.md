# experiments_gen3_d57_261016: plan-author decision 57, addendum 3 (Report D57-3; notes 261012 sec. 4.7.140)

##CHRIS 2026-10-16. Mac only; no KOA. The P2 data are KOA's, fetched into ../experiments_gen3_stage2_koa_261016/P2 (not tracked); the one Mac
run is here (data/). RUN_RECORD.txt holds the run's command, binary and file hashes, and the scripts' hashes.

| folder or file | what |
|---|---|
| tasks_d57.txt, data/, worker_0.out | THE ONE MAC RUN (task 4): P2's task line 21 (eta 0.72, G1 H 30, pin 0, alpha 5, s 0) through cluster/gen3_stage2_261016/p2_worker_koa.sh, unchanged except T_snap 10 (snapshots every 10 sigma-time over the whole run: the binary has no time window for them). Binary: experiments_gen3_d50_261014/bin_0544b7e8-dirty/00ALLINONE (D56-1's Mac build; 11 of its 12 sources equal P2's 8aaf5fef build inputs, kiss_fft_log.h differs by an #ifndef DEBUG guard). |
| scripts/d57_c_vs_e.py -> c_vs_e/ | task 1: E1 per cell against the relaxed strain e (the record's mean, both compartments), the 8 P2 cells and D56-1's 0.76 pilot; the slopes d ln c / d e. |
| scripts/d57_bistable.py -> bistable/ | task 2: the bistable run (0.78 G2 pin 0, alpha 2, run 9) against its nine siblings; the (cell, mass) E1 with and without it; the scan of all 800 compartments. |
| scripts/d57_geometry_cost.py -> geometry_cost/ | task 3: the registration's geometry table (12 cells turned, 2 cells pin 0; three heights each) and the design's cost. |
| scripts/d57_reorientation.py -> reorientation/ | task 4: the reorientation at 0.72 in P2 (psi6_t and snapshots of the two pin-0 cells) and frame by frame in the Mac run. |
| c_vs_e_c/, reorientation_c/ | THE OUTPUTS OF RECORD for tasks 1 and 4: the scripts after the independent verification (verification/). d57_c_vs_e.py adds a sensitivity row (point 1's e error = half the two G1 cells' e difference) and the G1 note; d57_reorientation.py reads the snapshot interval from the data (67-117 at 0.72, not 61-117) and keeps the Mac table's t = 150 row. Their data tables equal the first run's. |
| c_vs_e/, reorientation/ | the first run of tasks 1 and 4, before the verification (kept). |
| c_vs_e_b/, reorientation_b/ | a rerun of the UNCHANGED scripts after an edit that failed to compile (a quoting slip): identical to c_vs_e/ and reorientation/ (diff -r); kept, nothing is deleted. |
| verification/ | the independent verification: workflow d57-3-verify (4 sub-agents, read-only on the repo, their own code; workflow_d57-3-verify_script.js), its result (workflow_d57-3-verify_result.json: 89 checks, 89 agree) and each agent's scratch code and outputs (v1_c_vs_e, v2_bistable, v3_geometry, v4_reorientation), copied from the scratchpad. |
| *_stdout.txt | each script's printed output. |
