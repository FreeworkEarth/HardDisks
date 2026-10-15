##CHRIS 2026-10-12 (plan-author decision 39 items 5 a and 5 c; decision 20's rule): the Mac rehearsals of the stage 1 gate package.
- rehearse_s1_gate_mac_output_61a0a93.txt  cluster/gen3_s1_gate_261012/rehearse_s1_gate_mac.py at engine-gen3 61a0a93 (the gate job, the
                                           cross-node array and its compare, the STOP checks, the fetch and the gate verdict): PASS 12 of 12
- gate_job_log_rehearsal.out.txt, xnode_compare_log_rehearsal.out.txt   the two job logs of that run (block D prints these on KOA)
- rehearse_s1_build_mac_output_b63f650.txt, rehearse_s1_build_mac_output_688a377.txt   cluster/gen3_s1_gate_261012/rehearse_s1_build_mac.py:
                                           block B's build lines (build_s1_koa.sh, record_s1_build.sh) on a fresh f72d764 clone: PASS 8 of 8
                                           each (688a377: the record's first line names build_s1_koa.sh)
The rehearsal of the delivered commit's build lines is the first step of the stage 1 rehearsal (G4), which runs rehearse_s1_build_mac.py.
