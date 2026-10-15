##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.10, gate G4): the Mac rehearsal of stage 1, Part B
(cluster/gen3_stage1_261012/rehearse_s1_mac.py): the build rehearsal, the steps after the gate, a mini list of 37 registered task lines in 7
chunks (groups B100, D, S; chunks 3 and 4 into one folder), and the mini runsheet's blocks as typed (make_s1_runsheet.py on the mini table).
- run2_0969b0e/   engine-gen3 0969b0ed: REHEARSAL G4: PASS (24 of 24). rehearse_s1_mac_output.txt (the table), runsheet_stage1_mini.txt (the
                  blocks it ran), logs/ (every block's output, the three merge job logs, the verdict's --interim output, the sbatch shim's
                  record), SHA256_FILES.txt (every file of its Lustre stand-in and of the fetched folders; the data stay in the scratchpad).
- run1_631d8e3/   engine-gen3 631d8e3c: FAIL (19 of 24); every failure a fault of the rehearsal's own checks (macOS wc -l pads; the
                  expected counts were read from the echo command; block E counted the refused merges' exits). Its one finding for
                  the runsheet: block H's filter missed reduce_B's last folder line (a trailing space); fixed in 0969b0ed.
- run3_2c871b5/   engine-gen3 2c871b52 (plan-author decision 42 item 5): the runsheet's submit blocks guarded (every E, every G, block F's
                  resubmission line); each typed once (it submits) and a second time (it refuses); E-B100 after chunk 1; group D's
                  block E also with a pretended pending job and with a failing squeue; group S's guarded resubmission line once and
                  again; B100's merge job still refuses a second submission without its guard: REHEARSAL G4: PASS (34 of 34).
- manifest_check_at_2c871b5_output.txt   sha256sum -c of both manifests on the files of 2c871b52: 21 OK (decision 42 item 5 c)
- runsheet_command_diff_88d2eaf0_6398ae1e.txt   the stage 1 runsheet's command lines, 88d2eaf0 against 6398ae1e: 18 of 34 blocks
                  changed (E, F, G of each group), each holding its old command(s) unchanged and in their order (decision 42 item 5 d)
- run4_1bb9ae8/   ##CHRIS 2026-10-13 (plan-author decision 45 addendum 2): engine-gen3 1bb9ae8b, blocks F and G as ONE paste (block FG) for
                  D and S of the mini table (and B900, B1600 of the full one): REHEARSAL G4: PASS (40 of 40). Group D's block FG typed while F
                  is not clean -- an extra log of its chunk that did not finish (logs != chunks), its chunk's log without the done line (done
                  logs != chunks; the chunk check fails too), a FAILED line, its chunk folder moved aside (the chunk check exits 2) -- each
                  prints "F-D NOT CLEAN -- the merge was not submitted" and submits nothing (every change undone by moving files back); then
                  clean (the merge submitted inside it, exit 0) and a second time (G's own guard refuses: merged/ exists). Group S's block FG
                  prints NOT CLEAN after the faulty attempt (the resubmission lines; no "nothing to resubmit") and, by design, again after
                  the resubmission (4 logs for 2 chunks); G's guarded line alone (the then-branch of FG's last line) then stands in for the
                  plan author's go: once (exit 0) and a second time (refused); FG-S a second time: NOT CLEAN. Group B100 as in run 3.
- runsheet_diff_6398ae1e_a3eb77e4.txt   the stage 1 runsheet, 6398ae1e against a3eb77e4 (make_s1_runsheet.py at 1bb9ae8b, --commit 7bdfb67e):
                  blocks F and G of D, B900, S, B1600 replaced by block FG (F's command lines unchanged and first; G's guarded line inside
                  FG's last line byte for byte as before) and two header lines; every other block byte for byte unchanged
- run5_c1b2016/   ##CHRIS 2026-10-13 (plan-author decision 47 item 2): engine-gen3 c1b20169, the runsheet with the merges of B900, S and B1600 at
                  --time=24:00:00, the queue guards of E-S (no gen3-s1-B900 job) and E-B1600 (no gen3-s1-S or gen3-s1-B900 job), and the plan
                  author's timing line in every FG: REHEARSAL G4: PASS (41 of 41). The squeue shim filters a pretended queue by -n; E-S is
                  refused while a pretended gen3-s1-B900 array job is queued and submits beside a pretended gen3-s1-merge-B900 job; FG-D and
                  FG-S print the timing line; D's merge is submitted without --time, S's (G-S's line alone, after the NOT CLEAN by design)
                  with --time=24:00:00.
- runsheet_diff_a3eb77e4_19101796.txt, runsheet_command_diff_a3eb77e4_19101796.txt   the stage 1 runsheet a3eb77e4 against 19101796 (all
                  lines; command lines only): FG-D, FG-B900, FG-S, FG-B1600, E-S and E-B1600 changed, every other block unchanged
