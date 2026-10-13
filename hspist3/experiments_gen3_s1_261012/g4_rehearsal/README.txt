##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.10, gate G4): the Mac rehearsal of stage 1, Part B
(cluster/gen3_stage1_261012/rehearse_s1_mac.py): the build rehearsal, the steps after the gate, a mini list of 37 registered task lines in 7
chunks (groups B100, D, S; chunks 3 and 4 into one folder), and the mini runsheet's blocks as typed (make_s1_runsheet.py on the mini table).
- run2_0969b0e/   engine-gen3 0969b0ed: REHEARSAL G4: PASS (24 of 24). rehearse_s1_mac_output.txt (the table), runsheet_stage1_mini.txt (the
                  blocks it ran), logs/ (every block's output, the three merge job logs, the verdict's --interim output, the sbatch shim's
                  record), SHA256_FILES.txt (every file of its Lustre stand-in and of the fetched folders; the data stay in the scratchpad).
- run1_631d8e3/   engine-gen3 631d8e3c: FAIL (19 of 24); every failure a fault of the rehearsal's own checks (macOS wc -l pads; the
                  expected counts were read from the echo command; block E counted the refused merges' exits). Its one finding for
                  the runsheet: block H's filter missed reduce_B's last folder line (a trailing space); fixed in 0969b0ed.
