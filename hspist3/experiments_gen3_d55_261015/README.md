# The consolidation's checks (g)5 and (g)6 (plan-author decision 55, addendum 2)

<!-- ##CHRIS 2026-10-15 -->
Run by CC on 2026-10-15 in Repo/HardDisks after the merge commit 0544b7e8 (engine-gen3 into main). The steps and their output,
verbatim: `0000_PLAN_OVERALL/consolidation_d51_261014/d55_consolidation_log_261015.txt`; the report: notes 261012 sec. 4.7.132.

- The build: `hspist3/experiments_gen3_d50_261014/bin_0544b7e8-dirty/` (build_d50.sh's command from the one tree; "-dirty" is only
  the kissfft nested clone with the project's local kiss_fft_log.h change, as on engine-gen3). Its build_record.txt lists every
  source's SHA-256; all 12 equal engine-gen3's (kissfft's against bin_a71f394's record).
- `tests/g5_T1/`: check (g)5, Part B's six registered stage 1 task lines, the new build against the last engine-gen3 build
  (experiments_gen3_d52_261014/bin_a71f394): RESULTS.md, every run's cmd.txt and stdout.log. EQUAL in 6 of 6.
- `e0/`: check (g)6, rule 4 (gen2 byte-identical to 7b08827, ctrl_min and ctrl_leg), the default build and every command plus
  --engine=gen2 (bin_engine_gen2.sh), by cluster/resched_gate_261005/audit_runs_261007.py; the 7b08827 references copied next to
  them from experiments_gen3_m3_261009/evidence_f42befb/e0/<tag>/<case>/ref (as stageA_evidence.py did). REPORT_rule4.txt:
  "plain vs ref" IDENTICAL in both cases, both builds. Committed: the report, every run's command.txt and version.txt.
  The run outputs and the reference copies stay untracked (data).
