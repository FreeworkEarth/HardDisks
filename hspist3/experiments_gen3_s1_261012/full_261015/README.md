# Stage 1, Part B: the full report (plan-author decision 54, 2026-10-15)

<!-- ##CHRIS 2026-10-15 -->
All six groups fetched (B100, B400, D, B900, S, B1600 in `../../experiments_gen3_stage1_koa_261012/`); everything here is read from there, read only. The notes record it: Report D54-1 (sec. 4.7.123, the verdicts), D54-2 (4.7.124, B1600's INTERIM details, in `../interim_261013/B1600/`) and D54-3 (4.7.125, the figures).

- `a_verdict_full.txt`: `s1_verdict_261012.py --data experiments_gen3_stage1_koa_261012`, no --interim. The script as committed at gate G5 (b63f6502). VERDICT K: KR-INCONSISTENT (eta 0.60); Verdict S: START-INDEPENDENT in 20 of 20 cells.
- `b_verdict_noflag.txt`: `s1_verdict_noflag_261015.py`, what the verdict script's report() leaves out (Verdict K's z per density, the flags, every pooled result also without the flagged cells and densities), from the script's own analyse() and functions.
- `figures_261015b/`: the registered figures as reported (D54-3): eight figures, each as PNG and PDF with the CSV of the numbers drawn, by `s1_figures_261015.py`.
- `figures_261015/`: the first drawing of the same figures by an earlier layout of that script. Same numbers (the CSVs are identical), kept as it was. Its residual panels put the approach cells off scale, region names overlapped a note, and the per-trajectory psi6 dots hid under the cell markers. Only the layout changed in `figures_261015b/`.
- `c_figures.txt`, `c_figures_b.txt`: the two figure runs' output (the checks and the files written).
- `RUN_RECORD.txt`: every command, its time, exit code and seconds, the Python and library versions, and the SHA-256 of the scripts used.
