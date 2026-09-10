# Sound-speed structural analysis: completion and Claude Code handoff

Completed 2026-09-09. The output directory retains the date when this package was started (2026-09-08).

## Scope and status

**Completed:** existing-data structural screening, explicit exclusions, reproducible tables, four diagnostic figures in PNG and vector PDF, and permanent tests.

**Not completed or claimed:** final sound-speed validation/refits, temperature normalization, measured ψ4, initialization/hold-length experiments, equilibrium-divider fluctuations, or Paper 2's first-law gate. This closes the structural-screening part of Prompt 2, not all of Prompt 2 or Paper 1. No new simulation or model call was launched. No original data, physics code, existing plot, NeurAIpil file or research note was modified. Nothing was staged, committed or pushed.

Repository root:
`/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks`

Output directory relative to that root:
`hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/analysis_paper1_20260908/`

Open its `README.md` for definitions and counts. Start visually with `02_structure_change.png` and `03_routeB_distributions.png`.

## What the data establish

Read **18,938 unique trajectory records from 390 leaf datasets**. Merged copies are deliberately not counted again. The selected sources are route A, route B, strip ladders N=100/200/400, overnight N=1000, and fixed-aspect families A/C. The accelerated-core validation campaign is excluded from the selection.

The analysis retained **18,421 records** in paired structural summaries, comprising **122 density/size cells and 1,022 mass-resolved cells**. It excluded 517 distinct records from these summaries:

- **505 with explicit nonzero health warnings**, even though historical batch ledgers marked them valid.
- **25 with undefined ψ6 because no neighbours were present**, including 13 already excluded for health; hence 12 additional omissions. These missing structural values are not evidence of failed dynamics. They remain blank in the trajectory export, not zero.

### Important correction to previous status notes: health warnings

The single overnight clamp warning was not the only warning in the selected data:

| Series | Nominal η | Records with nonzero health |
|---|---:|---:|
| Route A | 0.019635 | 162 |
| Route A | 0.026180 | 180 |
| Route A | 0.039270 | 162 |
| Overnight N=1000 | 0.700 | 1 |

The first 504 have `wall_overdue=1` or `wall_overdue=2`. The overnight case is M=750, repeat=0, seed=2687575689, `wall_clamp_repairs=1`. Its original log is `overnight_N1000_20260826/eta_0p700/m_750/run.log`, line 6. All other selected families have no nonzero warning matched in their leaf logs.

**Do not infer the cause or size of a physical bias from the warning alone.** This package applies the already adopted strict exclusion rule. The meaning and provenance of the dilute overdue-wall warnings need a focused follow-up before recycling the previous low-density sound-speed fits. Nothing has been deleted, and this analysis has not yet re-extracted FFT frequencies or changed existing sound-speed summaries.

### Structural changes that matter for the paper

Values below are mean global |ψ6| across retained repeats and masses. Delta is the paired end-minus-hold change; its reported SEM is descriptive across those trajectories, not a confidence interval for an identical-condition population.

| Series | N | η (nominal) | Hold | End | Mean change ± SEM |
|---|---:|---:|---:|---:|---:|
| Route A | 100 | 0.72 | 0.90849 | 0.89779 | −0.01069 ± 0.00633 |
| Route B | 100 | 0.72 | 0.25890 | 0.69181 | +0.43292 ± 0.01730 |
| Fixed-aspect A | 1600 | 0.67 | 0.13286 | 0.15698 | +0.02412 ± 0.01961 |
| Fixed-aspect A | 1600 | 0.70 | 0.35107 | 0.46083 | +0.10976 ± 0.02087 |
| Fixed-aspect A | 1600 | 0.71 | 0.42200 | 0.52241 | +0.10040 ± 0.01907 |
| Fixed-aspect A | 1600 | 0.72 | 0.45882 | 0.53319 | +0.07437 ± 0.01687 |

Route B near 0.72 changes substantially during measurement. Some fixed-aspect runs also change, including N=1600 at 0.70 and 0.71 under the pre-existing absolute-change >0.1 screen. Larger N does not by itself establish stationarity. Conversely, a small mean change is not proof of equilibrium; opposing individual changes can cancel. The tables therefore also record the fraction of trajectories with absolute change >0.1.

The two-tail distribution screen is explicitly a heuristic, not a formal bimodality test. Low neighbour counts are flagged only at η≥0.6 and do not prove square symmetry. At very low density ψ6 is a poor structural diagnostic because too few neighbours lie inside the cutoff. The source neighbour count averages only over particles with neighbours; the source global ψ6 normalization uses all particles.

## Deliverables

```text
hspist3/validation/
  analyze_sound_structure.py
  test_sound_structure.py

.../00_eta_sweep_ROMAN/analysis_paper1_20260908/
  README.md
  trajectories.csv
  exclusions.csv
  structural_cells.csv
  structural_by_mass.csv
  source_ledgers.csv
  source_manifest.json
  01_structure_endpoints.png / .pdf
  02_structure_change.png / .pdf
  03_routeB_distributions.png / .pdf
  04_mass_dependence.png / .pdf

0000_PLAN_OVERALL/ALL_MARKDOWNS/
  260909_speedsound_structure_handoff_CODEX.md
```

The data-analysis workflow preserves raw observations, separates missing values from failures, exports explicit flags, and retains source paths/row numbers/seeds. The manifest fingerprints all 1,560 consumed CSV, command, ledger and log files, plus the analysis script. It is not a fingerprint of every simulation trace or the historical simulator binary.

## Verification and reproduction

- 19 permanent tests passed: paired statistics, missing-order handling, range checks, mass separation, source provenance, duplicate detection, ledger consistency, health-to-seed matching, accelerated-core exclusion and source-change detection.
- Compilation passed for both new Python files.
- All 1,560 input hashes were checked unchanged after generation.
- Independent calculation from Route B's original η=0.72 CSV reproduced the 225-record mean change: `0.4329183400186667`.
- All four PNG figures were visually inspected. Corresponding vector PDFs were exported.
- New package files pass whitespace checks. Repository-wide `git diff --check` returns 2 because **pre-existing, untouched** `hspist3/wall_x_FFT.py` has trailing whitespace at lines 1588, 1713 and 2015. Do not clean those unrelated edits as part of this package.
- Python used: `/opt/homebrew/bin/python3`, with NumPy 2.3.3 and Matplotlib 3.10.6 already installed. No dependency was installed. The bundled runtime did not provide Matplotlib.

From the repository root:

```sh
PYTHONPYCACHEPREFIX=/tmp/harddisks-structure-pyc \
/opt/homebrew/bin/python3 -m unittest discover \
  -s hspist3/validation -p test_sound_structure.py -v

PYTHONPYCACHEPREFIX=/tmp/harddisks-structure-pyc \
/opt/homebrew/bin/python3 -m py_compile \
  hspist3/validation/analyze_sound_structure.py \
  hspist3/validation/test_sound_structure.py

MPLCONFIGDIR=/tmp/harddisks-structure-mpl \
PYTHONPYCACHEPREFIX=/tmp/harddisks-structure-pyc \
/opt/homebrew/bin/python3 hspist3/validation/analyze_sound_structure.py \
  --root hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN \
  --out hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/analysis_paper1_20260908 \
  --refresh
```

Default behavior refuses an existing output directory. `--refresh` is only for this package's known generated outputs and refuses extra unknown files. Keep new handoff documents outside that output directory. No simulator is called by this command.

## Next bounded work for Claude Code

Do not repeat a broad audit or restart this package. Read this handoff and its output README, then continue the existing-data sound-speed work:

1. **Resolve the 504 dilute overdue-wall warnings.** Inspect their existing leaf logs and the relevant counter semantics/historical implementation. Report whether these are actual missed events, initialization bookkeeping or unresolved. Do not waive the strict rule or launch 504 reruns without evidence.
2. **Resolve temperature provenance.** The inspected structural CSVs and representative wall traces do not record measured temperature at release. Find whether another durable record or the actual historical normalization protocol establishes kBT/m per trajectory. Do not manufacture T=1 or reconstruct a measured temperature from nominal command flags alone.
3. **Build the trajectory-level fit-input manifest.** Join frequencies/FFT quality, health exclusions, source seed and structural flags. Remove peak-on-search-boundary results from eligible fits. Preserve originals and record unmatched joins instead of guessing. Missing ψ6 alone in dilute gas is not an automatic exclusion from a sound-speed fit.
4. **Recompute eligible sound-speed fits and plots where justified.** Named EOS and documented valid range; fixed-aspect family A finite-size fit with uncertainty/model sensitivity; per-mass residuals and heavy-divider subset. Flag structural-change cases and unresolved temperature cases explicitly. Do not force agreement with the earlier expected values or claim a bulk limit from only three sizes.
5. Report the exact remaining gaps. If temperature or initialization tests genuinely require new simulations, specify the smallest necessary test before launching a new campaign. Finish all unaffected analysis first.

For Paper 2, the next foundational task remains reconciling piston work, gas kinetic energy, moving-wall kinetic energy, spring energy and heat with explicit signs on the existing reference runs. The current package did not implement or validate that ledger. Driven ring-down is not automatically equilibrium data; thermal-wall and information-theory tasks are not prerequisites for every purely mechanical energy-transfer claim.

No new large campaigns, core physics changes, commits, NeurAIpil work or other integrations are part of the suggested next analysis package.
