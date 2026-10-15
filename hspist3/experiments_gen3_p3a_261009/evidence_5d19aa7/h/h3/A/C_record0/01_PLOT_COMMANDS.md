# Plot commands

Run from the repo root:

```sh
python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --dir hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_p3a_261009/evidence_5d19aa7/h/h3/A/C_record0 \
  --write-final \
  --roman-ref \
  --theory-cs simple
```

Combine the newest two separated eta simulation folders:

```sh
python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --separated-eta-simulations=true \
  --latest-count 2 \
  --write-final \
  --roman-ref \
  --theory-cs simple
```

Combined output:

```text
hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/combined_latest_2_analysis/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
```

Run the full high/mid/low eta split workflow and plot automatically:

```sh
python3 hspist3/run_speed_of_sound_eta_split_workflow.py
```

Main output:

```text
hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_p3a_261009/evidence_5d19aa7/h/h3/A/C_record0/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_p3a_261009/evidence_5d19aa7/h/h3/A/C_record0/analysis_by_eta/final_plots/combined_speed_of_sound_summary.csv
hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_p3a_261009/evidence_5d19aa7/h/h3/A/C_record0/analysis_by_eta/final_plots/speed_of_sound_freq_vs_kterm_combined.pdf
```

Notes:
- The analyzer stages each `L0`/eta group separately and uses eta-specific FFT windows.
- This is preferred over one global FFT window for all eta values.
- If the raw CSV `eta` column is historically wrong, the analyzer recomputes eta from `N`, `r`, `L0`, and `H`.
