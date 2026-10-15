# Command and Plot Commands

## Simulation Command

```sh
/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_m3_261009/bin_stageA/00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 --particles=400 --particles-boxes=200,200 --height=20 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=20.000000 --wall-masses=300 --repeats=1 --seed=20261091 --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=200 --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=60 --speed-sound-run-dir=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_m3_261009/evidence_f42befb/acc5/gen2/cell/m_300/r6 --speed-sound-exact-seed=3334348033 --engine=gen2
```

## Plot Command

Run from the repo root:

```sh
python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --dir hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_09_10_26_12_09_09 \
  --write-final \
  --roman-ref \
  --theory-cs simple
```

## Combine Newest Two Eta Runs

Use this after running separate dense/high-eta and low-eta simulation commands:

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

## One-Command High/Mid/Low Workflow

This runs three eta-specific simulation batches and then creates one combined plot:

```sh
python3 hspist3/run_speed_of_sound_eta_split_workflow.py
```

Output is written to a timestamped archive folder:

```text
hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/simulation_eta_split_<timestamp>/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
```

## Main Output

```text
hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_09_10_26_12_09_09/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
```
