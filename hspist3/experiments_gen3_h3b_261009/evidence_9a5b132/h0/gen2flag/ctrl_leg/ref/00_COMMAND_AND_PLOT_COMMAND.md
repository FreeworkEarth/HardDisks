# Command and Plot Commands

## Simulation Command

```sh
/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/gen3/bin_7b_O3 --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=10.0000 --wall-masses=50 --repeats=1 --seed=20261013 --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=25 --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=26 --speed-sound-run-dir=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out/ctrl_leg/ref --speed-sound-exact-seed=57831576 --legacy-resched
```

## Plot Command

Run from the repo root:

```sh
python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --dir hspist3//private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out/ctrl_leg/ref \
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
hspist3//private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out/ctrl_leg/ref/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
```
