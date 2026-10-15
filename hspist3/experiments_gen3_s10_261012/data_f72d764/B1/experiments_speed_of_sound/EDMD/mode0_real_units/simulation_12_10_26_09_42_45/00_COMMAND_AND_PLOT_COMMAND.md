# Command and Plot Commands

## Trace columns under --engine=gen3

- `Left_Count`, `Right_Count`: computed from a synchronised particle copy at the driver's validator steps (every --validator-every steps, default 60 = 1 sigma-time), at the first recorded row and at the last step; in the rows between they carry the last computed values. No disk can cross the divider, and the validator checks every disk's compartment against its initial one at each of its steps.
- `psi6_t_L0_<L>_wallmassfactor_<M>_run<r>.csv`: psi6(t) every 2.53271 sigma-time (--psi6-every) through hold and record.
- The run log holds each trajectory's run record `[EDMD3-HEALTH] ...: clean=...` (every health counter, the engine and the build line).

## Simulation Command

```sh
/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_d18_261010/bin_f72d764/00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 --particles=100 --particles-boxes=50,50 --height=10 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=13.089969389957473 --wall-masses=50 --repeats=1 --seed=7700001 --wall-hold-steps=60000 --fixed-dt=0.4 --target-oscillations=200 --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=46 --speed-sound-run-dir=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_s10_261012/data_f72d764/B1 --speed-sound-exact-seed=7700001 --engine=gen3 --gen3-seeding=lattice --gen3-exact-box --gen3-psi6-compartments --psi6-every=2.5327085535440355 --gen3-virial-blocks=50 --gen3-snapshots=60.654171070880714
```

## Plot Command

Run from the repo root:

```sh
python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --dir hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_12_10_26_09_42_45 \
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
hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_12_10_26_09_42_45/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
```
