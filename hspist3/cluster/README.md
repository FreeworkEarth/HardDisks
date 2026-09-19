# Running the speed-of-sound campaigns on KOA (UH Mānoa)

Written 2026-09-16. Nothing in here changes the physics: same sources, same flags, same
`--edmd-acc=0`, same estimator. What changes is where the runs execute and how they are gated.

## Why this works at all

Every trajectory is one core, one process, its own exact seed and its own output directory. There is
no coupling between runs, so the cluster shape is a **SLURM array**, one task per trajectory — not one
long multi-day job, which walltime limits punish and a single failure destroys.

The simulation never opens a display: `initSDL()` and `TTF_Init()` are behind `if (!cli_headless)`
and the render functions return early when headless. SDL2/GLEW are needed only to link.

## The one rule that matters: cluster results are a separate campaign

KOA is x86, the laptop is ARM. `-O3 -march=native`, FMA contraction and library differences change
the last bit; hard-disk dynamics is chaotic, so two architectures diverge into different trajectories
however correct both are. Therefore:

1. **Never mix laptop and cluster trajectories inside one cell or one figure.** Each campaign root is
   tagged by where it ran.
2. The cluster build is compiled with a **fixed, portable ISA** (`-O2 -march=x86-64-v3`, the `koa`
   target) so that all KOA nodes agree with each other, and `cluster/BUILD_KOA.txt` records compiler,
   flags, git commit and binary hash.
3. The cluster build is accepted only after the **mirror gate** below passes.

## Acceptance gate (mirror campaign)

Three cells that already exist on the laptop — (η = 0.30, N = 400), (η = 0.50, N = 1600),
(η = 0.65, N = 900) — rerun on KOA with the same masses, the same 10 seeds and the same 37.5-period
records. Accept when all four hold:

| check | criterion |
|---|---|
| health contract | every counter zero, exactly as on the laptop |
| per-mass frequency | agrees within the combined standard error, no \|t\| > 3 over the 15 (cell, mass) pairs |
| c_s per cell | agrees within the 1σ mass scatter |
| seeding invariants | `KE_left = KE_right = N_s` exactly, total momentum ≈ 0, ψ₆ distributions overlap |

Byte identity is **not** a criterion and must not be claimed.

## Procedure

```sh
# on the laptop: build the manifest (geometry and seeds come from validation/tests_20260913.py)
python3 cluster/make_manifest.py mirror     manifests/mirror.tsv
python3 cluster/make_manifest.py a2_long200 manifests/a2_long200.tsv
python3 cluster/make_manifest.py n900_sweep manifests/n900_sweep.tsv

# copy source + manifests to KOA
rsync -av --exclude 'experiments_*' --exclude '*.o' hspist3/ koa:~/harddisks/hspist3/

# on KOA
cd ~/harddisks/hspist3
bash cluster/build_koa.sh
mkdir -p logs manifests
sbatch --array=1-150%150 cluster/run_array.sbatch manifests/mirror.tsv

# when the gate passes, the science campaigns
sbatch --array=1-1250%200 cluster/run_array.sbatch manifests/a2_long200.tsv
sbatch --array=1-600%200  cluster/run_array.sbatch manifests/n900_sweep.tsv

# back on the laptop
bash cluster/fetch_results.sh koa ~/harddisks/runs
```

## Campaigns in the manifest generator

| campaign | what | runs |
|---|---|---|
| `mirror` | the acceptance gate, 3 existing cells, 37.5 periods | 150 |
| `a2_long200` | full A2 ladder, 200-period records (5 η × 5 N × 5 masses × 10 seeds) | 1250 |
| `n900_sweep` | second c_s(η) curve at N = 900, 12 densities, 200 periods | 600 |
| `edmd_acc` | the mirror cells on both backends, for the accelerated-backend test | 300 |

Cost on the laptop for comparison: `a2_long200` alone is ≈ 1400 core-hours ≈ 140 h at 10 slots; at
200 KOA cores it is ≈ 7 h.

## Still to fill in from KOA

`run_array.sbatch` carries `--partition=PARTITION_TBD` and `--time=04:00:00`. Replace both from:

```sh
sinfo -o "%P %l %D %c %m" | head -20
module avail 2>&1 | grep -iE "sdl|glew|mesa|gcc" | head
sacctmgr show assoc where user=$USER format=Account,Partition,MaxJobs,GrpTRES | head
lscpu | grep -E "Model name|Flags" | head -2      # confirm x86-64-v3 (avx2, bmi2) is supported
df -h $HOME "$SCRATCH" 2>/dev/null
```

The longest single run in `a2_long200` (η = 0.10, N = 2500, M = 2000, 200 periods) took ≈ 11 h on one
laptop core, so its array task needs a walltime above that or that cell must be split by seed.
