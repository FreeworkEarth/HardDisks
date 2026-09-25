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
2. The cluster build is compiled with a **fixed, portable ISA** (`-O2 -march=x86-64-v2`, the `koa`
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

## KOA facts (as of the account briefing, 2026-10-01)

| | |
|---|---|
| `sandbox` | 4 h walltime — interactive tests and the first-hour checks only |
| `shared` | **3 d** walltime, up to **93 cores per node** — the science partition |
| `kill-shared` | same but **preemptable**: jobs can be killed at any moment |
| home | **50 GiB, NOT backed up** — source and scripts only, never campaign output |
| `~/koa_scratch` | no quota, **purged after 90 days** — all campaign output lives here |
| transfer | DTN **`koa-dtn.its.hawaii.edu`** — use it for rsync, not the login node |
| modules | **Lmod** (`module avail`, `module load`) |
| CPUs | **heterogeneous, Ivy Bridge included** — see the ISA warning below |

**Preemption changes how you submit.** `kill-shared` is free capacity but a task can die mid-run, so
it is only usable because `run_array.sbatch` is **restart-safe**: it skips a trajectory whose trace
already reaches its planned duration, and `--requeue` is set. Submit long cells to `shared` and use
`kill-shared` for the wide, cheap, easily-redone arrays.

**The 90-day purge is a deadline, not a detail.** `cluster/fetch_results.sh` must be run before any
campaign output ages out; nothing on KOA is backed up anywhere.

## Three ways to lose a week

**1. ISA: use `-march=x86-64-v2`, not v3.** KOA's shared partition includes **Ivy Bridge** nodes.
`x86-64-v3` requires AVX2/FMA (Haswell and later), so a v3 binary dies with **SIGILL the moment
Slurm lands a task on an older node** — intermittently, on some array tasks only, which is the worst
possible way to discover it. `x86-64-v2` (SSE4.2/POPCNT, Nehalem and later) runs on every node in
the partition. The `koa` Makefile target was corrected from v3 to v2 on 2026-10-01. Do not raise it
without pinning `--constraint` to a node feature, and if you ever do, the mirror gate must be re-run.

**2. `make` builds AddressSanitizer, not science.** The default target is `debug`
(`-g -fsanitize=address -DDEBUG`), ≈ 5× slower. Science is `make release` on the laptop and
**`make koa`** on the cluster. A campaign on the debug build is not wrong, only slow — but it is not
byte-comparable with a release build and cannot be pooled with one.

**3. Trace size.** The energy-transfer trace writes one row per step and ignores `--output-dt`
(≈ 368 bytes/row). Always pass **`--trace-every=N`** for Level 4 cells — see the table in
`level4_massladder.sbatch`. Undecimated, the M_d = 200 cell below would be **14.5 GB per seed**.

## Provenance every node must record

Because acceptance across nodes is **statistical, never byte-identical**, each run has to say where
it ran. Every task appends to its `00_COMMAND.md`:

```sh
{ echo; echo "## node provenance";
  echo "- hostname: $(hostname)";
  echo "- cpu: $(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ *//')";
  echo "- slurm: job ${SLURM_JOB_ID:-none} task ${SLURM_ARRAY_TASK_ID:-none} on ${SLURM_JOB_PARTITION:-none}";
  echo "- binary: $(sha1sum "$BIN" | cut -d' ' -f1)";
} >> "$outdir/00_COMMAND.md"
```

If a cell's results ever split into two clusters, the CPU model is the first thing to check.

## Commands to confirm the above on first login

```sh
sinfo -o "%P %l %D %c %m" | head -20
sacctmgr show assoc where user=$USER format=Account,Partition,MaxJobs,GrpTRES | head
module avail 2>&1 | grep -iE "sdl|glew|mesa|gcc" | head
lscpu | grep -E "Model name|Flags" | head -2     # v2 needs sse4_2 + popcnt, NOT avx2
df -h $HOME ~/koa_scratch
```

The longest single run in `a2_long200` (η = 0.10, N = 2500, M = 2000, 200 periods) took ≈ 11 h on one
laptop core. That fits `shared`'s 3 d limit but not `sandbox`'s 4 h.

See `RUNBOOK_first_hour.md` for the sequence to run on day one.
