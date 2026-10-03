# KOA runsheet: first build, smoke test, method-A pilot (written 2026-10-02)

Every command below is typed by Chris. Each step has a stop condition; if the output does not match, stop and paste
the last ~30 lines to CC. Nothing here starts a campaign: the arrays come only after the pilot's gate and the go.

## Rules (read once)

- **Login node `login-0102`: light work only.** You may `git`, `ls`, `sbatch`, `squeue`, `tail`. Never compile or run
  the simulation there. Home is mounted `noexec` on the login node, so `./00ALLINONE` would not even start there.
- **`srun` / `sbatch` only from a `[charing@login-0102 ...]$` prompt, never from a `[charing@cn-...]$` prompt.** Typed
  inside a compute-node shell, `srun` hangs ("step creation temporarily disabled").
- **After `srun ... --pty /bin/bash`, wait until the prompt changes to `cn-...`** before pasting the next line. Lines
  pasted while the job is still queued are lost.
- **Data go to scratch** (`/mnt/lustre/koa/scratch/charing`, = `~/koa_scratch`), never to home. Scratch deletes files
  90 days after their last write: copy results back (step 7) well before then.
- **Nothing is overwritten.** Every script refuses an output directory that already exists.
- **To stop a job:** `scancel <jobid>`. To leave a compute-node shell: `exit` (the prompt returns to `login-0102`).
- **Onboarding quiz:** finish it at lamaku.hawaii.edu within 30 days of the account e-mail, or access is suspended.

## 1. Get the code (inside a sandbox session, ~10 s)

Do this once, after you pushed from the Mac. It fetches only the sources and scripts (about 11 MB), never binaries
or data. **The login node has no git** (found 2026-10-02: `git: command not found`); the compute nodes have
/usr/bin/git 2.52.0. So open the sandbox session of step 2 first, wait for the `cn-...` prompt, and clone there:

```sh
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

then, at the `cn-...` prompt:

```sh
cd ~
git clone --depth 1 --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks
cd harddisks
git sparse-checkout set --no-cone "/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"
git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft
git -C hspist3/kissfft -c advice.detachedHead=false checkout febd4caeed32e33ad8b2e0bb5ea77542c40f18ec
git log --oneline -1
cd hspist3 && mkdir -p logs
ls Makefile 00ALLINONE.c experiment_validation.c plot_speed_of_sound_edmd.py validation/tests_20260913.py cluster/koa_smoketest.sh
```

- `--depth 1 --sparse`: only the newest commit, only the source folders.
- kissfft is the FFT library. It is a separate public repository, pinned to the exact commit the Mac uses.
- **Expected:** `git log --oneline -1` prints the same hash and message as `git log --oneline -1` on the Mac after
  your push. That hash is what the binary will report as its `build_git`.
- **The `ls` line must list all six files without an error.** If `plot_speed_of_sound_edmd.py` is missing, it was not committed
  and pushed from the Mac (it is the analysis module with the equation of state; untracked until Chris decides): stop here.
  The smoke test would stop on it anyway, within seconds (its preflight imports the analysis modules before building).
- **If the clone fails** (no network from the login node), stop and tell CC. Do not copy files by hand.
- **Later updates:** `cd ~/harddisks && git pull`, also inside a sandbox session (after a push from the Mac). Then rebuild in step 2.

## 2. Build once, interactively (compute node, ~5 min)

In the same `cn-...` session as step 1 (if you left it, start it again with
`srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash` and wait for the `cn-...` prompt):

```sh
cd ~/harddisks/hspist3
bash cluster/build_koa.sh
```

- `srun ... --pty`: a shell on a compute node, 2 cores, at most 1 hour.
- `build_koa.sh`:
  - loads the compiler (module `compiler/GCC/14.3.0`) and the libraries from `~/envs/hd` (`cluster/koa_env.sh`);
  - compiles with `-O2 -march=x86-64-v2 -ffp-contract=off`;
  - checks that every library resolves (`ldd`) and that the binary's hash is the clean checkout's HEAD;
  - writes `logs/BUILD_KOA_<jobid>.txt`;
  - (added 2026-10-02, Task U4) records the binary's sha256 in `logs/BUILD_KOA_LAST.hash`, which every array checks.
- **Expected:**
  - the first line is `env: gcc (GCC) 14.3.0 | Python 3.12.14 at /home/charing/envs/hd/bin/python3 | git version ...`;
  - a line `== make compiler: gcc -> /.../gcc -> gcc (GCC) 14.3.0` (the compiler make really uses; if it says 11.5, stop);
  - then the compiler output (warnings are fine);
  - a line `recorded       logs/BUILD_KOA_LAST.hash: <sha256>  00ALLINONE` (from 2026-10-02 on);
  - the last line is `BUILD OK`.
- **Any `STOP:` line:** stop, `exit`, and paste the output.

## 3. Check the version (same compute-node shell)

```sh
./00ALLINONE --version
exit
```

- **Expected:**
  - `00ALLINONE  git <the hash from step 1>  target koa`;
  - `CFLAGS: -O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off`.
- No `-dirty`, no `unknown`.
- After `exit` the prompt is `login-0102` again.

## 4. Submit the smoke test (login node; ~15 min of run time plus queue wait)

```sh
cd ~/harddisks/hspist3
sbatch cluster/koa_smoketest.sh
squeue -u charing
```

- `sbatch` prints `Submitted batch job <jobid>`. Note the number.
- The job uses `sandbox`, 9 cores, at most 1 hour. It does four things:
  1. builds again, inside the job, so the log carries the provenance;
  2. runs the determinism test, the same seed twice as two separate steps;
  3. runs the π/8 pilot: nine masses, one seed each, 200 oscillations;
  4. checks the gates.
- **Watch it:**
  - `squeue -u charing`: `PD` = waiting, `R` = running, nothing listed = finished.
  - `tail -f logs/conf-smoke_<jobid>.out` shows the log live. **Ctrl+C stops only the watching, not the job.**
- **Then the cross-node check** (needs the binary from the smoke test; 2 cores, about 1 minute):

  ```sh
  sbatch cluster/koa_crossnode_det.sh
  ```

  Its log is `logs/det-xnode_<jobid>.out`.
- **Paste back to CC:** `tail -40 logs/conf-smoke_<jobid>.out` and `tail -8 logs/det-xnode_<jobid>.out`. Only the
  tails: the chat cuts pastes above 50 000 characters.

## 5. What PASS looks like

- **Smoke test log:**
  - `wall_x_positions_..._run0.csv: IDENTICAL` and `speed_of_sound_psi6.csv: IDENTICAL`, then
    `determinism self-test (... same node): IDENTICAL`.
  - `pilot cell: eta (trace) = [0.392699], L_0 (trace) = [10.0], ... L_eff = L_0 - 2r - t/2 = 8.975000`.
  - `health lines = 0`.
  - Gates: `eta: PASS`, `L_0: PASS`, `L_eff: PASS`, `health: PASS`, and
    `c_s within 0.06903 (2 sigma_diff) of 3.85886: PASS`.
  - The last line is `SMOKE TEST PASSED -- the method-A pilot may be submitted (after the go)`.
- **Cross-node log:** `determinism self-test (... run A on cn-X, run B on cn-Y, different nodes): IDENTICAL`.
- **Gate numbers** (pre-registered, 261012 § 1.10 and § 1.10.1):
  - determinism IDENTICAL;
  - η, L_0, H and L_eff equal to the Mac to 1e-6;
  - |c_s(KOA) − 3.85886| ≤ 0.06903.

  The Mac and KOA runs are not expected to be byte-identical (different CPU and maths library), so the c_s gate is
  statistical.
- **Any FAIL:** stop. No further submission. Paste the tails.

## 6. Method-A pilot (ONLY after CC/the plan author give the go on step 5)

```sh
cd ~/harddisks/hspist3
sbatch --array=1-1 cluster/confinement_20261013/conf_A_pilot.sbatch
```

- It runs on `sandbox`: 16 cores, at most 1 hour, 20 held-divider trajectories at the π/8 anchor. It does not build;
  it checks that `./00ALLINONE` is the clean koa build of HEAD.
- Log: `logs/conf-A_pilot_<jobid>_1.out`. **Expected last line:** `cell pilot_epi8_H_H10_L10 done; failures: 0`.

## 7. Copy results back (on the Mac, from the repo root)

Smoke-test pilot only (7.4 MB on the Mac run):

```sh
cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
rsync -av charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/ hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/
```

The whole confinement campaign (summaries, plus every file of the two pilot cells), once it exists:

```sh
bash hspist3/cluster/confinement_20261013/fetch_confinement.sh
```

- Both go through the data transfer node `koa-dtn` (password + Duo again).
- `rsync` only copies; it deletes nothing on either side.

## 8. Round plan (written 2026-10-02; nothing here is submitted until the pilot's gate 4 and the go)

**Rules for this step**: the rules at the top of this runsheet apply unchanged. There is one correction, dated 2026-10-02:
- The rule "Every script refuses an output directory that already exists" holds for the smoke test, the cross-node check and `confinement_pilot.py det1`.
- It does **not** hold for the five `conf_*.sbatch` arrays. They **resume**: `conf_worker.sh` skips a trajectory whose output already exists (B: `..._run<r>.csv`; A: a non-empty `red_<seed>.csv`) and runs the rest. A resubmitted array therefore continues; it does not start over and does not overwrite.
- All data go under `/mnt/lustre/koa/scratch/charing/harddisks/hspist3/` (`HD_DATA`). The Slurm logs go to `~/harddisks/hspist3/logs/` in home.
- Checked statically by `python3 hspist3/cluster/round_plan_261002.py` (Mac).

**Pull and rebuild rules (2026-10-02, Task U4; this replaces the earlier "do not `git pull` before Round 1").** The arrays no longer compare the binary with `git rev-parse HEAD`. They compare it with `logs/BUILD_KOA_LAST.hash`, which `cluster/build_koa.sh` writes. A later `git pull` therefore cannot stop a running campaign.
1. **After every `git pull`, rebuild (step 2) before any NEW submission.** Do the pull and the build in the same sandbox session.
2. **Never pull or rebuild while an array still has tasks pending or running** (`squeue -u charing` must be empty). A pending task reads the task files and the binary only when it starts. `conf_worker.sh` refuses to resume a directory that another build wrote (`FAILED build guard`), so such tasks would fail rather than mix builds.
3. **The first submission after this change needs a pull and a rebuild** (step 8e). The 70b2069 build did not write `logs/BUILD_KOA_LAST.hash`, so without a rebuild every array stops at once with `STOP: no logs/BUILD_KOA_LAST.hash`.

**The numbers in 8b–8d (planning seeds, Mac speed) are superseded by 8e** (gate 4 applied, KOA speed measured). Submit the arrays from 8e.

### 8a. Gate 4: the method-A pilot (job 14966594)

1. **On KOA, at the `login-0102` prompt.** The pilot needs no extra reduction: `conf_worker.sh` mode A already calls `reduce_A.py` per seed and writes `red_<seed>.csv`. Check the log's last line, then read three numbers:

   ```sh
   tail -15 logs/conf-A_pilot_14966594_1.out
   sacct -j 14966594 --format=JobID,Elapsed,TotalCPU,State
   lfs quota -h -u charing /mnt/lustre/koa
   ```

   - **Expected:** `cell pilot_epi8_H_H10_L10 done; failures: 0`.
   - `sacct` gives the pilot's elapsed time, which sets the KOA speed factor `--slow` below.
   - `lfs quota` gives the scratch quota, which is OPEN; it decides 8c.
   - All three are read-only and light, so they are fine on the login node.

2. **On the Mac, from the repo root.** Copy the pilot's summaries only: `red_*.csv`, `run_*.log`, `summary_*.csv`. The ~0.4 GB of event logs stay on scratch.

   ```sh
   cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
   rsync -av --prune-empty-dirs --include='*/' --include='red_*.csv' --include='run_*.log' --include='summary_*.csv' --exclude='*' charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks/hspist3/experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10/ hspist3/experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10/
   ```

3. **On the Mac, apply gate 4.** The script prints the verdict and the seeds per position for conf_A_0.39:

   ```sh
   python3 hspist3/cluster/gate4_pilot_261002.py
   ```

   - It first checks that the pre-registered rule, with the planning ε₀ = 0.2710, reproduces the committed task files: 9 of 9 cells, PASS on 2026-10-02 with `--dry`.
   - It then requires 20 seeds, health 0, and a window of 5000 ± 1 % σ-time per seed, read from the event log. The 1 % is this script's tolerance, not a pre-registered number.
   - It measures ε₀ at π/8 and recomputes T per position and seeds per position by the § 1.4 rule.
   - **conf_A_0.10 is not affected**: its ε₀ was measured at η = 0.10005 (Level 3), and § 1.8 lets the pilot change only (A) at π/8.

### 8b. Round 1: conf_B_0.10, conf_B_0.39, conf_A_0.10 (shared; all three may run together)

At most 64 cores at once, summed over the arrays that run together. Why 64:
- `shared` is used by every KOA user, and 64 cores is about three of its nodes;
- a systematic error then costs at most 64 cores × the time until Chris notices;
- the whole campaign is ~96 core-hours, so more cores would save little wall time.

KOA's own per-user limits are not read yet (OPEN; `sacctmgr -n show assoc user=charing format=account,partition,maxjobs,maxsubmit,grptres%30,maxtres%30` shows them).

Printed by `python3 hspist3/cluster/round_plan_261002.py` (cost model of the pre-registration, Mac speed; rerun it with `--slow F` once the pilot gives the KOA speed):

| round | array | partition | tasks | cores/task | traj. | core-h | longest cell (h) | --time (h) | throttle | cores at once | wall (h) | scratch GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Round 1 | conf_B_0.10 | shared | 10 | 8 | 2250 | 21.8 | 0.74 | 2 | %2 | 16 | 1.50 | 1.8 |
| Round 1 | conf_B_0.39 | shared | 9 | 8 | 2025 | 17.7 | 0.52 | 1 | %2 | 16 | 1.23 | 1.6 |
| Round 1 | conf_A_0.10 | shared | 10 | 16 | 7030 | 4.5 | 0.08 | 1 | %2 | 32 | 0.14 | 34.9 |
| Round 2 | conf_A_0.39 | shared | 9 | 16 | 10465 | 51.8 | 1.57 | 3 | %4 | 64 | 1.57 | 235.1 |


Wall time is from the cost model at Mac speed. "Longest cell" must stay below `--time`: the margins are ×2.7 (B_0.10), ×1.9 (B_0.39), ×12 (A_0.10) and ×1.9 (A_0.39). If the pilot shows KOA running more than ~1.8× slower than the Mac, `--time` of B_0.39 and A_0.39 must be raised before submission.

**On KOA, at the `login-0102` prompt, from `~/harddisks/hspist3`:**

```sh
cd ~/harddisks/hspist3
mkdir -p logs
sbatch --array=1-10%2 cluster/confinement_20261013/conf_B_0.10.sbatch
sbatch --array=1-9%2 cluster/confinement_20261013/conf_B_0.39.sbatch
sbatch --array=1-10%2 cluster/confinement_20261013/conf_A_0.10.sbatch
```

- `%2` means at most two cells of that array run at a time: 2 × 8 + 2 × 8 + 2 × 16 = 64 cores.
- Each log ends with `cell <id> done; failures: 0`.
- Check with `squeue -u charing`; stop an array with `scancel <jobid>`.

### 8c. Scratch: decision before Round 2 (OPEN)

The estimate comes from the E5 event log (DATA: 20424 rows over a wall length of 80 and 400 σ-time) and the contact theorem (the impact rate per unit length is $nZ/\sqrt{2\pi}$).
- The method-A event logs need about 35 GiB for conf_A_0.10 and **about 235 GiB for conf_A_0.39**. Method B needs under 2 GiB per array.
- If `lfs quota` (8a) shows less room than that, the event logs must be compressed or removed after each seed's reduction. That is a decision for Chris and the plan author, because it changes what stays regenerable on scratch.

### 8d. Round 2: conf_A_0.39 (only after gate 4 PASS and the go)

1. On the Mac, the conf_A_0.39 task files are regenerated with the pilot's ε₀; that is a separate CC task. It is committed, and Chris pushes.
2. On KOA, inside a sandbox session (step 2), run `git pull`, then rebuild with `bash cluster/build_koa.sh` (step 3).
3. Then, from `login-0102`:

   ```sh
   cd ~/harddisks/hspist3
   sbatch --array=1-9%4 cluster/confinement_20261013/conf_A_0.39.sbatch
   ```

   `%4` × 16 cores = 64. The row above uses the planning seeds. After gate 4, `python3 hspist3/cluster/round_plan_261002.py --eps0-ratio R` prints the updated row, with R = ε₀(pilot)/ε₀(plan) as gate 4 prints it.

### 8e. Launch plan after gate 4 (written 2026-10-02; Round 1 needs the go, Round 2 a second go)

**Gate 4 PASSED** (261012 § 1.12). ε₀(π/8) = 0.0820 ± 0.0106, against 0.2710 planned. The conf_A_0.39 seeds per position were recomputed by the pre-registered rule, and the generator rewrote the task files.

Printed by `python3 hspist3/cluster/round_plan_261002.py`:

### conf_A_0.39 task files: planning seeds (git 70b2069) vs gate-4 seeds (working tree)

| cell | seeds/position at 70b2069 | seeds/position now | lines now = first lines of each position at 70b2069 | seeds now |
|---|---|---|---|---|
| epi8_H_H5_L10 | 219 | 20 | yes | 9700..9719 |
| epi8_H_H10_L10 | 194 | 18 | yes | 9700..9717 |
| epi8_H_H20_L10 | 219 | 20 | yes | 9700..9719 |
| epi8_H_H40_L10 | 437 | 40 | yes | 9700..9739 |
| epi8_L_H10_L5 | 219 | 20 | yes | 9700..9719 |
| epi8_L_H10_L20 | 219 | 20 | yes | 9700..9719 |
| epi8_aspect_H7.08333_L14.125 | 218 | 20 | yes | 9700..9719 |
| epi8_aspect_H5_L20 | 194 | 18 | yes | 9700..9717 |
| epi8_aspect_H3.54167_L28.2917 | 174 | 16 | yes | 9700..9715 |

KOA speed (measured, pilot 14966594): 371.764 CPU-s / 20 = 18.59 CPU-s per trajectory; Mac cost model 10.30 -> factor 1.804
times below: cost model x 1.804 (--slow); trajectories per cell from the task files

| round | array | partition | tasks | cores/task | traj. | core-h (KOA) | longest cell (h) | --time (h) | rule 2 x longest (h) | throttle | cores at once | wall (h) | scratch GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Round 1 | conf_B_0.10 | shared | 10 | 8 | 2250 | 39.4 | 1.33 | 2.75 | 2.75 | %2 | 16 | 2.71 | 1.8 |
| Round 1 | conf_B_0.39 | shared | 9 | 8 | 2025 | 32.0 | 0.94 | 2 | 2 | %2 | 16 | 2.22 | 1.6 |
| Round 1 | conf_A_0.10 | shared | 10 | 16 | 7030 | 8.2 | 0.15 | 0.5 | 0.5 | %2 | 32 | 0.26 | 34.9 |
| Round 2 | conf_A_0.39 | shared | 9 | 16 | 960 | 8.6 | 0.27 | 0.75 | 0.75 | %4 | 64 | 0.27 | 21.6 |

**1. On the Mac:** push (`bash _commit_scripts/commit_20261007.sh`).

**2. On KOA, from a `login-0102` prompt:** check that nothing of yours is queued, then open a sandbox session.

```sh
squeue -u charing
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

`squeue` must list no jobs. Wait for the `cn-...` prompt.

**3. At the `cn-...` prompt:** pull, rebuild, check.

```sh
cd ~/harddisks/hspist3
git pull
bash cluster/build_koa.sh
cat logs/BUILD_KOA_LAST.hash
wc -l cluster/confinement_20261013/tasks_A_epi8_H_H10_L10.txt
exit
```

**Expected:**
- `version        00ALLINONE  git <the pushed commit>  target koa`, then `recorded       logs/BUILD_KOA_LAST.hash: ...`, then `BUILD OK`;
- `90 cluster/confinement_20261013/tasks_A_epi8_H_H10_L10.txt` (18 seeds × 5 positions).

Any `STOP:` line: stop and paste the output.

**4. Round 1, back at the `login-0102` prompt (after the go):**

```sh
cd ~/harddisks/hspist3
mkdir -p logs
sbatch --array=1-10%2 cluster/confinement_20261013/conf_B_0.10.sbatch
sbatch --array=1-9%2 cluster/confinement_20261013/conf_B_0.39.sbatch
sbatch --array=1-10%2 cluster/confinement_20261013/conf_A_0.10.sbatch
```

- 64 cores at once; about 2.7 h of wall time at the measured KOA speed.
- Every log starts with the binary's `--version` line and ends with `cell <id> done; failures: 0`.
- A `FAILED build guard` line means that a directory was written by another build. Stop and paste it.

**5. Round 2 (after its go; the binary is the same, so no pull and no rebuild in between):**

```sh
sbatch --array=1-9%4 cluster/confinement_20261013/conf_A_0.39.sbatch
```

- 64 cores at once; about 0.3 h.
- **Cross-check, free of charge:** the anchor cell `epi8_H_H10_L10` reruns the pilot's seeds 9700–9703 at all five positions. Its `red_970[0-3].csv` must equal the pilot's byte for byte (determinism). Their `red_*.csv` files come back with the summaries.
