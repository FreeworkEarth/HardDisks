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
4. **Exception (2026-10-03, Task W4): a pull may be followed by resubmission WITHOUT a rebuild if its diff touches no build input** (`*.c`, `*.h`, `Makefile`, `edmd_core/`, `kissfft`). This is safe for two reasons:
   - the arrays verify the binary by sha256 against `logs/BUILD_KOA_LAST.hash`, not by the commit;
   - the `.build_git` records hold the binary's own `--version` line, which a pull does not change.

   Check with `git diff --name-only <build commit>..HEAD`, where the build commit is the `git` hash in `./00ALLINONE --version`. Otherwise rebuild.
5. **Do not rebuild while any cell is incomplete.** A rebuild at a new commit changes the binary's `--version` line. Every directory already written by the old build is then REFUSED by the build guard (by design: no cell mixes builds), so its missing seeds could no longer be finished.

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

**Gate 4 PASSED** (261012 § 1.12). ε₀(π/8) = 0.0820 ± 0.0106, against 0.2710 planned. By **amendment C3** (261012 § 1.9, 2026-10-02), the conf_A_0.39 seeds per position come from the upper bound 0.0926 by the unchanged rule. The generator rewrote the task files.

Printed by `python3 hspist3/cluster/round_plan_261002.py` (re-run 2026-10-02 after amendment C3):

### conf_A_0.39 task files: seeds per position, plan (git 70b2069), gate 4 (git 303280d) -> now (working tree)

| cell | plan (70b2069) | gate 4 (303280d) | now | lines | nested in each earlier file | seeds now |
|---|---|---|---|---|---|---|
| epi8_H_H5_L10 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_H_H10_L10 | 194 | 18 | 23 | 115 | yes | 9700..9722 |
| epi8_H_H20_L10 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_H_H40_L10 | 437 | 40 | 51 | 255 | yes | 9700..9750 |
| epi8_L_H10_L5 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_L_H10_L20 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_aspect_H7.08333_L14.125 | 218 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_aspect_H5_L20 | 194 | 18 | 23 | 115 | yes | 9700..9722 |
| epi8_aspect_H3.54167_L28.2917 | 174 | 16 | 21 | 105 | yes | 9700..9720 |

KOA speed (measured, pilot 14966594): 371.764 CPU-s / 20 = 18.59 CPU-s per trajectory; Mac cost model 10.30 -> factor 1.804
times below: cost model x 1.804 (--slow); trajectories per cell from the task files

| round | array | partition | tasks | cores/task | traj. | core-h (KOA) | longest cell (h) | --time (h) | rule 2 x longest (h) | throttle | cores at once | wall (h) | scratch GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Round 1 | conf_B_0.10 | shared | 10 | 8 | 2250 | 39.4 | 1.33 | 2.75 | 2.75 | %2 | 16 | 2.71 | 1.8 |
| Round 1 | conf_B_0.39 | shared | 9 | 8 | 2025 | 32.0 | 0.94 | 2 | 2 | %2 | 16 | 2.22 | 1.6 |
| Round 1 | conf_A_0.10 | shared | 10 | 16 | 7030 | 8.2 | 0.15 | 0.5 | 0.5 | %2 | 32 | 0.26 | 34.9 |
| Round 2 | conf_A_0.39 | shared | 9 | 16 | 1240 | 11.0 | 0.33 | 0.75 | 0.75 | %4 | 64 | 0.33 | 27.8 |

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
- `115 cluster/confinement_20261013/tasks_A_epi8_H_H10_L10.txt` (23 seeds × 5 positions, amendment C3).

Any `STOP:` line: stop and paste the output.

**Launch rules (2026-10-02, Task V3):**
- **All four arrays are submitted from ONE build:** the one made in step 3.
- **No `git pull` between Round 1 and Round 2.** Round 2 needs nothing newer than that build, and a pull would require a rebuild (rule 1 above).
- **Round 2 goes in only after Round 1's first tasks have run cleanly for about 30 minutes.** Before submitting it, check:
  - `squeue -u charing` shows the Round-1 tasks running (`R`), not pending with an error;
  - `grep -l "STOP\|FAILED" logs/conf-*` prints nothing;
  - `head -3 logs/conf-B_0.10_<jobid>_1.out` shows the gcc 14.3.0 line and the binary's `--version`.

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

**5. Round 2 (after its go, and after Round 1 has run cleanly for about 30 min; the same build, so no pull and no rebuild in between):**

```sh
sbatch --array=1-9%4 cluster/confinement_20261013/conf_A_0.39.sbatch
```

- 64 cores at once while Round 1 is still running, so up to 128 cores in total for a short while; this is about 0.3 h and 11.0 core-h. If that is too much, wait until Round 1 has finished.
- **Determinism gate (261012 § 1.12, V2):** the anchor cell `epi8_H_H10_L10` reruns the pilot's seeds 9700–9703 at all five positions. Its `red_970[0-3].csv` must equal the pilot's byte for byte (determinism). Their `red_*.csv` files come back with the summaries.

### 8f. Round 1 repair (written 2026-10-03; nothing is submitted before the cell check is read and the go is given)

Round 1 ran on build `279282b target koa` (261012 § 1.13). Two things went wrong:
- the build guard's lock never locked, so some seeds were skipped;
- the three H = 40 cells hit their time limit.

The fix touches only scripts. **No build input changed, so there is no rebuild** (rule 4 above).

**1. On the Mac:** push (`bash _commit_scripts/commit_20261007.sh`).

**2. On KOA, from a `login-0102` prompt:** check that nothing is queued, then open a sandbox session.

```sh
squeue -u charing
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

`squeue` must list no jobs. Wait for the `cn-...` prompt.

**3. At the `cn-...` prompt:** pull, confirm that no build input changed, check the cells.

```sh
cd ~/harddisks/hspist3
git pull
git log --oneline -1
git diff --stat 279282b..HEAD
git diff --name-only 279282b..HEAD | grep -E '^hspist3/([^/]+\.(c|h)|Makefile|edmd_core/|kissfft)' || echo "no build input changed"
sha256sum --status -c logs/BUILD_KOA_LAST.hash && echo "binary = recorded build"
bash cluster/check_cells.sh
exit
```

**Expected:**
- `no build input changed`, then `binary = recorded build`;
- one block per array, with one verdict per cell, then a `SUMMARY:` line;
- `conf_A_0.39` and its cells are `NOT STARTED`.

**Paste everything from `git log` to `SUMMARY:`.** If `git diff 279282b..HEAD` says it does not know `279282b`, use `git diff --stat HEAD@{1}..HEAD` instead.

**4. After the go, from `~/harddisks/hspist3` on `login-0102`.**

The lines are printed by `python3 hspist3/cluster/round1_timing_261003.py` on the Mac, with `--time` from the measured Round 1 times:

Set 1 (now); at most 64 cores at once:
    sbatch --array=4 --time=18:30:00 cluster/confinement_20261013/conf_B_0.10.sbatch
    sbatch --array=1-3,5-10%1 cluster/confinement_20261013/conf_B_0.10.sbatch
    sbatch --array=4 --time=1-02:15:00 cluster/confinement_20261013/conf_B_0.39.sbatch
    sbatch --array=1-3,5-9%1 cluster/confinement_20261013/conf_B_0.39.sbatch
    sbatch --array=4 --time=03:00:00 cluster/confinement_20261013/conf_A_0.10.sbatch
    sbatch --array=1-3,5-10%1 cluster/confinement_20261013/conf_A_0.10.sbatch
Set 2 (Round 2, behind the A_0.10 resubmission); at most 48 cores at once:
    sbatch --array=1,2,5,7,8,9%1 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
    sbatch --array=3,6%1 --time=01:15:00 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
    sbatch --array=4 --time=17:30:00 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
Set 2 starts while the two long B H40 tasks of set 1 may still run: 16 + set 2 = 64 cores (cap 64).

- **Set 1** is all six lines, submitted together. Finished cells exit in seconds; only the missing seeds and the three H = 40 cells run.
- **Set 2** has three lines. Write down the two A_0.10 job numbers that `sbatch` prints in set 1: the line `--array=4 … conf_A_0.10` and the line `--array=1-3,5-10%1 … conf_A_0.10`. Put them into the `<…>` placeholders of set 2.
- **No `git pull` and no rebuild** between the sets, or while any of them runs (rule 5).

## 9. A-fixed (261012 § 3; written 2026-10-04; nothing is submitted before gate G1 and the go)

**What it is.** Method A again, but with the divider **held for the whole record**. Same cells, positions and seeds; new output directory `experiments_energy_transfer/paper1_confinement_Afix_261004/`; the method-A directories are never written. The binary is unchanged.

**Rules.** Those of step 8 apply:
- one build for everything;
- rule 4: a pull that touches no build input may be followed by a submission without a rebuild;
- rule 5: no rebuild while any cell is incomplete.

**1. On the Mac:** push (`bash _commit_scripts/commit_20261007.sh`).

**2. On KOA, from a `login-0102` prompt:** check that nothing is queued, then open a sandbox session.

```sh
squeue -u charing
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

`squeue` must list no jobs.

**3. At the `cn-...` prompt:** pull, confirm that no build input changed, confirm that the binary is the recorded build.

```sh
cd ~/harddisks/hspist3
git pull
git log --oneline -1
git diff --name-only 279282b..HEAD | grep -E '^hspist3/([^/]+\.(c|h)|Makefile|edmd_core/|kissfft)' || echo "no build input changed"
sha256sum --status -c logs/BUILD_KOA_LAST.hash && echo "binary = recorded build"
exit
```

**Expected:** `no build input changed`, then `binary = recorded build`.

**4. The A-fixed pilot (gate G1), from `login-0102`:**

```sh
cd ~/harddisks/hspist3
mkdir -p logs
sbatch --array=1-1 cluster/confinement_20261013/conf_Afix_pilot.sbatch
```

It runs on `sandbox`: 20 trajectories (the method-A pilot's tasks, held). The last log line must be `cell pilot_epi8_H_H10_L10 done; failures: 0`.

**5. On the Mac:** fetch the pilot cell in full and run gate G1.

```sh
cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
bash hspist3/cluster/confinement_20261013/fetch_afix.sh pilot
cd hspist3 && python3 cluster/afix_pilot_check_261004.py
```

Paste the output. Its last line must be `**GATE G1: PASS**`. The fetch copies about 0.4 GB (the 20 event logs).

**6. After G1 and the go, from `login-0102` (at most 64 cores at once):**

```sh
cd ~/harddisks/hspist3
sbatch --array=4 --time=1:10:00 cluster/confinement_20261013/conf_Afix_0.10.sbatch
sbatch --array=1,2,3,5,6,7,8,9,10%1 cluster/confinement_20261013/conf_Afix_0.10.sbatch
sbatch --array=4 --time=3:05:00 cluster/confinement_20261013/conf_Afix_0.39.sbatch
sbatch --array=1,2,3,5,6,7,8,9%1 cluster/confinement_20261013/conf_Afix_0.39.sbatch
```

(Corrected 2026-10-04: the two rest arrays run at `%1`. With `%2`, the four jobs could hold 96 cores at once, not the 64 stated; the H40 lines are separate jobs. Now: 4 × 16 = 64.) (Amended 2026-10-04.) Every `--time` is 2 × the measured time of the identical method-A cell, rounded up to 5 min, with a minimum of 0:30. Round 1 and Round 2 ran these cells with the released divider: the same steps and the same collisions. In total the arrays use about 49 core-hours.

Printed by `python3 hspist3/cluster/gen_afix_sbatch_261004.py`:

| group | task | cell | N_s | trajectories | same (position, seed) set as method A | measured wall (h) | source | --time = 2 x measured | core-h (wall x 16) |
|---|---|---|---|---|---|---|---|---|---|
| Afix_0.10 | 1 | e0p10_H_H5_L39.25 | 25 | 710 | yes | 0.032 | Round 1 sacct | 0:30:00 | 0.5 |
| Afix_0.10 | 2 | e0p10_H_H10_L39.25 | 50 | 720 | yes | 0.051 | Round 1 sacct | 0:30:00 | 0.8 |
| Afix_0.10 | 3 | e0p10_H_H20_L39.25 | 100 | 675 | yes | 0.112 | Round 1 sacct | 0:30:00 | 1.8 |
| Afix_0.10 | 4 | e0p10_H_H40_L39.25 | 200 | 720 | yes | 0.548 | Round 1 32:06 (707/720) + repair 0:47 | 1:10:00 | 8.8 |
| Afix_0.10 | 5 | e0p10_L_H10_L19.625 | 25 | 750 | yes | 0.034 | Round 1 sacct | 0:30:00 | 0.6 |
| Afix_0.10 | 6 | e0p10_L_H10_L78.5 | 100 | 675 | yes | 0.108 | Round 1 sacct | 0:30:00 | 1.7 |
| Afix_0.10 | 7 | e0p10_aspect_H19.7917_L19.7917 | 50 | 730 | yes | 0.063 | Round 1 sacct | 0:30:00 | 1.0 |
| Afix_0.10 | 8 | e0p10_aspect_H14_L28 | 50 | 685 | yes | 0.055 | Round 1 sacct | 0:30:00 | 0.9 |
| Afix_0.10 | 9 | e0p10_aspect_H9.91667_L39.625 | 50 | 680 | yes | 0.053 | Round 1 sacct | 0:30:00 | 0.8 |
| Afix_0.10 | 10 | e0p10_aspect_H7_L56.0417 | 50 | 685 | yes | 0.052 | Round 1 sacct | 0:30:00 | 0.8 |
| Afix_0.10 | -- | sbatch default --time = 2 x the largest non-H40 cell = 0:30:00; H40 override 1:10:00 | | | | | | | |
| Afix_0.39 | 1 | epi8_H_H5_L10 | 25 | 130 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | 2 | epi8_H_H10_L10 | 50 | 115 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | 3 | epi8_H_H20_L10 | 100 | 130 | yes | 0.142 | Round 2 sacct | 0:30:00 | 2.3 |
| Afix_0.39 | 4 | epi8_H_H40_L10 | 200 | 255 | yes | 1.533 | Round 2 sacct | 3:05:00 | 24.5 |
| Afix_0.39 | 5 | epi8_L_H10_L5 | 25 | 130 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | 6 | epi8_L_H10_L20 | 100 | 130 | yes | 0.095 | Round 2 sacct | 0:30:00 | 1.5 |
| Afix_0.39 | 7 | epi8_aspect_H7.08333_L14.125 | 50 | 130 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | 8 | epi8_aspect_H5_L20 | 50 | 115 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | 9 | epi8_aspect_H3.54167_L28.2917 | 50 | 105 | yes | 0.035 | Round 2 sacct, range 0:55-2:07 (upper end) | 0:30:00 | 0.6 |
| Afix_0.39 | -- | sbatch default --time = 2 x the largest non-H40 cell = 0:30:00; H40 override 3:05:00 | | | | | | | |

total (measured wall x 16 cores, an upper bound): Afix_0.10 17.7 core-h, Afix_0.39 31.7 core-h, together 49.4 core-h

**7. When `squeue -u charing` is empty:** run `bash cluster/check_cells.sh` in a sandbox session (it covers the A-fixed groups `Afix_*`) and paste its `== Afix_*` blocks and `SUMMARY:` line. Then, on the Mac, fetch the summaries:

```sh
bash hspist3/cluster/confinement_20261013/fetch_afix.sh
```

The analysis is a separate task (261012 § 3.4).

## 10. Engine cost profile (written 2026-10-05; one short sandbox job, no build, no code change)

**What it is.** One held-divider trajectory at π/8 with N_s = 200 (H = 40) and one with N_s = 50 (H = 10), 700 σ-time each, run with the recorded build. It measures where the time goes, to decide how large the melting size sweep can be. The code itself already shows the expensive step (`edmd_core/edmd.c`):
- **Every disk–disk collision** rebuilds the grid and re-schedules both disks against all N partners (`:1529–1532`, `:791–794`). That is O(N) per event.
- **Every divider collision** calls `reschedule_all_internal(S)` (`:1537–1539`), which re-schedules all N²/2 pairs (`:812–816`). That is O(N²) per divider event, even though a held divider cannot move.

**1. Mac:** push (`bash _commit_scripts/commit_20261007.sh`). On KOA nothing needs a rebuild (rule 4), but the script must be pulled: in a sandbox session, `cd ~/harddisks/hspist3 && git pull && exit`.

**2. KOA, from a `login-…` prompt:**

```sh
cd ~/harddisks/hspist3
sbatch cluster/profile_edmd_koa.sh
squeue -u charing
```

**3. When `squeue` is empty** (a few minutes):

```sh
cat logs/profile-edmd_*.out
```

Paste the whole output. It shows, for each run, the wall time, the divider and wall events per second, and, if `perf` is usable on KOA, the top 15 functions.

## 11. Engine gate for the minimal divider rescheduling (written 2026-10-05; branch `engine-divider-resched`; 261012 § 4.4)

**What it is.** A new build generation: after a disk–divider collision the engine now re-schedules only that disk (plus every disk's divider event if the divider moved), instead of all N²/2 pairs. It runs from a **second clone**, `~/harddisks_resched`, and writes to a **second data root**, `/mnt/lustre/koa/scratch/charing/harddisks_resched/`. `~/harddisks`, its 279282b binary and its data are not touched. Nothing here may be merged into main before the plan author reads the gate.

**Rules for this step.**
- Before every `sbatch`: run `squeue -u charing` and `sacct -u charing -S today --format=JobID,JobName%20,State,Elapsed`, and never submit the same job twice.
- At most about 32 cores at once: the replay uses 3 × 8 = 24.
- On the login node, only `cd`, `ls`, `cat`, `tail`, `squeue`, `sacct` and `sbatch`.

**0. Mac (repo root):** push main and the branch, and note the branch hash.

```sh
git push origin main
git push origin engine-divider-resched
git log --oneline -1 engine-divider-resched
```

**1. KOA, sandbox session (`cn-…` prompt): clone the branch and build it.** Do the clone exactly as in step 1, into a new folder:

```sh
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

then, at the `cn-...` prompt:

```sh
cd ~
git clone --depth 1 --branch engine-divider-resched --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_resched
cd harddisks_resched
git sparse-checkout set --no-cone "/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"
git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft
git -C hspist3/kissfft -c advice.detachedHead=false checkout febd4caeed32e33ad8b2e0bb5ea77542c40f18ec
git log --oneline -1
cd hspist3 && mkdir -p logs
bash cluster/build_koa.sh
./00ALLINONE --version
exit
```

- **Expected:**
  - `git log` prints the hash of step 0;
  - the build ends with `BUILD OK`;
  - `--version` prints `00ALLINONE  git <that hash>  target koa`, with no `-dirty`.

**2. KOA, login node: the smoke test (G-E1, 0.06903 gate).** It also rebuilds, and that is fine.

```sh
cd ~/harddisks_resched/hspist3
sbatch cluster/koa_smoketest.sh
```

When `squeue -u charing` is empty, run `tail -30 logs/conf-smoke_*.out`.

- **Expected:**
  - `determinism self-test ...: IDENTICAL`;
  - five `PASS` gate lines;
  - `SMOKE TEST PASSED`.

**3. Cross-node determinism (G-E1).**

```sh
sbatch cluster/koa_crossnode_det.sh
```

Then run `tail -5 logs/det-xnode_*.out`.

- **Expected:** `... different nodes): IDENTICAL`.

**4. Minimal vs legacy, same binary and same seed (G-E2).** About 5 min, 6 cores (it also runs the 279282b binary of `~/harddisks`, read only).

```sh
sbatch cluster/resched_gate_261005/ge2.sbatch
```

Then run `cat logs/resched-ge2_*.out`.

- **Expected:** a table and the last line `G-E2: energy PASS; ledger PASS; health PASS; policy PASS; contact PASS`.
- Byte identity is expected to say `no` for minimal vs legacy, and `IDENTICAL` for legacy vs 279282b (261012 § 4.4).
- **Stop here and paste the output** if any part says FAIL, or if an `afix_` run shows an `[EDMD-HEALTH]` line. Do not go on to step 6.

**5. Profile, both policies on one node (G-E5).** About 5 min, 1 core.

```sh
sbatch cluster/profile_edmd_koa.sh
```

Then run `tail -12 logs/profile-edmd_*.out`.

- **Expected:** the last table, with `factor` and `p` columns.

**6. The three replayed cells (G-E3, G-E4).** Shared partition, 3 × 8 cores, up to 1 h.

```sh
sbatch --array=1-3 cluster/resched_gate_261005/replay.sbatch
```

When done, run `grep -h "done; failures" logs/resched-replay_*.out`.

- **Expected:** three lines, each ending `failures: 0`.
- **If a line shows failures:** stop and paste `grep -h "FAILED" logs/resched-replay_*.out`. Do not resubmit. CC will propose rerunning the failed seed with `--legacy-resched` on the same binary, which tells an old counter from a new defect.

**7. Mac (repo root, on main): copy back.** The fetch script and the analysis are on main too (identical copies), so no branch switch is needed.

```sh
bash hspist3/cluster/resched_gate_261005/fetch_resched.sh
```

Then tell CC. The analysis runs `cd hspist3 && python3 validation/resched_gate_261005.py`. It prints G-E2 to G-E5 and checks that every output comes from one clean build.

**If the old root is ever written to again:** after a merge, the worker's root guard refuses a data root that holds data but has no `.build_generation` record. For the 279282b root the record is one line, written by hand on purpose:

```sh
printf '%s\n' "00ALLINONE  git 279282b  target koa" > /mnt/lustre/koa/scratch/charing/harddisks/hspist3/.build_generation
```

## 12. Engine gate, version 2 (written 2026-10-06 HST = 2026-10-07 on the plan author's clock; 261012 § 4.4.10–4.4.11)

**What it is.** The second, cleaner exam of the faster engine.
- **The build:** 73fc07f plus one debugging switch (`--resched-audit`). Nothing else in the engine changed.
- **Where it runs:** from a **third clone**, `~/harddisks_resched2`, into its own data root, `/mnt/lustre/koa/scratch/charing/harddisks_resched2/`.
- **Untouched:** `~/harddisks` (279282b) and `~/harddisks_resched` (73fc07f). Never `git pull` in either.
- **The plan author's go** applies once the amendments of § 4.4.11 are committed and pushed.

**Rules for this step.**
- Run `squeue -u charing` before every `sbatch`, and never submit the same job twice.
- At most 32 cores at once.
- On the login node, only `cd`, `ls`, `cat`, `tail`, `grep`, `squeue`, `sacct` and `sbatch`.
- If anything below is not what is expected: stop and paste it.

**0. Mac (repo root).**

```sh
git push origin main
git push origin engine-divider-resched
git log --oneline -1 engine-divider-resched
```

Note the hash from the last line; KOA must show the same hash.

**1. KOA, sandbox session: clone and build.**

```sh
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

Then, at the `cn-...` prompt:

```sh
cd ~
git clone --depth 1 --branch engine-divider-resched --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_resched2
cd harddisks_resched2
git sparse-checkout set --no-cone "/hspist3/*.c" "/hspist3/*.h" "/hspist3/*.py" "/hspist3/Makefile" "/hspist3/edmd_core/" "/hspist3/cluster/" "/hspist3/validation/" "/hspist3/kissfft"
git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft
git -C hspist3/kissfft -c advice.detachedHead=false checkout febd4caeed32e33ad8b2e0bb5ea77542c40f18ec
git log --oneline -1
cd hspist3 && mkdir -p logs
bash cluster/build_koa.sh
./00ALLINONE --version
exit
```

- **Expected:**
  - `git log` prints the hash of step 0;
  - the build ends with `BUILD OK`;
  - `--version` prints `00ALLINONE  git <that hash>  target koa`.

**2. Same-node determinism (E1).** From `~/harddisks_resched2/hspist3`:

```sh
sbatch cluster/koa_smoketest.sh
```

When `squeue -u charing` is empty:

```sh
tail -30 logs/conf-smoke_*.out
```

- **Expected:** `determinism self-test ...: IDENTICAL`.
- **Expected, and NOT a failure of gate v2:** `c_s = 3.74424 +- 0.06756`, the same per-mass table as the 73fc07f smoke test (M = 50: 0.07040729 … M = 2000: 0.01506197), and `SMOKE TEST FAILED -- STOP`. It is the identical trajectory, because the audit switch is off. Carry on to step 3.
- **If c_s is any other number:** E0 is violated. Stop and paste.
- **If determinism is not IDENTICAL:** stop.

**3. Cross-node determinism (E1).**

```sh
sbatch cluster/koa_crossnode_det.sh
```

Then:

```sh
tail -5 logs/det-xnode_*.out
```

- **Expected:** `... different nodes): IDENTICAL`.

**4. E0 and E2.** Sandbox, 8 cores, up to 1 h.

```sh
sbatch cluster/resched_gate_261005/e0e2.sbatch
```

Then:

```sh
cat logs/resched-e0e2_*.out
```

- **Expected:** a table whose judged rows all end `ok`, a line of minimum sizes all `yes`, and the last line `E2 (amended): PASS; E0 (...): PASS`.
- Non-zero numbers in the column `abs(dt) > 1e-9 (engine)` are rounding and are allowed by the amended rule.

**5. Profile (E4).**

```sh
sbatch cluster/profile_edmd_koa.sh
```

Then:

```sh
tail -14 logs/profile-edmd_*.out
```

- **Expected:** a table with rows `held`, `free` and `dense`, and no `INVALID` line.

**6. Test T (E3).** Shared partition, at most 32 cores, about 8 core-hours.

```sh
sbatch --array=1-9%4 cluster/resched_gate_261005/testT.sbatch
```

When `squeue -u charing` is empty:

```sh
grep -h "task list SHA-256\|done; failures" logs/resched-testT_*.out
```

- **Expected:**
  - nine `task list SHA-256 60a00704e168d95bde3bac28ecd7fd4257428d94989c724d3ca294d5076e6be1 (expected 60a00704…)` lines;
  - nine `testT M=… done; failures: 0` lines.
- **If any line shows failures:** stop, paste `grep -h FAILED logs/resched-testT_*.out`, and do not resubmit.

**7. Mac (repo root, on main): copy back.**

```sh
bash hspist3/cluster/resched_gate_261005/fetch_resched2.sh
```

Then tell CC. The analysis is `python3 validation/resched_testT_261007.py` plus the audit report.

## 13. Engine gate, version 3: Test T-prime and ASan (written 2026-10-07 HST; 261012 § 4.4.13)

**What it is.** The last test of the faster engine, and a memory checker.
- **Test T-prime:** the same comparison as Test T, but only at M = 300 and M = 1500, with 400 fresh seeds per path. That is 1600 trajectories, about 11 core-hours.
- **ASan:** a separate debug build of the same code checks every memory access in four short runs.
- **Where:** the same clone (`~/harddisks_resched2`) and the same binary as Test T. **Nothing is rebuilt.** The data go to the same data root.
- **Untouched:** `~/harddisks` (279282b) and `~/harddisks_resched` (73fc07f). Never `git pull` in either.

**Rules for this step.**
- Run `squeue -u charing` before every `sbatch`, and never submit the same job twice.
- At most 32 cores at once. This step uses 16 + 4.
- On the login node, only `cd`, `ls`, `cat`, `tail`, `grep`, `squeue`, `sacct` and `sbatch`.
- If anything below is not what is expected: stop and paste it, and do not resubmit.

**0. Mac (repo root).**

```sh
git push origin main
git push origin engine-divider-resched
git log --oneline -1 engine-divider-resched
```

Note the hash from the last line; KOA must show the same hash in step 1.

**1. KOA, sandbox session: bring the clone up to date; check that the binary is still Test T's.**

```sh
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

Then, at the `cn-...` prompt:

```sh
cd ~/harddisks_resched2
git pull --ff-only origin engine-divider-resched
git log --oneline -1
git diff --stat 7b08827 HEAD -- hspist3/00ALLINONE.c hspist3/edmd_core hspist3/experiment_validation.c hspist3/experiment_validation.h hspist3/Makefile
cd hspist3
sha256sum -c logs/BUILD_KOA_LAST.hash
./00ALLINONE --version | head -1
sha256sum cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt
exit
```

- **Expected:**
  - `git log` prints the hash of step 0;
  - `git diff --stat` prints **nothing**: the engine is unchanged since 7b08827;
  - `00ALLINONE: OK`;
  - `00ALLINONE  git 7b08827  target koa`;
  - `cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f  cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt`.
- **Do not run `make` or `build_koa.sh`.** The binary of Test T must stay.

**2. Test T-prime (login node).** Shared partition, 2 tasks × 8 cores.

```sh
squeue -u charing
cd ~/harddisks_resched2/hspist3
sbatch --array=1-2%2 cluster/resched_gate_261005/testTprime.sbatch
squeue -u charing
```

- **Expected:**
  - the first `squeue` shows only its header line;
  - `Submitted batch job <id>`;
  - the second `squeue` shows `<id>_1` and `<id>_2`.

**3. ASan (login node).** Sandbox partition, 4 cores, up to 1 h.

```sh
squeue -u charing
sbatch cluster/resched_gate_261005/asan.sbatch
```

- **Expected:** `squeue` shows only the two T-prime tasks; then `Submitted batch job <id>`.

**4. Checks, when `squeue -u charing` is empty.** M = 1500 takes about 1 h; ASan about 10–20 min.

```sh
grep -h "task list SHA-256\|done; failures\|STOP" logs/resched-testTprime_*.out
tail -n 40 logs/resched-asan_*.out
```

- **Expected, T-prime:**
  - two `task list SHA-256 cdef566b… (expected cdef566b…)` lines;
  - `testTprime M=300 done; failures: 0`;
  - `testTprime M=1500 done; failures: 0`;
  - no `STOP`.
- **Expected, ASan:**
  - `sanitizer build: 00ALLINONE  git 7b08827  target asan-scratch`;
  - `recorded build:  00ALLINONE  git 7b08827  target koa`;
  - four lines `smoke_min: exit 0 (… s)` … `tT300_leg: exit 0 (… s)`;
  - `recorded binary and hash file unchanged`;
  - a table whose rows all end `yes`;
  - the last line `ASAN (decision 3, item 3): CLEAN -- zero sanitizer reports, exit 0, audit missing = extra = 0 in all four runs`.
- **If any line shows a failure, `STOP`, `NOT CLEAN` or `WARNING`:** stop and paste it. Do not resubmit; any sanitizer report decides nothing until the plan author has read it.

**5. Mac (repo root, on main): copy back.**

```sh
bash hspist3/cluster/resched_gate_261005/fetch_resched2.sh
```

Then tell CC. The verdict is printed by `python3 validation/resched_testTprime_261007.py`, by the registered rule, including the ASan report.

## 14. KOA cap (plan author's decision of 2026-10-09 on the plan author's clock, written 2026-10-08 HST; 261012 § 4.7.1 item 3)

The cap is our own choice. KOA's limits (printed by Chris on 2026-10-08):
- no per-user CPU or job limit;
- MaxSubmit 60001;
- MaxArraySize 25001;
- `shared` has 3,012 CPUs on 90 nodes.

The cap:
- **64 cores standing.** This is the sum over every running job of `charing`.
- **128 cores only while** `sinfo -p shared -s` shows at least 10 % of the nodes idle. In `NODES(A/I/O/T)`, I / T ≥ 0.10; with T = 90 that is I ≥ 9.
- **One trajectory per array task,** one core per task.
- **Arrays chunked below 25001 tasks.** Throttle with `%` so that the running tasks stay within the cap, e.g. `--array=1-20000%64` for one-core tasks.
- **Before every sbatch:** run `squeue -u charing` (never submit twice). For the 128 cap, also run `sinfo -p shared -s`.
- **The sshare value is not interpreted** (plan author).
- **Generation 3:** no KOA runs until M6, the gate (261012 § 4.7.1 item 5).
