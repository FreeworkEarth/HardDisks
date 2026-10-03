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

## 1. Get the code (login node, ~10 s)

Do this once, after you pushed from the Mac. It fetches only the sources and scripts (about 11 MB), never binaries
or data.

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
- **Later updates:** `cd ~/harddisks && git pull` (after a push from the Mac). Then rebuild in step 2.

## 2. Build once, interactively (compute node, ~5 min)

```sh
srun -p sandbox -t 1:00:00 -c 2 --mem=4G --pty /bin/bash
```

Wait for the `cn-...` prompt, then:

```sh
cd ~/harddisks/hspist3
bash cluster/build_koa.sh
```

- `srun ... --pty`: a shell on a compute node, 2 cores, at most 1 hour.
- `build_koa.sh`:
  - loads the compiler (module `compiler/GCC/14.3.0`) and the libraries from `~/envs/hd` (`cluster/koa_env.sh`);
  - compiles with `-O2 -march=x86-64-v2 -ffp-contract=off`;
  - checks that every library resolves (`ldd`) and that the binary's hash is the clean checkout's HEAD;
  - writes `logs/BUILD_KOA_<jobid>.txt`.
- **Expected:**
  - the first line is `env: gcc (GCC) 14.3.0 | Python 3.12.14 at /home/charing/envs/hd/bin/python3 | git version ...`;
  - then the compiler output (warnings are fine);
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
