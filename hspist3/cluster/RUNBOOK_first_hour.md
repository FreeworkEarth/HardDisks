# KOA, the first hour

Written 2026-10-01, before the account exists. Each step has a **stop condition**: if it does not
hold, fix it before the next step rather than pushing on. Nothing here runs science — the point of
the first hour is to prove the path works end to end on something cheap.

Facts assumed: `sandbox` 4 h, `shared` 3 d / ≤ 93 cores per node, `kill-shared` preemptable,
home 50 GiB unbacked, `~/koa_scratch` unquota'd and purged at 90 days, DTN
`koa-dtn.its.hawaii.edu`, Lmod, heterogeneous CPUs including Ivy Bridge. Confirm them in step 2;
they came from a briefing, not from the machine.

---

## 1. Get on

```sh
ssh <uhid>@koa.its.hawaii.edu
```

**Stop condition:** a shell prompt. If Duo/2FA blocks the DTN later, it is the same credential —
sort it out now, not mid-transfer.

## 2. Confirm the facts, and write down what you find

```sh
sinfo -o "%P %l %D %c %m" | head -20
sacctmgr show assoc where user=$USER format=Account,Partition,MaxJobs,GrpTRES | head
lscpu | grep -E "Model name|Flags" | head -2
df -h $HOME ~/koa_scratch
module avail 2>&1 | grep -iE "sdl|glew|mesa|gcc" | head
```

**Stop condition:** you can name the partition to submit to, your core limit, and the oldest CPU in
it. Paste the output into `cluster/KOA_FACTS_<date>.txt` — later confusion about which node a
number came from is answered from that file.

> Look for `sse4_2` and `popcnt` in Flags. If `avx2` is **absent** on any node, that is the Ivy
> Bridge case and confirms `-march=x86-64-v2` was the right call.

## 3. Hello, scheduler

```sh
mkdir -p ~/harddisks/logs && cd ~/harddisks
sbatch cluster/hello.sbatch
squeue -u $USER
cat logs/hello_*.out
```

**Stop condition:** the job runs on `sandbox` and prints its hostname and CPU model. This proves
submission, the log path, and that `logs/` exists — the three things that silently break a first
array.

## 4. Move the source over the DTN, not the login node

```sh
# FROM THE LAPTOP
rsync -av --exclude 'experiments_*' --exclude '*.o' --exclude '00ALLINONE*' \
      hspist3/ <uhid>@koa-dtn.its.hawaii.edu:~/harddisks/hspist3/
```

**Stop condition:** `ls ~/harddisks/hspist3/00ALLINONE.c` on KOA. Note what is *excluded*: no
binaries (built there), no campaign output (lives in scratch). Source and scripts belong in home;
**nothing you cannot regenerate should ever be only on KOA**, which is unbacked.

## 5. Build for the cluster

```sh
cd ~/harddisks/hspist3
module load gcc          # whatever step 2 showed
make koa                 # -O2 -march=x86-64-v2, portable across the partition
bash cluster/build_koa.sh   # records compiler, flags, git commit, binary hash
./00ALLINONE --help | head -3
```

**Stop condition:** `--help` prints. If it dies with **Illegal instruction**, the ISA is too high
for this node — that is exactly the failure `-march=x86-64-v2` exists to prevent, so check
`CFLAGS_KOA` before doing anything else.

> `make` alone builds the AddressSanitizer debug target and is ≈ 5× slower. Never use it here.

## 6. One cell, interactively, on sandbox

```sh
srun -p sandbox -t 0:30:00 --cpus-per-task=1 --pty bash
cd ~/harddisks/hspist3
mkdir -p ~/koa_scratch/smoke
./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
  --edmd-acc=0 --seed-drift-order=drift-first \
  --energy-transfer-summary=$HOME/koa_scratch/smoke/summary.csv \
  --energy-transfer-trace=$HOME/koa_scratch/smoke/trace.csv --trace-every=200 \
  --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
  --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
  --wall-mass-factors=10 --eff-output=wall-ke \
  --wall-hold-steps=12000 --steps=120000 --fixed-dt=0.4 --kbt1 --seed=9200
ls -la ~/koa_scratch/smoke/
exit
```

**Stop condition:** a trace of a few hundred rows, and `grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'` returning **0**.
The health contract is the same on every machine; it is not relaxed for the cluster.

## 7. Ten tasks as an array

```sh
sbatch --array=1-10%10 --partition=shared cluster/level4_massladder.sbatch 10
squeue -u $USER
seff <jobid>            # after it finishes: CPU efficiency and peak memory
```

**Stop condition:** 10/10 succeed, `seff` shows CPU efficiency near 100 % and memory well under the
requested `--mem`. If efficiency is low, the job is waiting on I/O — lower the trace cadence
(raise `--trace-every`) rather than asking for more cores; each task is one core by construction.

Now repeat on `kill-shared` and kill one task by hand (`scancel <jobid>_3`), then resubmit the same
array. **Stop condition:** the resubmission skips the nine completed trajectories and redoes only
the third. That is the restart-safety that makes the preemptable partition usable; if it does not
hold, do not use `kill-shared`.

## 8. Bring one cell home

```sh
# FROM THE LAPTOP
bash cluster/fetch_results.sh <uhid>@koa-dtn.its.hawaii.edu ~/koa_scratch/level4
```

**Stop condition:** the files are on the laptop and readable by the analysis scripts. Do this once
on day one — not for the first time in month three, when the purge clock matters.

---

## What day one does NOT do

No science campaign, and **no comparison of a KOA number with a laptop number**. Cluster results are
a separate campaign: different architecture, different FMA contraction, and chaotic dynamics
diverging from the last bit. Acceptance is the statistical mirror gate in `README.md`, run after the
first hour succeeds — never byte identity, which must not be claimed across architectures.
