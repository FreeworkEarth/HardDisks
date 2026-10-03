#!/bin/bash
# ##CHRIS 2026-10-01: the smallest thing that proves submission works. Step 3 of RUNBOOK_first_hour.
#
#   sbatch cluster/hello.sbatch
#
# It prints the three facts every later provenance line depends on: which node, which CPU, which
# partition. If the CPU model here lacks avx2, that is the Ivy Bridge case and -march=x86-64-v2 in
# the `koa` Makefile target is doing its job.
#
#SBATCH --job-name=hello
#SBATCH --partition=sandbox
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

set -u
echo "hostname   : $(hostname)"
echo "cpu        : $(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ *//')"
echo "avx2       : $(grep -qm1 avx2 /proc/cpuinfo && echo yes || echo NO -- Ivy Bridge class, v2 required)"
echo "sse4_2     : $(grep -qm1 sse4_2 /proc/cpuinfo && echo yes || echo NO -- v2 will NOT run here)"
echo "partition  : ${SLURM_JOB_PARTITION:-none}"
echo "job        : ${SLURM_JOB_ID:-none}"
echo "cores seen : $(nproc)"
echo "home free  : $(df -h "$HOME" | awk 'NR==2{print $4}')"
echo "scratch    : $(df -h "$HOME/koa_scratch" 2>/dev/null | awk 'NR==2{print $4}' || echo 'MISSING -- create ~/koa_scratch')"
echo "date       : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
