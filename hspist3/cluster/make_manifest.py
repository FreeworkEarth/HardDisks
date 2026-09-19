#!/usr/bin/env python3
"""##CHRIS 2026-09-16: build a run manifest for the KOA cluster (one line per trajectory).

Run this on the laptop, where the geometry and seed helpers live; the cluster then needs nothing but
bash and the binary. Every line is a complete, independent run with its own exact seed and output
directory, which is what makes a SLURM array the right shape: one array task per line, restartable,
and skippable when its trace already exists.

Campaigns
  mirror      three cells that already exist on the laptop, same seeds, same 37.5-period records.
              This is the acceptance gate for the cluster build (see README.md) -- not science.
  a2_long200  the full A2 ladder with 200-period records: 5 densities x 5 sizes x 5 masses x 10 seeds.
  n900_sweep  a second c_s(eta) curve at N = 900, 12 densities x 5 masses x 10 seeds, 200 periods.
  edmd_acc    the mirror cells run on both backends, for the accelerated-backend validation.

usage: make_manifest.py CAMPAIGN OUT.tsv [--root REMOTE_ROOT]
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "validation"))
import tests_20260913 as T  # noqa: E402  (k_root, kr_cs, d_stride, run_seed)

MASSES = (50, 200, 500, 1000, 2000)
COLS = ("idx", "campaign", "outdir", "particles", "boxes", "height", "radius", "wall",
        "lengths", "wall_mass", "target_osc", "safety", "stride", "max_steps", "edmd_acc", "seed_base", "exact_seed")


def a2_geometry(eta, N):
    """famB convention: L0 and H both scale as sqrt(N/100), so the aspect ratio is fixed."""
    f = math.sqrt(N / 100.0)
    return f"{3.926990816987241 / eta * f:.6f}", f"{10 * f:.6f}"


def x_of(M, N, L0):
    return T.k_root(M / float(N)) / (2 * math.pi * (L0 - 1.0))


def stride_for(eta, N, L0, M, target):
    """~32 samples per predicted period, as in A1 v2 (TEST C showed 8x coarser changes nothing)."""
    nu = T.kr_cs(eta) * x_of(M, N, L0) if eta <= 0.69 else 0.0
    return T.d_stride(nu) if nu > 0 else "auto"


def cells(campaign):
    """-> list of (eta, N, target_oscillations, safety, seed_base, edmd_acc, tag)"""
    if campaign == "mirror":
        # exactly the famB settings of three existing laptop cells
        return [(0.30, 400, 25, 1.5, 24010000, 0, "mirror"),
                (0.50, 1600, 25, 1.5, 24010001, 0, "mirror"),
                (0.65, 900, 25, 1.5, 24010002, 0, "mirror")]
    if campaign == "edmd_acc":
        out = []
        for i, (eta, N) in enumerate(((0.30, 400), (0.50, 1600), (0.65, 900))):
            for acc in (0, 1):
                out.append((eta, N, 25, 1.5, 24020000 + 10 * i + acc, acc, f"acc{acc}"))
        return out
    if campaign == "a2_long200":
        return [(eta, N, 200, 1.0, 24030000 + 100 * i + j, 0, "long200")
                for i, eta in enumerate((0.10, 0.30, 0.50, 0.60, 0.65))
                for j, N in enumerate((100, 400, 900, 1600, 2500))]
    if campaign == "n900_sweep":
        etas = (0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.68)
        return [(eta, 900, 200, 1.0, 24040000 + i, 0, "n900") for i, eta in enumerate(etas)]
    raise SystemExit(f"unknown campaign {campaign}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign")
    ap.add_argument("out")
    ap.add_argument("--root", default="$SCRATCH/harddisks/runs",
                    help="campaign root ON THE CLUSTER; expanded by the sbatch script, not here")
    ap.add_argument("--seeds", type=int, default=10)
    a = ap.parse_args()

    lines, idx = [], 0
    for eta, N, target, safety, base, acc, tag in cells(a.campaign):
        L0, H = a2_geometry(eta, N)
        for mi, M in enumerate(MASSES):
            stride = stride_for(eta, N, float(L0), M, target)
            for r in range(a.seeds):
                idx += 1
                out = (f"{a.root}/{a.campaign}/eta_{('%.2f' % eta).replace('.', 'p')}/N{N}"
                       f"{'' if acc == 0 else '_acc1'}/m_{M}/r{r}")
                lines.append((idx, a.campaign, out, N, f"{N // 2},{N // 2}", H, "0.5", "0.05",
                              L0, M, target, f"{safety:.1f}", stride, 400000000, acc, base,
                              T.run_seed(base, 0, mi, r)))
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for ln in lines:
            fh.write("\t".join(str(v) for v in ln) + "\n")
    cell_count = len(cells(a.campaign))
    print(f"{a.out}: {len(lines)} runs over {cell_count} cells "
          f"({a.seeds} seeds x {len(MASSES)} masses), campaign {a.campaign}")
    print(f"submit with:  sbatch --array=1-{len(lines)}%200 cluster/run_array.sbatch {os.path.basename(a.out)}")


if __name__ == "__main__":
    main()
