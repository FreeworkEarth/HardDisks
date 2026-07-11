# EDMD Accelerated Core (`edmd_acc_*`)

This repo contains **two EDMD backends**:

- **Default EDMD**: `hspist3/edmd_core/edmd.c` (baseline, simpler scheduler)
- **Accelerated EDMD**: `hspist3/edmd_core/edmd_accelerated.c` (faster scheduler for larger `N`)

Both are compiled into `00ALLINONE` by default, and you choose at runtime via a CLI flag.

## Why an accelerated core?

In EDMD (event-driven molecular dynamics) you advance the system from **collision to collision**.
The expensive part is **finding the next particle–particle collision events**.

If you naïvely schedule particle–particle events for all pairs, that’s **O(N²)** work when (re)building the event schedule.
At `N≈600` this becomes painfully slow.

The accelerated core replaces the “all-pairs” scheduling with a **uniform grid + neighbor-only scheduling** approach that is closer to the “collision window” intuition from time-stepping codes.

## Key idea (uniform grid + cell-crossing events)

Particles are binned into a **uniform grid** (spatial hashing).

For each particle `a`, we only schedule candidate collisions with particles in the **3×3 neighborhood** of `a`’s grid cell:

- same cell
- the 8 adjacent cells

This drastically reduces the number of candidate pairs in typical dense-ish gases.

### Correctness problem

If you only look at the current neighbor cells, you can miss collisions when a particle moves into a new cell.

### Fix: `EV_CC` (cell crossing) events

The accelerated core schedules an additional event type:

- **`EV_CC`**: “cell crossing” event = the next time a particle would cross a grid cell boundary.

When an `EV_CC` event fires, we:

1. advance positions to that time
2. rebuild the grid occupancy
3. reschedule that particle’s local events (neighbors + walls/dividers/pistons)

This is the mechanism that keeps neighbor-only scheduling consistent over time.

## What is actually accelerated?

The accelerated core mainly changes **event scheduling**, not the physical collision rules:

- particle–particle collisions (`EV_AB`)
- wall collisions (`EV_WL`, `EV_WR`, `EV_WB`, `EV_WT`)
- divider/piston collisions (`EV_DL`, `EV_DR`, `EV_PL`, `EV_PR`)
- plus `EV_CC` for grid maintenance

In practice, the big win is avoiding repeated global **N²** rebuilds.

## Exported API / symbol names

To allow both backends to be linked into the same binary, the accelerated core exports **prefixed symbols**:

- `edmd_acc_create`, `edmd_acc_advance_to`, ...

Header:

- `hspist3/edmd_core/edmd_accelerated.h`

The rest of `00ALLINONE.c` calls the canonical `edmd_*` names, but these calls are dispatched to either backend by a small wrapper layer (selected via `--edmd-acc`).

## How to use it (runtime flag)

Use the accelerated backend only in EDMD modes:

- `--mode=edmd`
- `--mode=edmd-hybrid`

Enable accelerated backend:

```bash
./00ALLINONE --mode=edmd --edmd-acc=1 ...
```

Disable (default backend):

```bash
./00ALLINONE --mode=edmd --edmd-acc=0 ...
```

You can also pass just `--edmd-acc` (equivalent to `--edmd-acc=1`).

`00ALLINONE` prints the selected backend in the “Simulation Info” section:

- `EDMD backend: accelerated` or `EDMD backend: default`

## How to build

Default: build both EDMD backends (so runtime can choose):

```bash
cd hspist3
make
```

If you want to build only one backend (useful for debugging):

```bash
make clean
make EDMD_SRCS=edmd_core/edmd.c
```

or:

```bash
make clean
make EDMD_SRCS=edmd_core/edmd_accelerated.c
```

## Notes / limitations

- This is a **scheduler acceleration**, not a change in the physical model.
- **Harmonic (spring) divider walls:** EDMD now supports a divider whose center position follows a
  harmonic oscillator between collisions: `m ẍ = -k (x - x_eq)`. This requires **root-finding** to
  compute particle–divider collision times because the divider is not constant-velocity.
  - Enabled by setting `EDMD_Params.divider_k[d] > 0` and `EDMD_Params.divider_mass[d] > 0`.
  - Works in both backends (default + accelerated), but is slower than constant-velocity dividers.
- Grid sizing matters: the implementation chooses a `cell_size` intended to keep neighbor checks small while still being correct with `EV_CC` maintenance.
- The accelerated core is newer: validate by running small-`N` cases where you can compare default vs accelerated outputs statistically.
- If you hit unexpected behavior, the fastest sanity check is: rerun the same config with `--edmd-acc=0` and see if the issue disappears.

## Where to look in code

- Accelerated core implementation: `hspist3/edmd_core/edmd_accelerated.c`
  - neighbor scheduling: `schedule_for(...)`
  - cell crossing events: `EV_CC` handling in `edmd_acc_advance_to(...)`
- Backend selection flag + dispatch: `hspist3/00ALLINONE.c`
  - CLI flag parsing: `--edmd-acc`
  - dispatch wrappers: `edmd_backend_*` helpers and `#define edmd_* ...` mapping
