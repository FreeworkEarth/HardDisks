#!/usr/bin/env python3
"""Screen existing sound-speed structure measurements; never run the simulator.

Outputs are diagnostic, not a certification of equilibrium or bulk sound speed.
Original leaf records are used instead of their merged copies. Source hashes
are checked again before outputs are written to detect concurrent modification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

CAMPAIGNS = {
    "campaign_r25_psi6_20260823": "Route A",
    "routeB_radius_N100_L0_20_20260825": "Route B",
    "ladder_N100_20260825": "Strip N100",
    "ladder_N200_20260825": "Strip N200",
    "ladder_N400_20260825": "Strip N400",
    "overnight_N1000_20260826": "Strip N1000",
    "finitesize_aspect_20260826": "Fixed aspect",
}
METRICS = ("psi6_global_hold", "psi6_global_end", "psi6_local_hold",
           "psi6_local_end", "neighbors_hold", "neighbors_end")
HEALTH = re.compile(r"\[EDMD-HEALTH\].*?M=(\d+) run=(\d+) seed=(\d+): (.*)")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_metrics(raw):
    """Allow only the documented undefined-order case: no neighbours at all."""
    row = {k: float(raw[k]) for k in METRICS}
    missing = []
    for phase in ("hold", "end"):
        names = [f"psi6_global_{phase}", f"psi6_local_{phase}"]
        if row[f"neighbors_{phase}"] == 0 and all(math.isnan(row[k]) for k in names):
            missing.append(f"undefined_psi6_{phase}_no_neighbors")
            for k in names:
                row[k] = None
    if any(v is not None and not math.isfinite(v) for v in row.values()):
        raise ValueError("Unexpected nonfinite structural measurement")
    for key, value in row.items():
        if value is not None and (value < 0 or (key.startswith("psi6") and value > 1.000001)):
            raise ValueError(f"Out-of-range structural measurement: {key}={value}")
    return row, "; ".join(missing)


def screen(rows):
    """Paired endpoint screen. Dual tails are NOT a bimodality test."""
    n = len(rows)
    result = {"n": n}
    for metric in METRICS:
        vals = [float(r[metric]) for r in rows]
        result[metric + "_mean"] = st.mean(vals)
        result[metric + "_sd"] = st.stdev(vals) if n > 1 else None
    delta = [r["psi6_global_end"] - r["psi6_global_hold"] for r in rows]
    result["delta_mean"] = st.mean(delta)
    result["delta_sem"] = st.stdev(delta) / math.sqrt(n) if n > 1 else None
    result["fraction_abs_delta_gt_0p1"] = sum(abs(x) > .1 for x in delta) / n
    result["end_lt_0p4"] = sum(r["psi6_global_end"] < .4 for r in rows)
    result["end_gt_0p8"] = sum(r["psi6_global_end"] > .8 for r in rows)
    result["endpoint_shift_flag"] = abs(result["delta_mean"]) > .1
    result["dual_tail_screen"] = (n >= 10 and
        min(result["end_lt_0p4"], result["end_gt_0p8"]) >= max(2, math.ceil(.1*n)))
    # Low neighbour counts are expected in dilute gas: not a square-order test.
    result["dense_low_neighbors_flag"] = (rows[0]["eta"] >= .6 and
        min(result["neighbors_hold_mean"], result["neighbors_end_mean"]) < 4.3)
    return result


def grouped(rows, keys):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    return [dict(zip(keys, key), **screen(vals))
            for key, vals in sorted(groups.items())]


def load_sources(root):
    manifest, records, ledgers = {}, [], []

    def read(path):
        raw = path.read_bytes()
        manifest[str(path.relative_to(root))] = hashlib.sha256(raw).hexdigest()
        return raw.decode()

    for folder, label in CAMPAIGNS.items():
        base = root / folder
        if not base.is_dir():
            raise ValueError(f"Missing required campaign: {base}")
        for path in sorted(base.rglob("speed_of_sound_psi6.csv")):
            if any(part in {"merged", "analysis", "analysis_pereta", "an"}
                   for part in path.relative_to(base).parts):
                continue
            parent = path.parent
            command = read(parent / "00_COMMAND.md")
            match = re.search(r"--particles=(\d+)\b", command)
            if not match:
                raise ValueError(f"Missing particle count: {parent}")
            n_total = int(match[1])
            acc = re.search(r"--edmd-acc=(\d+)", command)
            if acc and acc[1] != "0":
                raise ValueError(f"Accelerated core unexpectedly in scope: {parent}")
            family = label
            if label == "Fixed aspect":
                family += " " + path.relative_to(base).parts[0]
            status = json.loads(read(parent / "speed_of_sound_batch_status.json"))
            complete = (status.get("complete", True) and status["invalid_runs"] == 0
                        and status["valid_runs"] == status["requested_runs"]
                        and status.get("in_progress_runs", 0) == 0)
            log = read(parent / "run.log")
            health = {}
            for m in HEALTH.finditer(log):
                counters = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", m[4])}
                if any(counters.values()):
                    health[(int(m[1]), int(m[2]), int(m[3]))] = m[4]
            source_rows = list(csv.DictReader(io.StringIO(read(path))))
            if len(source_rows) != status["valid_runs"]:
                raise ValueError(f"Structure/ledger count mismatch: {path}")
            ledgers.append(dict(source=str(path.relative_to(root)), family=family,
                                N_total=n_total, records=len(source_rows),
                                complete=complete, health_warnings=len(health)))
            for index, raw in enumerate(source_rows, 2):
                try:
                    row, missing = parse_metrics(raw)
                except ValueError as exc:
                    raise ValueError(f"{path}:{index}: {exc}") from exc
                mass, repeat, seed = (int(raw[k]) for k in ("wall_mass_factor", "repeat", "seed"))
                reason = health.pop((mass, repeat, seed), "")
                row.update(family=family, N_total=n_total, eta=round(float(raw["eta"]), 6),
                           M=mass, repeat=repeat, seed=seed, L0=float(raw["L0"]),
                           source=str(path.relative_to(root)), source_line=index,
                           health_reason=reason, missing_structure_reason=missing,
                           ledger_complete=complete,
                           eligible_structure=complete and not reason and not missing,
                           measured_temperature=None)
                records.append(row)
            if health:
                raise ValueError(f"Unmatched health warning: {parent}: {health}")
    seen = set()
    for r in records:
        key = tuple(r[k] for k in ("family", "N_total", "eta", "M", "repeat", "seed"))
        if key in seen:
            raise ValueError(f"Duplicate source trajectory: {key}")
        seen.add(key)
    for name, sha in manifest.items():
        if digest(root / name) != sha:
            raise ValueError(f"Source changed during analysis: {name}; retry after writer stops")
    return records, ledgers, manifest


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing empty table: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plots(out, records, cells, masses):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180,
                         "pdf.fonttype": 42})
    families = sorted({r["family"] for r in cells})

    def save(fig, name):
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"{name}.{suffix}", bbox_inches="tight")
        plt.close(fig)

    def tidy_axis(ax, family):
        ax.ticklabel_format(axis="x", useOffset=False)
        if family == "Fixed aspect famC":
            # Nominally one density, not a density sweep at the 1e-6 level.
            ax.set_xlim(.705, .715)
            ax.set_xticks([.71])

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.4), sharey=True, layout="constrained")
    for ax, family in zip(axes.flat, families):
        rs = [r for r in cells if r["family"] == family]
        for n in sorted({r["N_total"] for r in rs}):
            rr = sorted((r for r in rs if r["N_total"] == n), key=lambda x: x["eta"])
            line, = ax.plot([r["eta"] for r in rr], [r["psi6_global_end_mean"] for r in rr],
                            "o-", ms=3, label=f"N={n}, end")
            ax.plot([r["eta"] for r in rr], [r["psi6_global_hold_mean"] for r in rr],
                    "x--", color=line.get_color(), lw=1.1, ms=4)
        ax.set(title=family, xlabel=r"Packing fraction $\eta$", ylim=(-.04, 1.04))
        tidy_axis(ax, family)
        ax.grid(alpha=.18); ax.legend(fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Mean global $|\psi_6|$")
    fig.suptitle("Structure before release and after measurement\nDashed: end of hold; solid: end of run. Endpoint agreement alone does not prove equilibrium.", fontsize=12)
    save(fig, "01_structure_endpoints")

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.4), sharey=True, layout="constrained")
    for ax, family in zip(axes.flat, families):
        rs = [r for r in cells if r["family"] == family]
        for n in sorted({r["N_total"] for r in rs}):
            rr = sorted((r for r in rs if r["N_total"] == n), key=lambda x: x["eta"])
            ax.errorbar([r["eta"] for r in rr], [r["delta_mean"] for r in rr],
                        yerr=[r["delta_sem"] or 0 for r in rr], fmt="o-", ms=3,
                        capsize=2, label=f"N={n}")
        ax.axhspan(-.1, .1, color="#aaaaaa", alpha=.12)
        ax.axhline(0, color="gray", lw=.7)
        ax.set(title=family, xlabel=r"Packing fraction $\eta$")
        tidy_axis(ax, family)
        ax.grid(alpha=.18); ax.legend(fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Paired $|\psi_6|_{end}-|\psi_6|_{hold}$")
    fig.suptitle("Structural change during the measurement\nError bars: descriptive SEM across trajectories; shaded band: ±0.1 screen, not an acceptance test.", fontsize=12)
    save(fig, "02_structure_change")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for ax, eta in zip(axes.flat, (.63, .70, .72, .75)):
        rr = [r for r in records if r["family"] == "Route B" and abs(r["eta"]-eta)<1e-5]
        ax.hist([r["psi6_global_hold"] for r in rr], bins=np.linspace(0,1,21),
                histtype="step", linewidth=1.8, label="End of hold")
        ax.hist([r["psi6_global_end"] for r in rr], bins=np.linspace(0,1,21),
                alpha=.5, label="End of run")
        ax.set(title=f"Route B: η={eta:.2f}, n={len(rr)}", xlabel=r"Global $|\psi_6|$", ylabel="Trajectories")
        ax.legend(); ax.grid(alpha=.15)
    fig.suptitle("Endpoint distributions in the confined geometry\nMasses pooled for this diagnostic; see the per-mass table before interpreting mixtures.", fontsize=12)
    save(fig, "03_routeB_distributions")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    for ax, family in zip(axes, ("Route A", "Route B")):
        for eta in (.63,.70,.72,.75):
            rr = sorted((r for r in masses if r["family"]==family and abs(r["eta"]-eta)<1e-5), key=lambda r:r["M"])
            ax.errorbar([r["M"] for r in rr], [r["delta_mean"] for r in rr],
                        yerr=[r["delta_sem"] or 0 for r in rr], fmt="o-", ms=4,
                        capsize=2, label=f"η={eta:.2f}")
        ax.axhspan(-.1,.1,color="gray",alpha=.1); ax.axhline(0,color="gray",lw=.7)
        ax.set(title=family, xscale="log", xlabel="Divider mass / particle mass",
               ylabel=r"Paired change in global $|\psi_6|$")
        ax.legend(); ax.grid(alpha=.15)
    fig.suptitle("Does structural change depend on the divider mass?\nError bars: seed SEM; mass also changes the planned measurement duration.",fontsize=12)
    save(fig, "04_mass_dependence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--refresh", action="store_true", help="Refresh only this script's known outputs")
    args = parser.parse_args()
    root, out = args.root.resolve(), args.out.resolve()
    if out.exists():
        known = {"README.md", "trajectories.csv", "structural_cells.csv", "structural_by_mass.csv",
                 "source_ledgers.csv", "source_manifest.json", "exclusions.csv"}
        known.update(f"{name}.{ext}" for name in ("01_structure_endpoints", "02_structure_change",
                     "03_routeB_distributions", "04_mass_dependence") for ext in ("png", "pdf"))
        if (not args.refresh or not (out / "source_manifest.json").is_file()
                or any(p.name not in known for p in out.iterdir())):
            raise SystemExit("Existing output is not refreshable; choose a new directory")
    records, ledgers, manifest = load_sources(root)
    eligible = [r for r in records if r["eligible_structure"]]
    cells = grouped(eligible, ("family", "N_total", "eta"))
    masses = grouped(eligible, ("family", "N_total", "eta", "M"))
    out.mkdir(parents=True, exist_ok=args.refresh)
    write_csv(out / "trajectories.csv", records)
    write_csv(out / "structural_cells.csv", cells)
    write_csv(out / "structural_by_mass.csv", masses)
    write_csv(out / "source_ledgers.csv", ledgers)
    plots(out, eligible, cells, masses)
    manifest["analysis_script_sha256"] = digest(Path(__file__))
    (out / "source_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    excluded = [r for r in records if not r["eligible_structure"]]
    if excluded:
        write_csv(out / "exclusions.csv", excluded)
    lines = ["# Sound-speed structural screening — completed 2026-09-09", "",
        "Status: analysis checkpoint, not final sound-speed validation.", "",
        f"Read {len(records)} unique trajectories from {len(ledgers)} leaf datasets. "
        f"Excluded {len(excluded)} from structural summaries. Merged copies are not counted again.", "",
        "## Results", "",
        "| Series | N | Cells | Mean-shift flags | Dual-tail screens | Dense low-neighbour flags |",
        "|---|---:|---:|---:|---:|---:|"]
    for family, n in sorted({(r["family"],r["N_total"]) for r in cells}):
        rr=[r for r in cells if r["family"]==family and r["N_total"]==n]
        lines.append(f"| {family} | {n} | {len(rr)} | {sum(r['endpoint_shift_flag'] for r in rr)} | {sum(r['dual_tail_screen'] for r in rr)} | {sum(r['dense_low_neighbors_flag'] for r in rr)} |")
    lines += ["", "## Explicit exclusions", "",
        f"{sum(bool(r['health_reason']) for r in records)} records have nonzero logged health counters. "
        f"{sum(bool(r['missing_structure_reason']) for r in records)} have undefined structure (overlap possible). "
        "Every excluded record, source row, seed and reason is in `exclusions.csv`.", "",
        "| Health reason | Records |", "|---|---:|"]
    for reason, count in sorted(Counter(r['health_reason'] for r in records if r['health_reason']).items()):
        lines.append(f"| {reason} | {count} |")
    lines += ["", "Known overnight clamp-repair exclusion:", ""]
    for r in excluded:
        if 'wall_clamp_repairs=1' not in r['health_reason']:
            continue
        lines.append(f"- {r['family']}, N={r['N_total']}, η={r['eta']}, M={r['M']}, seed={r['seed']}: {r['health_reason'] or r['missing_structure_reason'] or 'incomplete ledger'}; `{r['source']}:{r['source_line']}`.")
    lines += ["", "## Definitions and limitations", "",
        "- Paired differences are computed per trajectory, then averaged. SEM is across repeats; it is not a time-series uncertainty or proof of independent seeds.",
        "- Density-cell summaries pool divider masses with their observed repeat counts. Their SEM is descriptive, not an iid confidence interval for one mass. Use the per-mass table for mass-specific inference.",
        "- Source mean-neighbour counts average over particles with at least one neighbour; the source global order uses all particles in its normalization. No definitions were changed here.",
        "- Absolute mean endpoint change >0.1 is the pre-existing descriptive screen. Per-trajectory large-change fractions are also supplied so opposing changes cannot silently cancel.",
        "- Dual-tail screen: at least 10 samples and at least max(2, ceil(0.1 n)) end values in each tail (<0.4 and >0.8). This requests distribution inspection; it does not establish bimodality. Per-mass tables separate mass mixtures.",
        "- Mean neighbour count <4.3 is flagged only for η≥0.6. It is not proof of square symmetry. ψ4 requires particle configurations; the structural CSVs do not contain them.",
        "- Undefined ψ6 when there are no neighbours is retained as missing in trajectories.csv and excluded from paired structural summaries. This is not a failed simulation. Dilute-gas ψ6 with few neighbours is not an equilibrium/phase diagnostic.",
        "- Complete historical batch ledgers plus no explicit health warning are not a proof that every possible health event was logged. The known clamp warning is excluded regardless of its historical valid flag.",
        "- The structural records do not contain measured release temperature. No temperature normalization, equilibrium certification, FFT extraction or sound-speed refit is performed here.",
        "- Absence of endpoint change does not establish stationarity within the measurement, ergodicity, or equilibrium. No causal claim about the initializer follows from this screen alone.",
        "- Accelerated-core validation datasets are outside the explicit campaign selection. Raw data, old plots, simulation core and NeurAIpil are unchanged.",
        "", "## Figures", "",
        "1. `01_structure_endpoints`: hold/end comparisons across all selected series.",
        "2. `02_structure_change`: paired changes with seed SEM.",
        "3. `03_routeB_distributions`: diagnostic distributions at four densities.",
        "4. `04_mass_dependence`: whether endpoint drift varies with divider mass.",
        "", "All figures are provided as PNG and vector PDF. Tables retain source paths, row numbers, seeds and flags. Source SHA-256 hashes are in `source_manifest.json`.",
        "", "## Next bounded work", "",
        "1. Resolve measured-temperature provenance and join trajectory-level FFT quality flags to these records before final sound-speed refits.",
        "2. Refit fixed-aspect family A and mass-subset sensitivity using explicitly eligible records; retain structural flags and avoid interpreting a three-size fit as proof of a bulk limit.",
        "3. Only if needed for the chosen claim, design the small initialization/hold-duration comparison. No additional simulation was launched here.",
        "4. For Paper 2, reconcile the existing driven-run energy ledger before a larger energy-transfer campaign. Do not treat driven ring-down as equilibrium fluctuation data.", ""]
    (out / "README.md").write_text("\n".join(lines))
    print("\n".join(lines[:len(set((r['family'],r['N_total']) for r in cells))+12]))


if __name__ == "__main__":
    main()
