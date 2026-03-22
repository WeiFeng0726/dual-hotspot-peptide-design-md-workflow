from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_xvg(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for line in path.read_text().splitlines():
        if not line or line[0] in "#@":
            continue
        fields = line.split()
        xs.append(float(fields[0]))
        ys.append(float(fields[1]))
    return xs, ys


def mean_in_window(xs: list[float], ys: list[float], start_ns: float, end_ns: float) -> float:
    window = [y for x, y in zip(xs, ys) if start_ns <= x <= end_ns]
    if not window:
        window = ys
    return sum(window) / len(window)


def moving_average(values: list[float], window: int) -> list[float]:
    if window <= 1 or len(values) <= 2:
        return values[:]
    window = max(1, window)
    if window % 2 == 0:
        window += 1
    radius = window // 2
    smoothed: list[float] = []
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    for idx in range(len(values)):
        start = max(0, idx - radius)
        end = min(len(values), idx + radius + 1)
        smoothed.append((prefix[end] - prefix[start]) / (end - start))
    return smoothed


def build_overlay_plot(
    output_path: Path,
    title: str,
    ylabel: str,
    series: list[tuple[list[float], list[float], str, str]],
) -> None:
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(10, 6))
    for xs, ys, label, color in series:
        axis.plot(xs, ys, label=label, color=color, linewidth=1.6)
    axis.set_xlabel("Time (ns)")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--systems", nargs="+", required=True)
    parser.add_argument("--smooth-window", type=int, default=101)
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # High-contrast palette for multi-series overlays.
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    summary_rows: list[dict[str, str | float]] = []
    raw_series: list[tuple[list[float], list[float], str, str]] = []
    smooth_series: list[tuple[list[float], list[float], str, str]] = []

    for idx, system_arg in enumerate(args.systems):
        run_dir = Path(system_arg)
        system_name = run_dir.parent.name
        analysis_dir = run_dir / "09_interaction_energy"
        lj_x, lj_y = parse_xvg(analysis_dir / "lj_sr_protein_lig.xvg")
        coul_x, coul_y = parse_xvg(analysis_dir / "coul_sr_protein_lig.xvg")

        if lj_x != coul_x:
            raise ValueError(f"Time grids do not match for {system_name}")

        total_y = [lj + coul for lj, coul in zip(lj_y, coul_y)]
        data_tsv = analysis_dir / "protein_ligand_interaction_energy.tsv"
        with data_tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["time_ns", "lj_sr_kj_per_mol", "coul_sr_kj_per_mol", "total_kj_per_mol"])
            for time_ns, lj_val, coul_val, total_val in zip(lj_x, lj_y, coul_y, total_y):
                writer.writerow(
                    [
                        f"{time_ns:.3f}",
                        f"{lj_val:.6f}",
                        f"{coul_val:.6f}",
                        f"{total_val:.6f}",
                    ]
                )

        color = colors[idx % len(colors)]
        raw_series.append((lj_x, total_y, system_name, color))
        smooth_series.append((lj_x, moving_average(total_y, args.smooth_window), system_name, color))

        summary_rows.append(
            {
                "system": system_name,
                "mean_0_50_ns_kj_per_mol": sum(total_y) / len(total_y),
                "mean_30_50_ns_kj_per_mol": mean_in_window(lj_x, total_y, 30.0, 50.0),
                "min_kj_per_mol": min(total_y),
                "max_kj_per_mol": max(total_y),
                "points": len(total_y),
            }
        )

    build_overlay_plot(
        output_dir / "protein_ligand_interaction_energy_overlay_raw.png",
        "Protein-Ligand Interaction Energy Across Recent MD Systems",
        "Protein-Ligand Interaction Energy (kJ/mol)",
        raw_series,
    )
    build_overlay_plot(
        output_dir / "protein_ligand_interaction_energy_overlay_smoothed.png",
        f"Protein-Ligand Interaction Energy Across Recent MD Systems (Smoothed, window={args.smooth_window})",
        "Protein-Ligand Interaction Energy (kJ/mol)",
        smooth_series,
    )
    build_overlay_plot(
        output_dir / "protein_ligand_interaction_energy_overlay.png",
        "Protein-Ligand Interaction Energy Across Recent MD Systems",
        "Protein-Ligand Interaction Energy (kJ/mol)",
        raw_series,
    )

    with (output_dir / "protein_ligand_interaction_energy_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "system",
                "mean_0_50_ns_kj_per_mol",
                "mean_30_50_ns_kj_per_mol",
                "min_kj_per_mol",
                "max_kj_per_mol",
                "points",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
