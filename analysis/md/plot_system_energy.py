from __future__ import annotations

import argparse
import csv
from pathlib import Path

SYSTEM_COLOR_MAP = {
    "s_cdtc-2_c_model": "#1b9e77",
    "s_cdtc-2_2c_model": "#d95f02",
    "s_cdtc-2_3c-2_model": "#7570b3",
    "s_cdtc-2_214a_c_model": "#e7298a",
    "s_cdtc-2_214a_2c_model": "#66a61e",
    "s_cdtc-2_214g_c_model": "#e6ab02",
    "s_cdtc-2_214g_2c_model": "#a6761d",
    "s_cdtc-2_214g_3c-2_model": "#1f78b4",
}

SYSTEM_LABEL_MAP = {
    "s_cdtc-2_214g_3c-2_model": "s_cdtc-2_214g_3c_model",
}


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


def derive_raw_system_name(run_dir: Path) -> str:
    if run_dir.name.startswith("gmx_"):
        return run_dir.parent.name
    return run_dir.name


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
    figure.savefig(output_path.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--systems", nargs="+", required=True)
    parser.add_argument("--smooth-window", type=int, default=101)
    parser.add_argument("--xvg-relpath", default="10_system_energy/system_potential_energy.xvg")
    parser.add_argument("--tsv-filename", default="system_potential_energy.tsv")
    parser.add_argument("--summary-filename", default="system_potential_energy_summary.tsv")
    parser.add_argument("--output-prefix", default="system_potential_energy_overlay")
    parser.add_argument("--title", default="Whole-System Potential Energy Across Recent MD Systems")
    parser.add_argument("--ylabel", default="Potential Energy (kJ/mol)")
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fallback_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    summary_rows: list[dict[str, str | float]] = []
    raw_series: list[tuple[list[float], list[float], str, str]] = []
    smooth_series: list[tuple[list[float], list[float], str, str]] = []

    for idx, system_arg in enumerate(args.systems):
        run_dir = Path(system_arg)
        raw_system_name = derive_raw_system_name(run_dir)
        system_name = SYSTEM_LABEL_MAP.get(raw_system_name, raw_system_name)
        xvg_path = run_dir / args.xvg_relpath
        xs, ys = parse_xvg(xvg_path)

        data_tsv = xvg_path.parent / args.tsv_filename
        with data_tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["time_ps", "potential_kj_per_mol"])
            for time_ps, potential_val in zip(xs, ys):
                writer.writerow([f"{time_ps:.3f}", f"{potential_val:.6f}"])

        xs_ns = [x / 1000.0 for x in xs]
        color = SYSTEM_COLOR_MAP.get(raw_system_name, fallback_colors[idx % len(fallback_colors)])
        raw_series.append((xs_ns, ys, system_name, color))
        smooth_series.append((xs_ns, moving_average(ys, args.smooth_window), system_name, color))

        summary_rows.append(
            {
                "system": system_name,
                "mean_0_50_ns_kj_per_mol": sum(ys) / len(ys),
                "mean_30_50_ns_kj_per_mol": mean_in_window(xs_ns, ys, 30.0, 50.0),
                "min_kj_per_mol": min(ys),
                "max_kj_per_mol": max(ys),
                "points": len(ys),
            }
        )

    build_overlay_plot(
        output_dir / f"{args.output_prefix}_raw.png",
        args.title,
        args.ylabel,
        raw_series,
    )
    build_overlay_plot(
        output_dir / f"{args.output_prefix}_smoothed.png",
        f"{args.title} (Smoothed, window={args.smooth_window})",
        args.ylabel,
        smooth_series,
    )
    build_overlay_plot(
        output_dir / f"{args.output_prefix}.png",
        args.title,
        args.ylabel,
        raw_series,
    )

    with (output_dir / args.summary_filename).open("w", newline="") as handle:
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
