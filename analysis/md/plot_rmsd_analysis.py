from __future__ import annotations

import math
import statistics
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def read_xvg(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for line in path.read_text().splitlines():
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        xs.append(float(parts[0]))
        ys.append(float(parts[1]))
    return xs, ys


def linear_slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    num = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    den = sum((x - x_mean) ** 2 for x in xs)
    if den == 0:
        return 0.0
    return num / den


def window_stats(xs: list[float], ys: list[float], start_ns: float) -> dict[str, float]:
    pairs = [(x, y) for x, y in zip(xs, ys) if x >= start_ns]
    if not pairs:
        pairs = list(zip(xs, ys))
    wx = [x for x, _ in pairs]
    wy = [y for _, y in pairs]
    return {
        "start_ns": wx[0],
        "end_ns": wx[-1],
        "mean": statistics.fmean(wy),
        "std": statistics.pstdev(wy) if len(wy) > 1 else 0.0,
        "min": min(wy),
        "max": max(wy),
        "slope_nm_per_ns": linear_slope(wx, wy),
    }


def convergence_label(stats_20: dict[str, float], stats_10: dict[str, float], kind: str) -> str:
    slope_abs = abs(stats_20["slope_nm_per_ns"])
    std_10 = stats_10["std"]
    span_10 = stats_10["max"] - stats_10["min"]

    if kind == "protein_backbone":
        if slope_abs < 0.002 and std_10 < 0.03 and span_10 < 0.12:
            return "Protein backbone RMSD is broadly plateaued in the last 10-20 ns."
        if slope_abs < 0.005 and std_10 < 0.05:
            return "Protein backbone RMSD looks near-stationary, with moderate fluctuations."
        return "Protein backbone RMSD still shows noticeable drift or broad fluctuations."

    if slope_abs < 0.003 and std_10 < 0.05 and span_10 < 0.20:
        return "Ligand RMSD is broadly stable relative to the protein backbone."
    if slope_abs < 0.008 and std_10 < 0.10:
        return "Ligand RMSD is moderately stable but still samples multiple poses."
    return "Ligand RMSD remains quite mobile relative to the protein backbone."


def main() -> int:
    if len(sys.argv) not in {2, 3}:
        print("Usage: python plot_rmsd_analysis.py <md_dir> [prefix]")
        return 1

    md_dir = Path(sys.argv[1])
    prefix = sys.argv[2] if len(sys.argv) == 3 else ""
    protein_xvg = md_dir / f"rmsd_backbone{prefix}.xvg"
    ligand_xvg = md_dir / f"rmsd_ligand_heavy_fit_backbone{prefix}.xvg"

    px, py = read_xvg(protein_xvg)
    lx, ly = read_xvg(ligand_xvg)

    p20 = window_stats(px, py, 30.0)
    p10 = window_stats(px, py, 40.0)
    l20 = window_stats(lx, ly, 30.0)
    l10 = window_stats(lx, ly, 40.0)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    axes[0].plot(px, py, color="#355070", linewidth=1.6)
    axes[0].axvspan(30, 50, color="#355070", alpha=0.08)
    axes[0].axhline(p20["mean"], color="#355070", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("RMSD (nm)")
    axes[0].set_title("Protein Backbone RMSD")
    axes[0].grid(True, alpha=0.25)
    axes[0].text(
        0.02,
        0.95,
        f"30-50 ns mean {p20['mean']:.3f} nm\n40-50 ns sd {p10['std']:.3f} nm\n30-50 ns slope {p20['slope_nm_per_ns']:.4f} nm/ns",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "#355070"},
    )

    axes[1].plot(lx, ly, color="#b56576", linewidth=1.6)
    axes[1].axvspan(30, 50, color="#b56576", alpha=0.08)
    axes[1].axhline(l20["mean"], color="#b56576", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("Time (ns)")
    axes[1].set_ylabel("RMSD (nm)")
    axes[1].set_title("Ligand Heavy-Atom RMSD (Fit on Protein Backbone)")
    axes[1].grid(True, alpha=0.25)
    axes[1].text(
        0.02,
        0.95,
        f"30-50 ns mean {l20['mean']:.3f} nm\n40-50 ns sd {l10['std']:.3f} nm\n30-50 ns slope {l20['slope_nm_per_ns']:.4f} nm/ns",
        transform=axes[1].transAxes,
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "#b56576"},
    )

    fig.tight_layout()
    png_name = "rmsd_overview.png" if not prefix else f"rmsd_overview{prefix}.png"
    summary_name = "rmsd_summary.txt" if not prefix else f"rmsd_summary{prefix}.txt"
    fig.savefig(md_dir / png_name, dpi=200)
    plt.close(fig)

    summary_lines = [
        "RMSD analysis summary",
        "",
        f"Protein backbone, 30-50 ns: mean={p20['mean']:.4f} nm, std={p20['std']:.4f} nm, min={p20['min']:.4f} nm, max={p20['max']:.4f} nm, slope={p20['slope_nm_per_ns']:.5f} nm/ns",
        f"Protein backbone, 40-50 ns: mean={p10['mean']:.4f} nm, std={p10['std']:.4f} nm, min={p10['min']:.4f} nm, max={p10['max']:.4f} nm, slope={p10['slope_nm_per_ns']:.5f} nm/ns",
        f"Ligand heavy atom, 30-50 ns: mean={l20['mean']:.4f} nm, std={l20['std']:.4f} nm, min={l20['min']:.4f} nm, max={l20['max']:.4f} nm, slope={l20['slope_nm_per_ns']:.5f} nm/ns",
        f"Ligand heavy atom, 40-50 ns: mean={l10['mean']:.4f} nm, std={l10['std']:.4f} nm, min={l10['min']:.4f} nm, max={l10['max']:.4f} nm, slope={l10['slope_nm_per_ns']:.5f} nm/ns",
        "",
        convergence_label(p20, p10, "protein_backbone"),
        convergence_label(l20, l10, "ligand"),
    ]
    (md_dir / summary_name).write_text("\n".join(summary_lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
