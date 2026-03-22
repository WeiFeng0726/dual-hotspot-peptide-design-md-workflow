from __future__ import annotations

import argparse
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--titles", nargs=8, required=True)
    parser.add_argument("--xvgs", nargs=8, required=True)
    parser.add_argument(
        "--suptitle",
        default="Protein-Ligand Complex All-Atom RMSD",
    )
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True, sharey=True)
    colors = [
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#1f78b4",
    ]

    for ax, title, xvg_path, color in zip(axes.flat, args.titles, args.xvgs, colors):
        xs, ys = parse_xvg(Path(xvg_path))
        ax.plot(xs, ys, color=color, linewidth=1.5)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, alpha=0.25)

    for ax in axes[1]:
        ax.set_xlabel("Time (ns)")
    for ax in axes[:, 0]:
        ax.set_ylabel("RMSD (nm)")

    fig.suptitle(args.suptitle, fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    output_png = Path(args.output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=250)
    fig.savefig(output_png.with_suffix(".pdf"))
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
