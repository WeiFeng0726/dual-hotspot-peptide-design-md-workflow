import csv
import os
import shutil
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/mnt/d/Project/cw/ly"))
PATTERN = os.environ.get("PATTERN", "Complexa_B348_secondary_screening*")
OUTPUT_DIR = Path(os.environ.get("OUTPUT_DIR", str(ROOT / "Complexa_B348_aggregate_top10")))
TOP_N = int(os.environ.get("TOP_N", "10"))


def as_float(row, key, default=-1e9):
    try:
        return float(row[key])
    except Exception:
        return default


def load_rows():
    rows = []
    for screen_dir in sorted(ROOT.glob(PATTERN)):
        shortlist = screen_dir / "secondary_shortlist.csv"
        if not shortlist.exists():
            continue
        with shortlist.open(newline="") as handle:
            for row in csv.DictReader(handle):
                row["screen_dir"] = str(screen_dir)
                rows.append(row)
    return rows


def main():
    rows = load_rows()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    final_struct_dir = OUTPUT_DIR / "final_complexes"
    if final_struct_dir.exists():
        shutil.rmtree(final_struct_dir)
    final_struct_dir.mkdir(parents=True, exist_ok=True)

    unique = {}
    for row in rows:
        key = row.get("pdb_path") or row.get("design_name")
        prev = unique.get(key)
        if prev is None or as_float(row, "secondary_score") > as_float(prev, "secondary_score"):
            unique[key] = row

    ranked = sorted(
        unique.values(),
        key=lambda row: (
            as_float(row, "secondary_score"),
            as_float(row, "binder_plddt"),
            as_float(row, "i_ptm"),
        ),
        reverse=True,
    )

    top_rows = ranked[:TOP_N]
    for idx, row in enumerate(top_rows, start=1):
        src = Path(row["pdb_path"])
        dst = final_struct_dir / f"{idx:02d}_{src.name}"
        if src.exists():
            shutil.copy2(src, dst)
            row["final_pdb_path"] = str(dst)
        row["aggregate_rank"] = idx

    all_csv = OUTPUT_DIR / "all_ranked.csv"
    if ranked:
        with all_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(ranked[0].keys()))
            writer.writeheader()
            writer.writerows(ranked)

    final_csv = OUTPUT_DIR / "top10.csv"
    if top_rows:
        with final_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(top_rows[0].keys()))
            writer.writeheader()
            writer.writerows(top_rows)

    summary = OUTPUT_DIR / "summary.txt"
    with summary.open("w") as handle:
        handle.write(f"Pattern: {PATTERN}\n")
        handle.write(f"Total shortlisted rows loaded: {len(rows)}\n")
        handle.write(f"Unique shortlisted candidates: {len(ranked)}\n")
        handle.write(f"Top N requested: {TOP_N}\n")
        handle.write(f"Top N available: {len(top_rows)}\n")
        handle.write(f"Output dir: {OUTPUT_DIR}\n")

    print(f"loaded_rows={len(rows)}")
    print(f"unique_candidates={len(ranked)}")
    print(f"top_n_available={len(top_rows)}")
    print(f"final_csv={final_csv}")


if __name__ == "__main__":
    main()
