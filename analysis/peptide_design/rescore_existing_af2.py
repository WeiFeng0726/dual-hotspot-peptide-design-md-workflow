import csv
import os
from pathlib import Path

from proteinfoundation.rewards.alphafold2_reward import AF2RewardModel
from proteinfoundation.utils.pdb_utils import get_chain_ids_from_pdb


TARGET_LEN = int(os.environ.get("TARGET_LEN", "156"))
INPUT_ROOT = Path(
    os.environ.get(
        "INPUT_ROOT",
        "/mnt/d/Project/cw/ly/Proteina-Complexa-runs/"
        "ly_b404_b405_pep8_12_bulk64_noreward/inference/"
        "search_binder_local_pipeline_LY_B404_B405_PEP8_12_ly_b404_b405_pep8_12_bulk64_noreward",
    )
)
OUTPUT_ROOT = Path(
    os.environ.get(
        "OUTPUT_ROOT",
        "/mnt/d/Project/cw/ly/Proteina-Complexa-runs/"
        "ly_b404_b405_pep8_12_bulk64_af2_rescore",
    )
)
AF2_DIR = os.environ.get("AF2_DIR", "/home/fatcat/Proteina-Complexa/community_models/ckpts/AF2")


def binder_length_from_name(name: str) -> int | None:
    parts = name.split("_")
    try:
        total_len = int(parts[3])
    except Exception:
        return None
    return total_len - TARGET_LEN


def to_float(value):
    try:
        return float(value)
    except Exception:
        try:
            return float(value.item())
        except Exception:
            return None


def main() -> int:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    progress_csv = OUTPUT_ROOT / "af2_rescore_progress.csv"
    ranked_csv = OUTPUT_ROOT / "af2_rescore_ranked.csv"
    top50_csv = OUTPUT_ROOT / "af2_top50.csv"

    existing = {}
    if progress_csv.exists():
        with progress_csv.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                existing[row["pdb_path"]] = row

    pdb_paths = sorted(INPUT_ROOT.glob("job_*/*.pdb"))
    print(f"found_pdbs={len(pdb_paths)}")
    print(f"resume_existing={len(existing)}")

    fieldnames = [
        "rank",
        "pdb_path",
        "design_name",
        "binder_length",
        "total_reward",
        "i_pae",
        "min_ipae",
        "plddt",
        "ptm",
        "i_ptm",
        "i_con",
        "status",
        "error",
    ]

    model = AF2RewardModel(
        protocol="binder",
        af_params_dir=AF2_DIR,
        reward_weights={
            "con": 0.0,
            "i_pae": -1.0,
            "plddt": 0.0,
            "dgram_cce": 0.0,
            "min_ipae": 0.0,
            "min_ipsae": 0.0,
            "avg_ipsae": 0.0,
            "max_ipsae": 0.0,
            "min_ipsae_10": 0.0,
            "max_ipsae_10": 0.0,
            "avg_ipsae_10": 0.0,
        },
        use_multimer=True,
        num_recycles=3,
        use_initial_guess=True,
        use_initial_atom_pos=False,
        seed=0,
        device_id=0,
    )

    rows = list(existing.values())
    try:
        for idx, pdb_path in enumerate(pdb_paths, start=1):
            pdb_path_str = str(pdb_path)
            if pdb_path_str in existing:
                print(f"[{idx}/{len(pdb_paths)}] skip {pdb_path.name}")
                continue

            design_name = pdb_path.parent.name
            binder_length = binder_length_from_name(design_name)
            target_chain, binder_chain = get_chain_ids_from_pdb(pdb_path_str)

            row = {
                "rank": "",
                "pdb_path": pdb_path_str,
                "design_name": design_name,
                "binder_length": binder_length,
                "total_reward": "",
                "i_pae": "",
                "min_ipae": "",
                "plddt": "",
                "ptm": "",
                "i_ptm": "",
                "i_con": "",
                "status": "ok",
                "error": "",
            }

            try:
                result = model.score(
                    pdb_path_str,
                    requires_grad=False,
                    binder_chain=binder_chain,
                    target_chain=target_chain,
                )
                reward = result["reward"]
                row["total_reward"] = to_float(result["total_reward"])
                row["i_pae"] = to_float(reward.get("i_pae"))
                row["min_ipae"] = to_float(reward.get("min_ipae"))
                row["plddt"] = to_float(reward.get("plddt"))
                row["ptm"] = to_float(result.get("ptm"))
                row["i_ptm"] = to_float(reward.get("i_ptm_log"))
                row["i_con"] = to_float(reward.get("i_con"))
                print(
                    f"[{idx}/{len(pdb_paths)}] ok {design_name} "
                    f"reward={row['total_reward']:.4f} i_pae={row['i_pae']:.4f} plddt={row['plddt']:.4f}"
                )
            except Exception as exc:
                row["status"] = "error"
                row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"[{idx}/{len(pdb_paths)}] error {design_name} {row['error']}")

            existing[pdb_path_str] = row
            rows = list(existing.values())
            with progress_csv.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
    finally:
        model.cleanup()

    ok_rows = [row for row in rows if row["status"] == "ok" and row["total_reward"] not in ("", None)]
    ok_rows.sort(key=lambda row: float(row["total_reward"]), reverse=True)
    for rank, row in enumerate(ok_rows, start=1):
        row["rank"] = rank

    with ranked_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ok_rows)

    with top50_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ok_rows[:50])

    print(f"ok_rows={len(ok_rows)}")
    print(f"ranked_csv={ranked_csv}")
    print(f"top50_csv={top50_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
