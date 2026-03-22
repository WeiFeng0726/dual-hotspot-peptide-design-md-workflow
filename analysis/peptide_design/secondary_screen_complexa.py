import csv
import math
import os
import shutil
from collections import defaultdict
from pathlib import Path


SOURCE_CSV = Path(
    os.environ.get(
        "SOURCE_CSV",
        "/mnt/d/Project/cw/ly/Proteina-Complexa-runs/"
        "ly_b404_b405_pep8_12_bulk64_af2_rescore/af2_rescore_ranked.csv",
    )
)
OUTPUT_DIR = Path(
    os.environ.get(
        "OUTPUT_DIR",
        "/mnt/d/Project/cw/ly/Complexa_B404_B405_secondary_screening",
    )
)
FINAL_STRUCT_DIR = OUTPUT_DIR / "final_complexes"
TARGET_INPUT_START = int(os.environ.get("TARGET_INPUT_START", "265"))
HOTSPOT_SOURCE_IDS = [
    int(token)
    for token in os.environ.get("HOTSPOT_SOURCE_IDS", "404,405").split(",")
    if token.strip()
]
HOTSPOT_RENUMBERED = [rid - TARGET_INPUT_START + 1 for rid in HOTSPOT_SOURCE_IDS]
REQUIRED_HOTSPOT_COVERAGE = int(os.environ.get("REQUIRED_HOTSPOT_COVERAGE", str(len(HOTSPOT_SOURCE_IDS))))
MIN_CONTACT_BINDER_RESIDUES = int(os.environ.get("MIN_CONTACT_BINDER_RESIDUES", "2"))
MAX_HOTSPOT_DISTANCE = float(os.environ.get("MAX_HOTSPOT_DISTANCE", "5.0"))
MIN_BINDER_PLDDT = float(os.environ.get("MIN_BINDER_PLDDT", "0.58"))
MIN_I_PTM = float(os.environ.get("MIN_I_PTM", "0.20"))
MAX_MIN_IPAE = float(os.environ.get("MAX_MIN_IPAE", "0.45"))
FINAL_SELECTION_LIMIT = int(os.environ.get("FINAL_SELECTION_LIMIT", "20"))


def parse_pdb_atoms(pdb_path: Path):
    residues = defaultdict(list)
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            chain = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residues[(chain, resseq, icode)].append((atom_name, x, y, z))
    return residues


def residue_distance(res_atoms, hotspot_atoms):
    best = math.inf
    for _, x1, y1, z1 in res_atoms:
        for _, x2, y2, z2 in hotspot_atoms:
            d = math.dist((x1, y1, z1), (x2, y2, z2))
            if d < best:
                best = d
    return best


def choose_chains(residues):
    chain_counts = defaultdict(set)
    for chain, resseq, icode in residues:
        chain_counts[chain].add((resseq, icode))
    ordered = sorted(chain_counts.items(), key=lambda item: len(item[1]), reverse=True)
    target_chain = ordered[0][0]
    binder_chain = ordered[-1][0]
    return target_chain, binder_chain, {chain: len(res_ids) for chain, res_ids in chain_counts.items()}


def compute_hotspot_metrics(pdb_path: Path):
    residues = parse_pdb_atoms(pdb_path)
    target_chain, binder_chain, chain_sizes = choose_chains(residues)
    binder_residues = {
        key: atoms
        for key, atoms in residues.items()
        if key[0] == binder_chain
    }
    hotspot_residues = {
        hotspot_id: residues.get((target_chain, hotspot_id, ""), [])
        for hotspot_id in HOTSPOT_RENUMBERED
    }

    min_dist_any = math.inf
    contact_residues_any = set()
    per_hotspot_contacts = {}
    per_hotspot_min_dist = {}

    for source_hotspot_id, hotspot_id in zip(HOTSPOT_SOURCE_IDS, HOTSPOT_RENUMBERED):
        hotspot_atoms = hotspot_residues.get(hotspot_id, [])
        hotspot_contact_residues = set()
        hotspot_min_dist = math.inf
        if hotspot_atoms:
            for (chain, resseq, icode), res_atoms in binder_residues.items():
                d = residue_distance(res_atoms, hotspot_atoms)
                hotspot_min_dist = min(hotspot_min_dist, d)
                min_dist_any = min(min_dist_any, d)
                if d <= 5.0:
                    hotspot_contact_residues.add(resseq)
                    contact_residues_any.add(resseq)
        per_hotspot_contacts[source_hotspot_id] = len(hotspot_contact_residues)
        per_hotspot_min_dist[source_hotspot_id] = hotspot_min_dist if hotspot_min_dist < math.inf else None

    hotspot_coverage = sum(1 for source_id in HOTSPOT_SOURCE_IDS if per_hotspot_contacts.get(source_id, 0) > 0)
    metrics = {
        "target_chain": target_chain,
        "binder_chain": binder_chain,
        "target_residue_count": chain_sizes.get(target_chain, 0),
        "binder_residue_count": chain_sizes.get(binder_chain, 0),
        "hotspot_coverage": hotspot_coverage,
        "contact_binder_residues": len(contact_residues_any),
        "min_hotspot_distance": min_dist_any if min_dist_any < math.inf else None,
    }
    for source_hotspot_id in HOTSPOT_SOURCE_IDS:
        metrics[f"contacts_hotspot_{source_hotspot_id}"] = per_hotspot_contacts.get(source_hotspot_id, 0)
        metrics[f"min_dist_hotspot_{source_hotspot_id}"] = per_hotspot_min_dist.get(source_hotspot_id)
    return metrics


def as_float(row, key):
    try:
        return float(row[key])
    except Exception:
        return None


def binder_plddt_from_loss(plddt_loss):
    if plddt_loss is None:
        return None
    return 1.0 - plddt_loss


def normalized_pae_to_angstrom(norm_pae):
    if norm_pae is None:
        return None
    return norm_pae * 31.0


def secondary_score(row):
    # Higher is better. Reward is less negative when better.
    reward = as_float(row, "total_reward") or -999.0
    binder_plddt = as_float(row, "binder_plddt") or 0.0
    i_ptm = as_float(row, "i_ptm") or 0.0
    min_ipae = as_float(row, "min_ipae") or 999.0
    hotspot_coverage = int(row["hotspot_coverage"])
    contact_res = int(row["contact_binder_residues"])
    min_dist = as_float(row, "min_hotspot_distance")
    min_dist_bonus = 0.0 if min_dist is None else max(0.0, 6.0 - min_dist) * 0.05
    return (
        reward
        + 0.35 * binder_plddt
        + 0.40 * i_ptm
        - 0.15 * min_ipae
        + 0.20 * hotspot_coverage
        + 0.03 * contact_res
        + min_dist_bonus
    )


def passes_secondary_filter(row):
    return (
        int(row["hotspot_coverage"]) >= REQUIRED_HOTSPOT_COVERAGE
        and int(row["contact_binder_residues"]) >= MIN_CONTACT_BINDER_RESIDUES
        and (as_float(row, "min_hotspot_distance") or 999.0) <= MAX_HOTSPOT_DISTANCE
        and (as_float(row, "binder_plddt") or 0.0) >= MIN_BINDER_PLDDT
        and (as_float(row, "i_ptm") or 0.0) >= MIN_I_PTM
        and (as_float(row, "min_ipae") or 999.0) <= MAX_MIN_IPAE
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if FINAL_STRUCT_DIR.exists():
        shutil.rmtree(FINAL_STRUCT_DIR)
    FINAL_STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    with SOURCE_CSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    enriched = []
    for row in rows:
        pdb_path = Path(row["pdb_path"])
        metrics = compute_hotspot_metrics(pdb_path)
        merged = dict(row)
        merged.update(metrics)
        plddt_loss = as_float(merged, "plddt")
        i_pae = as_float(merged, "i_pae")
        min_ipae = as_float(merged, "min_ipae")
        merged["plddt_loss"] = plddt_loss
        merged["binder_plddt"] = binder_plddt_from_loss(plddt_loss)
        merged["mean_interface_pae_A"] = normalized_pae_to_angstrom(i_pae)
        merged["best_interface_pae_A"] = normalized_pae_to_angstrom(min_ipae)
        merged["secondary_score"] = secondary_score(merged)
        merged["secondary_pass"] = "yes" if passes_secondary_filter(merged) else "no"
        enriched.append(merged)

    enriched_fields = list(enriched[0].keys())
    all_csv = OUTPUT_DIR / "all_ranked_with_hotspot_metrics.csv"
    with all_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=enriched_fields)
        writer.writeheader()
        writer.writerows(enriched)

    passing = [row for row in enriched if row["secondary_pass"] == "yes"]
    passing.sort(key=lambda row: row["secondary_score"], reverse=True)
    for idx, row in enumerate(passing, start=1):
        row["secondary_rank"] = idx

    pass_fields = list(passing[0].keys()) if passing else enriched_fields + ["secondary_rank"]
    shortlist_csv = OUTPUT_DIR / "secondary_shortlist.csv"
    with shortlist_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pass_fields)
        writer.writeheader()
        writer.writerows(passing)

    # Also provide a capped final bundle for convenience.
    final_selection = passing[:FINAL_SELECTION_LIMIT] if len(passing) > FINAL_SELECTION_LIMIT else passing
    for idx, row in enumerate(final_selection, start=1):
        src = Path(row["pdb_path"])
        dst = FINAL_STRUCT_DIR / f"{idx:02d}_{src.name}"
        shutil.copy2(src, dst)
        row["final_pdb_path"] = str(dst)

    final_fields = list(final_selection[0].keys()) if final_selection else pass_fields + ["final_pdb_path"]
    final_csv = OUTPUT_DIR / "final_selection.csv"
    with final_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=final_fields)
        writer.writeheader()
        writer.writerows(final_selection)

    summary_txt = OUTPUT_DIR / "selection_summary.txt"
    with summary_txt.open("w") as handle:
        hotspot_notes = [
            f"source {source_id} -> generated complex residue A{renumbered_id}"
            for source_id, renumbered_id in zip(HOTSPOT_SOURCE_IDS, HOTSPOT_RENUMBERED)
        ]
        handle.write("Target hotspots: " + ", ".join(hotspot_notes) + "\n")
        handle.write("Important note: 'plddt' in the AF2 rescore CSV is a loss term (1 - binder mean pLDDT), not raw pLDDT.\n")
        handle.write("Secondary filter criteria:\n")
        handle.write(f"- hotspot_coverage >= {REQUIRED_HOTSPOT_COVERAGE}\n")
        handle.write(f"- contact_binder_residues >= {MIN_CONTACT_BINDER_RESIDUES}\n")
        handle.write(f"- min_hotspot_distance <= {MAX_HOTSPOT_DISTANCE:.2f} A\n")
        handle.write(f"- binder_plddt >= {MIN_BINDER_PLDDT:.2f}\n")
        handle.write(f"- i_ptm >= {MIN_I_PTM:.2f}\n")
        handle.write(f"- min_ipae <= {MAX_MIN_IPAE:.2f}\n")
        handle.write(f"\nAll candidates scored: {len(enriched)}\n")
        handle.write(f"Secondary shortlist size: {len(passing)}\n")
        handle.write(f"Final copied structures: {len(final_selection)}\n")
        handle.write(f"All ranked table: {all_csv}\n")
        handle.write(f"Shortlist table: {shortlist_csv}\n")
        handle.write(f"Final selection table: {final_csv}\n")
        handle.write(f"Final structures dir: {FINAL_STRUCT_DIR}\n")

    print(f"all_csv={all_csv}")
    print(f"shortlist_count={len(passing)}")
    print(f"final_selection_count={len(final_selection)}")
    print(f"final_csv={final_csv}")
    print(f"final_struct_dir={FINAL_STRUCT_DIR}")


if __name__ == "__main__":
    main()
