#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True)
class AtomRecord:
    record_name: str
    atom_name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    coord: np.ndarray
    occupancy: float
    tempfactor: float
    element: str
    charge: str

    @property
    def is_hydrogen(self) -> bool:
        if self.element:
            return self.element.upper() == "H"
        return self.atom_name.strip().upper().startswith("H")


@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    resseq: int
    icode: str
    resname: str
    atoms: Tuple[AtomRecord, ...]


@dataclass
class ParsedStructure:
    path: Path
    chains: Dict[str, List[ResidueRecord]]


@dataclass
class HotspotMetrics:
    coverage_count: int
    min_overall_distance: float
    min_distance_by_hotspot: Dict[int, float]
    contact_residues_by_hotspot: Dict[int, int]
    total_contact_residues: int


@dataclass
class Candidate:
    group: str
    row: Dict[str, str]
    path: Path
    rank: int
    design_name: str
    sequence: str
    target_chain: str
    peptide_chain: str
    structure: ParsedStructure
    target_residues: List[ResidueRecord]
    peptide_residues: List[ResidueRecord]
    baseline_metrics: HotspotMetrics
    secondary_score: float
    binder_plddt: float
    i_ptm: float
    min_ipae: float
    mean_interface_pae_a: float
    binder_length: int


def as_float(value: str | None, default: float = float("nan")) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except Exception:
        return default


def as_int(value: str | None, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except Exception:
        return default


def resolve_path(path_str: str) -> Path:
    raw = (path_str or "").strip()
    if raw.startswith("/mnt/") and len(raw) > 6:
        drive = raw[5].upper()
        rest = raw[6:].replace("/", "\\")
        return Path(f"{drive}:{rest}")
    return Path(raw)


def candidate_path_from_row(row: Dict[str, str]) -> Path:
    for key in (
        "organized_pdb_path",
        "final_pdb_path_windows",
        "final_pdb_path",
        "priority_pdb_path",
        "pdb_path",
    ):
        if row.get(key):
            return resolve_path(row[key])
    raise ValueError(f"Could not resolve pdb path from row keys: {sorted(row.keys())}")


def parse_pdb(path: Path) -> ParsedStructure:
    chains: Dict[str, List[ResidueRecord]] = {}
    residue_buckets: Dict[Tuple[str, int, str, str], List[AtomRecord]] = {}
    residue_order: List[Tuple[str, int, str, str]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom = AtomRecord(
                record_name=line[0:6].strip() or "ATOM",
                atom_name=line[12:16],
                altloc=line[16].strip(),
                resname=line[17:20].strip(),
                chain=(line[21].strip() or "A"),
                resseq=int(line[22:26].strip()),
                icode=line[26].strip(),
                coord=np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                ),
                occupancy=as_float(line[54:60].strip(), 1.0),
                tempfactor=as_float(line[60:66].strip(), 0.0),
                element=line[76:78].strip(),
                charge=line[78:80].strip(),
            )
            key = (atom.chain, atom.resseq, atom.icode, atom.resname)
            if key not in residue_buckets:
                residue_buckets[key] = []
                residue_order.append(key)
            residue_buckets[key].append(atom)

    for chain, resseq, icode, resname in residue_order:
        residue = ResidueRecord(
            chain=chain,
            resseq=resseq,
            icode=icode,
            resname=resname,
            atoms=tuple(residue_buckets[(chain, resseq, icode, resname)]),
        )
        chains.setdefault(chain, []).append(residue)

    return ParsedStructure(path=path, chains=chains)


def detect_target_and_peptide_chains(structure: ParsedStructure) -> Tuple[str, str]:
    chain_lengths = {chain: len(residues) for chain, residues in structure.chains.items()}
    if len(chain_lengths) < 2:
        raise ValueError(f"{structure.path} does not have at least two chains")
    ranked = sorted(chain_lengths.items(), key=lambda item: item[1], reverse=True)
    return ranked[0][0], ranked[-1][0]


def residue_sequence(residues: Sequence[ResidueRecord]) -> str:
    return "".join(AA3_TO_1.get(res.resname.upper(), "X") for res in residues)


def heavy_atoms(
    residues: Sequence[ResidueRecord],
) -> Tuple[np.ndarray, List[Tuple[int, int, str]], List[AtomRecord]]:
    coords: List[np.ndarray] = []
    residue_refs: List[Tuple[int, int, str]] = []
    atoms_out: List[AtomRecord] = []
    for res_idx, residue in enumerate(residues):
        for atom in residue.atoms:
            if atom.is_hydrogen:
                continue
            coords.append(atom.coord)
            residue_refs.append((res_idx, residue.resseq, residue.icode))
            atoms_out.append(atom)
    if not coords:
        return np.zeros((0, 3), dtype=float), residue_refs, atoms_out
    return np.vstack(coords), residue_refs, atoms_out


def residue_atoms_by_local_index(residues: Sequence[ResidueRecord]) -> Dict[int, np.ndarray]:
    out: Dict[int, np.ndarray] = {}
    for idx, residue in enumerate(residues, start=1):
        atoms = [atom.coord for atom in residue.atoms if not atom.is_hydrogen]
        if atoms:
            out[idx] = np.vstack(atoms)
    return out


def target_ca_map(residues: Sequence[ResidueRecord]) -> Dict[int, np.ndarray]:
    out: Dict[int, np.ndarray] = {}
    for residue in residues:
        for atom in residue.atoms:
            if atom.atom_name.strip() == "CA":
                out[residue.resseq] = atom.coord
                break
    return out


def pairwise_distances(coords_a: np.ndarray, coords_b: np.ndarray) -> np.ndarray:
    if len(coords_a) == 0 or len(coords_b) == 0:
        return np.zeros((len(coords_a), len(coords_b)), dtype=float)
    delta = coords_a[:, None, :] - coords_b[None, :, :]
    return np.linalg.norm(delta, axis=2)


def compute_hotspot_metrics(
    peptide_residues: Sequence[ResidueRecord],
    hotspot_coords: Dict[int, np.ndarray],
    threshold: float = 5.0,
) -> HotspotMetrics:
    pep_coords, pep_res_refs, _ = heavy_atoms(peptide_residues)
    min_distance_by_hotspot: Dict[int, float] = {}
    contact_residues_by_hotspot: Dict[int, int] = {}
    union_contacts: set[Tuple[int, int, str]] = set()
    coverage_count = 0

    for hotspot_id, coords in hotspot_coords.items():
        dist = pairwise_distances(pep_coords, coords)
        if dist.size == 0:
            min_dist = float("inf")
            contact_count = 0
        else:
            min_dist = float(dist.min())
            contact_ids = {pep_res_refs[i] for i, _ in zip(*np.where(dist <= threshold))}
            contact_count = len(contact_ids)
            union_contacts |= contact_ids
        min_distance_by_hotspot[hotspot_id] = min_dist
        contact_residues_by_hotspot[hotspot_id] = contact_count
        if min_dist <= threshold:
            coverage_count += 1

    min_overall_distance = min(min_distance_by_hotspot.values()) if min_distance_by_hotspot else float("inf")
    return HotspotMetrics(
        coverage_count=coverage_count,
        min_overall_distance=min_overall_distance,
        min_distance_by_hotspot=min_distance_by_hotspot,
        contact_residues_by_hotspot=contact_residues_by_hotspot,
        total_contact_residues=len(union_contacts),
    )


def local_hotspot_map(
    source_target_pdb: Path,
    source_chain: str,
    source_hotspots: Iterable[int],
    design_reference_pdb: Path,
) -> Dict[int, int]:
    source_structure = parse_pdb(source_target_pdb)
    source_residues = source_structure.chains[source_chain]
    source_numbers = [res.resseq for res in source_residues]
    source_sequence = residue_sequence(source_residues)

    design_structure = parse_pdb(design_reference_pdb)
    design_target_chain, _ = detect_target_and_peptide_chains(design_structure)
    design_sequence = residue_sequence(design_structure.chains[design_target_chain])

    start = source_sequence.find(design_sequence)
    if start < 0:
        raise ValueError(
            f"Could not locate design target sequence from {design_reference_pdb} inside {source_target_pdb}"
        )
    if source_sequence.find(design_sequence, start + 1) != -1:
        raise ValueError(
            f"Design target sequence from {design_reference_pdb} is not unique inside {source_target_pdb}"
        )

    mapped: Dict[int, int] = {}
    source_hotspots = set(source_hotspots)
    index_lookup = {resseq: idx for idx, resseq in enumerate(source_numbers)}
    for hotspot in source_hotspots:
        if hotspot not in index_lookup:
            raise ValueError(f"Missing hotspot {hotspot} in {source_target_pdb}")
        local_idx = index_lookup[hotspot] - start + 1
        if local_idx < 1 or local_idx > len(design_sequence):
            raise ValueError(
                f"Hotspot {hotspot} falls outside the design target window inferred from {design_reference_pdb}"
            )
        mapped[hotspot] = local_idx
    return mapped


def first_candidate_pdb(csv_path: Path) -> Path:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        first = next(reader, None)
    if first is None:
        raise ValueError(f"{csv_path} is empty")
    return candidate_path_from_row(first)


def build_candidate(row: Dict[str, str], group: str, hotspot_ids: Sequence[int], hotspot_map: Dict[int, int]) -> Candidate:
    pdb_path = candidate_path_from_row(row)
    structure = parse_pdb(pdb_path)
    target_chain, peptide_chain = detect_target_and_peptide_chains(structure)
    target_residues = structure.chains[target_chain]
    peptide_residues = structure.chains[peptide_chain]
    local_hotspots = {hotspot: hotspot_map[hotspot] for hotspot in hotspot_ids}
    target_atoms_by_local_idx = residue_atoms_by_local_index(target_residues)
    hotspot_coords = {hotspot: target_atoms_by_local_idx[local_idx] for hotspot, local_idx in local_hotspots.items()}

    return Candidate(
        group=group,
        row=row,
        path=pdb_path,
        rank=as_int(
            row.get("all_pass_rank")
            or row.get("aggregate_rank")
            or row.get("priority_rank")
            or row.get("rank"),
            0,
        ),
        design_name=row.get("design_name") or pdb_path.stem,
        sequence=row.get("binder_sequence") or residue_sequence(peptide_residues),
        target_chain=target_chain,
        peptide_chain=peptide_chain,
        structure=structure,
        target_residues=target_residues,
        peptide_residues=peptide_residues,
        baseline_metrics=compute_hotspot_metrics(peptide_residues, hotspot_coords),
        secondary_score=as_float(row.get("secondary_score"), 0.0),
        binder_plddt=as_float(row.get("binder_plddt"), 0.0),
        i_ptm=as_float(row.get("i_ptm"), 0.0),
        min_ipae=as_float(row.get("min_ipae"), 1.0),
        mean_interface_pae_a=as_float(row.get("mean_interface_pae_A"), float("nan")),
        binder_length=as_int(row.get("binder_length"), len(peptide_residues)),
    )


def load_candidates(csv_path: Path, group: str, hotspot_ids: Sequence[int], hotspot_map: Dict[int, int]) -> List[Candidate]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [build_candidate(row, group, hotspot_ids, hotspot_map) for row in rows]


def kabsch_align(reference: np.ndarray, mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    if reference.shape != mobile.shape:
        raise ValueError(f"Shape mismatch for alignment: {reference.shape} vs {mobile.shape}")
    if len(reference) < 3:
        raise ValueError("Need at least 3 points for stable alignment")

    ref_centroid = reference.mean(axis=0)
    mob_centroid = mobile.mean(axis=0)
    ref_centered = reference - ref_centroid
    mob_centered = mobile - mob_centroid
    covariance = mob_centered.T @ ref_centered
    v, _, wt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(v @ wt))
    rotation = v @ correction @ wt
    translation = ref_centroid - mob_centroid @ rotation
    transformed = mobile @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - reference) ** 2, axis=1))))
    return rotation, translation, rmsd


def transform_coords(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    if len(coords) == 0:
        return coords.copy()
    return coords @ rotation + translation


def common_target_ca_points(reference: Candidate, mobile: Candidate) -> Tuple[np.ndarray, np.ndarray]:
    ref_map = target_ca_map(reference.target_residues)
    mob_map = target_ca_map(mobile.target_residues)
    shared = sorted(set(ref_map).intersection(mob_map))
    if len(shared) < 3:
        raise ValueError(f"Not enough shared target CA atoms between {reference.design_name} and {mobile.design_name}")
    return np.vstack([ref_map[idx] for idx in shared]), np.vstack([mob_map[idx] for idx in shared])


def peptide_peptide_metrics(
    reference_peptide: Sequence[ResidueRecord],
    mobile_peptide: Sequence[ResidueRecord],
    rotation: np.ndarray,
    translation: np.ndarray,
    hard_cutoff: float,
    soft_cutoff: float,
    contact_cutoff: float,
) -> Dict[str, float]:
    ref_coords, _, _ = heavy_atoms(reference_peptide)
    mob_coords, _, _ = heavy_atoms(mobile_peptide)
    mob_transformed = transform_coords(mob_coords, rotation, translation)
    dist = pairwise_distances(ref_coords, mob_transformed)
    if dist.size == 0:
        min_dist = float("inf")
        hard = 0
        soft = 0
        contacts = 0
    else:
        min_dist = float(dist.min())
        hard = int(np.count_nonzero(dist < hard_cutoff))
        soft = int(np.count_nonzero((dist >= hard_cutoff) & (dist < soft_cutoff)))
        contacts = int(np.count_nonzero(dist < contact_cutoff))
    ref_com = ref_coords.mean(axis=0) if len(ref_coords) else np.zeros(3)
    mob_com = mob_transformed.mean(axis=0) if len(mob_transformed) else np.zeros(3)
    com_dist = float(np.linalg.norm(ref_com - mob_com)) if len(ref_coords) and len(mob_transformed) else float("inf")
    return {
        "peptide_peptide_min_heavy_distance_A": min_dist,
        "peptide_peptide_hard_clash_count_lt_2A": hard,
        "peptide_peptide_soft_clash_count_2_to_2_5A": soft,
        "peptide_peptide_contact_pairs_lt_4_5A": contacts,
        "peptide_peptide_com_distance_A": com_dist,
    }


def transformed_hotspot_metrics(
    mobile_candidate: Candidate,
    target_atoms_by_local_idx: Dict[int, np.ndarray],
    hotspot_ids: Sequence[int],
    hotspot_map: Dict[int, int],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> HotspotMetrics:
    transformed_residues: List[ResidueRecord] = []
    for residue in mobile_candidate.peptide_residues:
        new_atoms = tuple(
            AtomRecord(
                record_name=atom.record_name,
                atom_name=atom.atom_name,
                altloc=atom.altloc,
                resname=atom.resname,
                chain=atom.chain,
                resseq=atom.resseq,
                icode=atom.icode,
                coord=transform_coords(atom.coord.reshape(1, 3), rotation, translation)[0],
                occupancy=atom.occupancy,
                tempfactor=atom.tempfactor,
                element=atom.element,
                charge=atom.charge,
            )
            for atom in residue.atoms
        )
        transformed_residues.append(
            ResidueRecord(
                chain=residue.chain,
                resseq=residue.resseq,
                icode=residue.icode,
                resname=residue.resname,
                atoms=new_atoms,
            )
        )

    hotspot_coords = {
        hotspot: target_atoms_by_local_idx[hotspot_map[hotspot]]
        for hotspot in hotspot_ids
    }
    return compute_hotspot_metrics(transformed_residues, hotspot_coords)


def cross_hotspot_min_distance(
    peptide_residues: Sequence[ResidueRecord],
    target_atoms_by_local_idx: Dict[int, np.ndarray],
    hotspot_ids: Sequence[int],
    hotspot_map: Dict[int, int],
    rotation: np.ndarray | None = None,
    translation: np.ndarray | None = None,
) -> float:
    coords, _, _ = heavy_atoms(peptide_residues)
    if rotation is not None and translation is not None:
        coords = transform_coords(coords, rotation, translation)
    hotspot_arrays = [target_atoms_by_local_idx[hotspot_map[hotspot]] for hotspot in hotspot_ids]
    all_hotspot_coords = np.vstack(hotspot_arrays)
    dist = pairwise_distances(coords, all_hotspot_coords)
    return float(dist.min()) if dist.size else float("inf")


def compatibility_score(row: Dict[str, float]) -> float:
    base_terms = [
        row["secondary_score_404405"],
        row["binder_plddt_404405"],
        row["i_ptm_404405"],
        row["secondary_score_348"],
        row["binder_plddt_348"],
        row["i_ptm_348"],
        row["hotspot_coverage_404405"] / 2.0,
        1.0 if row["hotspot_coverage_348"] >= 1 else 0.0,
    ]
    base = float(np.mean(base_terms))
    penalty = (
        0.15 * row["peptide_peptide_hard_clash_count_lt_2A"]
        + 0.02 * row["peptide_peptide_soft_clash_count_2_to_2_5A"]
        + 0.002 * row["peptide_peptide_contact_pairs_lt_4_5A"]
        + 0.05 * max(0.0, row["delta_min_hotspot_distance_348_A"])
        + 0.05 * max(0.0, row["delta_min_hotspot_distance_404405_A"])
        + 0.10 * max(0.0, 2.5 - row["peptide_peptide_min_heavy_distance_A"])
    )
    return base - penalty


def recommendation_label(row: Dict[str, float]) -> str:
    recommended = (
        row["target_alignment_ca_rmsd_A"] <= 1.5
        and row["peptide_peptide_hard_clash_count_lt_2A"] == 0
        and row["peptide_peptide_soft_clash_count_2_to_2_5A"] <= 5
        and row["peptide_peptide_min_heavy_distance_A"] >= 2.5
        and row["hotspot_coverage_404405"] == 2
        and row["hotspot_coverage_348"] >= 1
        and row["cross_404405_peptide_to_348_hotspot_min_A"] > 5.0
        and row["cross_348_peptide_to_404405_hotspot_min_A"] > 5.0
        and row["delta_min_hotspot_distance_348_A"] <= 1.5
        and row["delta_min_hotspot_distance_404405_A"] <= 1.5
    )
    if recommended:
        return "recommended"

    borderline = (
        row["peptide_peptide_hard_clash_count_lt_2A"] == 0
        and row["hotspot_coverage_404405"] == 2
        and row["hotspot_coverage_348"] >= 1
        and row["cross_404405_peptide_to_348_hotspot_min_A"] > 4.0
        and row["cross_348_peptide_to_404405_hotspot_min_A"] > 4.0
    )
    if borderline:
        return "borderline"
    return "reject"


def priority_tier(row: Dict[str, float]) -> str:
    if row["joint_use_recommendation"] == "reject":
        return "low"
    if (
        row["compatibility_score"] >= 0.73
        and row["binder_plddt_348"] >= 0.60
        and row["i_ptm_348"] >= 0.58
        and row["binder_plddt_404405"] >= 0.66
        and row["i_ptm_404405"] >= 0.73
    ):
        return "high"
    if (
        row["compatibility_score"] >= 0.68
        and row["binder_plddt_348"] >= 0.56
        and row["i_ptm_348"] >= 0.52
    ):
        return "medium"
    return "low"


def write_pdb_atom(serial: int, atom: AtomRecord, chain: str, coord: np.ndarray) -> str:
    element = (atom.element or atom.atom_name.strip()[:1]).upper()
    return (
        f"{atom.record_name:<6}{serial:>5} "
        f"{atom.atom_name:<4}{atom.altloc:1}{atom.resname:>3} {chain:1}"
        f"{atom.resseq:>4}{atom.icode:1}   "
        f"{coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}"
        f"{atom.occupancy:>6.2f}{atom.tempfactor:>6.2f}          "
        f"{element:>2}{atom.charge:>2}"
    )


def write_combo_pdb(
    out_path: Path,
    target_candidate: Candidate,
    pep_404_candidate: Candidate,
    pep_348_candidate: Candidate,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> None:
    serial = 1
    lines = ["MODEL     1"]

    for residue in target_candidate.target_residues:
        for atom in residue.atoms:
            lines.append(write_pdb_atom(serial, atom, "A", atom.coord))
            serial += 1
    lines.append(
        f"TER   {serial:>5}      {target_candidate.target_residues[-1].resname:>3} A{target_candidate.target_residues[-1].resseq:>4}"
    )
    serial += 1

    for residue in pep_404_candidate.peptide_residues:
        for atom in residue.atoms:
            lines.append(write_pdb_atom(serial, atom, "B", atom.coord))
            serial += 1
    lines.append(
        f"TER   {serial:>5}      {pep_404_candidate.peptide_residues[-1].resname:>3} B{pep_404_candidate.peptide_residues[-1].resseq:>4}"
    )
    serial += 1

    for residue in pep_348_candidate.peptide_residues:
        for atom in residue.atoms:
            coord = transform_coords(atom.coord.reshape(1, 3), rotation, translation)[0]
            lines.append(write_pdb_atom(serial, atom, "C", coord))
            serial += 1
    lines.append(
        f"TER   {serial:>5}      {pep_348_candidate.peptide_residues[-1].resname:>3} C{pep_348_candidate.peptide_residues[-1].resseq:>4}"
    )
    lines.extend(["ENDMDL", "END"])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_row_count(csv_path: Path) -> int:
    if not csv_path.exists():
        return 0
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def wait_for_row_target(csv_path: Path, expected_rows: int, poll_seconds: int, timeout_hours: float) -> int:
    deadline = time.time() + timeout_hours * 3600.0
    while True:
        rows = load_row_count(csv_path)
        if rows >= expected_rows:
            return rows
        if time.time() > deadline:
            raise TimeoutError(
                f"Timed out waiting for {csv_path} to reach {expected_rows} rows; current={rows}"
            )
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[wait] {timestamp} current_rows={rows}, target={expected_rows}")
        time.sleep(poll_seconds)


def evaluate_combos(
    candidates_404: Sequence[Candidate],
    candidates_348: Sequence[Candidate],
    hotspot_map: Dict[int, int],
) -> List[Dict[str, float | str]]:
    rows: List[Dict[str, float | str]] = []
    for cand_404 in candidates_404:
        target_atoms_by_local_idx = residue_atoms_by_local_index(cand_404.target_residues)
        for cand_348 in candidates_348:
            ref_ca, mob_ca = common_target_ca_points(cand_404, cand_348)
            rotation, translation, target_rmsd = kabsch_align(ref_ca, mob_ca)

            pep_metrics = peptide_peptide_metrics(
                cand_404.peptide_residues,
                cand_348.peptide_residues,
                rotation=rotation,
                translation=translation,
                hard_cutoff=2.0,
                soft_cutoff=2.5,
                contact_cutoff=4.5,
            )

            metrics_404 = cand_404.baseline_metrics
            metrics_348 = transformed_hotspot_metrics(
                mobile_candidate=cand_348,
                target_atoms_by_local_idx=target_atoms_by_local_idx,
                hotspot_ids=[348],
                hotspot_map=hotspot_map,
                rotation=rotation,
                translation=translation,
            )

            cross_404_to_348 = cross_hotspot_min_distance(
                cand_404.peptide_residues,
                target_atoms_by_local_idx,
                hotspot_ids=[348],
                hotspot_map=hotspot_map,
            )
            cross_348_to_404405 = cross_hotspot_min_distance(
                cand_348.peptide_residues,
                target_atoms_by_local_idx,
                hotspot_ids=[404, 405],
                hotspot_map=hotspot_map,
                rotation=rotation,
                translation=translation,
            )

            row: Dict[str, float | str] = {
                "combo_id": f"{cand_404.design_name}__{cand_348.design_name}",
                "design_404405": cand_404.design_name,
                "design_348": cand_348.design_name,
                "rank_404405": cand_404.rank,
                "rank_348": cand_348.rank,
                "sequence_404405": cand_404.sequence,
                "sequence_348": cand_348.sequence,
                "binder_length_404405": cand_404.binder_length,
                "binder_length_348": cand_348.binder_length,
                "secondary_score_404405": cand_404.secondary_score,
                "binder_plddt_404405": cand_404.binder_plddt,
                "i_ptm_404405": cand_404.i_ptm,
                "min_ipae_404405": cand_404.min_ipae,
                "mean_interface_pae_A_404405": cand_404.mean_interface_pae_a,
                "secondary_score_348": cand_348.secondary_score,
                "binder_plddt_348": cand_348.binder_plddt,
                "i_ptm_348": cand_348.i_ptm,
                "min_ipae_348": cand_348.min_ipae,
                "mean_interface_pae_A_348": cand_348.mean_interface_pae_a,
                "target_alignment_ca_rmsd_A": target_rmsd,
                "hotspot_coverage_404405": metrics_404.coverage_count,
                "hotspot_coverage_348": metrics_348.coverage_count,
                "min_hotspot_distance_404405_A": metrics_404.min_overall_distance,
                "min_hotspot_distance_348_A": metrics_348.min_overall_distance,
                "delta_min_hotspot_distance_404405_A": metrics_404.min_overall_distance - cand_404.baseline_metrics.min_overall_distance,
                "delta_min_hotspot_distance_348_A": metrics_348.min_overall_distance - cand_348.baseline_metrics.min_overall_distance,
                "contact_residues_404405_total": metrics_404.total_contact_residues,
                "contact_residues_348_total": metrics_348.total_contact_residues,
                "min_dist_404pep_to_404_A": metrics_404.min_distance_by_hotspot[404],
                "min_dist_404pep_to_405_A": metrics_404.min_distance_by_hotspot[405],
                "min_dist_348pep_to_348_A": metrics_348.min_distance_by_hotspot[348],
                "contacts_404pep_to_404": metrics_404.contact_residues_by_hotspot[404],
                "contacts_404pep_to_405": metrics_404.contact_residues_by_hotspot[405],
                "contacts_348pep_to_348": metrics_348.contact_residues_by_hotspot[348],
                "cross_404405_peptide_to_348_hotspot_min_A": cross_404_to_348,
                "cross_348_peptide_to_404405_hotspot_min_A": cross_348_to_404405,
                "pdb_404405": str(cand_404.path),
                "pdb_348": str(cand_348.path),
            }
            row.update(pep_metrics)
            row["compatibility_score"] = compatibility_score(row)  # type: ignore[arg-type]
            row["joint_use_recommendation"] = recommendation_label(row)  # type: ignore[arg-type]
            row["joint_priority_tier"] = priority_tier(row)  # type: ignore[arg-type]
            rows.append(row)

    rows.sort(
        key=lambda r: (
            {"recommended": 0, "borderline": 1, "reject": 2}[str(r["joint_use_recommendation"])],
            {"high": 0, "medium": 1, "low": 2}[str(r["joint_priority_tier"])],
            -float(r["compatibility_score"]),
            float(r["peptide_peptide_hard_clash_count_lt_2A"]),
            float(r["peptide_peptide_soft_clash_count_2_to_2_5A"]),
        )
    )
    for idx, row in enumerate(rows, start=1):
        row["combo_rank"] = idx
    return rows


def write_csv(path: Path, rows: Sequence[Dict[str, float | str]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summary_text(
    output_dir: Path,
    hotspot_map: Dict[int, int],
    rows: Sequence[Dict[str, float | str]],
    csv_404: Path,
    csv_348: Path,
    wrote_complexes: int,
) -> str:
    counts = {"recommended": 0, "borderline": 0, "reject": 0}
    tiers = {"high": 0, "medium": 0, "low": 0}
    for row in rows:
        counts[str(row["joint_use_recommendation"])] += 1
        tiers[str(row["joint_priority_tier"])] += 1

    lines = [
        "Joint peptide compatibility screen",
        f"404/405 source csv: {csv_404}",
        f"348 source csv: {csv_348}",
        f"Output dir: {output_dir}",
        "",
        "Static evaluation logic:",
        "- align the 348-complex target onto the 404/405-complex target via target-chain CA atoms",
        "- keep the 404/405 peptide fixed and transform the 348 peptide into the shared target frame",
        "- measure peptide-peptide clashes, hotspot retention, and cross-hotspot interference without MD",
        "",
        "Target hotspot mapping (original -> design local index):",
        f"- 348 -> {hotspot_map[348]}",
        f"- 404 -> {hotspot_map[404]}",
        f"- 405 -> {hotspot_map[405]}",
        "",
        "Recommendation thresholds:",
        "- recommended: no hard clashes, <=5 soft clashes, min heavy-atom peptide distance >=2.5 A, 404/405 retained, 348 retained, cross-hotspot distances >5 A, target alignment RMSD <=1.5 A",
        "- borderline: no hard clashes and both hotspot zones still look accessible, but geometry is less clean",
        "- reject: clashes or obvious cross-hotspot interference",
        "- priority_tier: combines static compatibility with the single-peptide confidence of both members; use this for which combos to test first",
        "",
        "Counts:",
        f"- total_combos: {len(rows)}",
        f"- recommended: {counts['recommended']}",
        f"- borderline: {counts['borderline']}",
        f"- reject: {counts['reject']}",
        f"- high_priority: {tiers['high']}",
        f"- medium_priority: {tiers['medium']}",
        f"- low_priority: {tiers['low']}",
        f"- merged_complexes_written: {wrote_complexes}",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate static compatibility between 404/405 peptide candidates and B348 peptide candidates."
    )
    parser.add_argument(
        "--csv-404405",
        default=r"D:\Project\cw\ly\Complexa_B404_B405_all_pass\all_pass_sequences.csv",
        help="CSV containing 404/405 candidates.",
    )
    parser.add_argument(
        "--csv-348",
        default=r"D:\Project\cw\ly\Complexa_B348_aggregate_top10\top10.csv",
        help="CSV containing current B348 aggregate candidates.",
    )
    parser.add_argument(
        "--source-target-pdb",
        default=r"D:\Project\cw\ly\BindCraft\chainB.pdb",
        help="Original target pdb used to map original hotspot numbering onto design-local indices.",
    )
    parser.add_argument(
        "--output-dir",
        default=r"D:\Project\cw\ly\Joint_peptide_combo_screen",
        help="Output directory for combo tables and merged structures.",
    )
    parser.add_argument(
        "--required-b348-count",
        type=int,
        default=10,
        help="Number of B348 aggregate candidates to wait for before evaluating.",
    )
    parser.add_argument(
        "--wait-for-b348",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Wait until the B348 aggregate csv reaches the required row count.",
    )
    parser.add_argument(
        "--poll-seconds",
        type=int,
        default=300,
        help="Polling interval while waiting for B348 to finish.",
    )
    parser.add_argument(
        "--timeout-hours",
        type=float,
        default=48.0,
        help="Maximum wait time before giving up.",
    )
    parser.add_argument(
        "--write-top-complexes",
        type=int,
        default=30,
        help="Write merged combo pdbs for the top N scored non-reject combos.",
    )
    args = parser.parse_args()

    csv_404 = Path(args.csv_404405)
    csv_348 = Path(args.csv_348)
    source_target_pdb = Path(args.source_target_pdb)
    output_dir = Path(args.output_dir)

    if args.wait_for_b348:
        rows_ready = wait_for_row_target(
            csv_path=csv_348,
            expected_rows=args.required_b348_count,
            poll_seconds=args.poll_seconds,
            timeout_hours=args.timeout_hours,
        )
        print(f"[wait] B348 ready with {rows_ready} rows")

    output_dir.mkdir(parents=True, exist_ok=True)
    merged_dir = output_dir / "merged_complexes"
    if merged_dir.exists():
        shutil.rmtree(merged_dir)
    merged_dir.mkdir(parents=True, exist_ok=True)

    design_reference_pdb = first_candidate_pdb(csv_404)
    hotspot_map = local_hotspot_map(
        source_target_pdb,
        source_chain="B",
        source_hotspots=[348, 404, 405],
        design_reference_pdb=design_reference_pdb,
    )
    candidates_404 = load_candidates(csv_404, group="404405", hotspot_ids=[404, 405], hotspot_map=hotspot_map)
    candidates_348 = load_candidates(csv_348, group="348", hotspot_ids=[348], hotspot_map=hotspot_map)

    rows = evaluate_combos(candidates_404, candidates_348, hotspot_map=hotspot_map)
    recommended = [row for row in rows if row["joint_use_recommendation"] == "recommended"]
    borderline = [row for row in rows if row["joint_use_recommendation"] == "borderline"]
    high_priority = [row for row in rows if row["joint_priority_tier"] == "high"]
    medium_priority = [row for row in rows if row["joint_priority_tier"] == "medium"]

    lookup_404 = {cand.design_name: cand for cand in candidates_404}
    lookup_348 = {cand.design_name: cand for cand in candidates_348}
    to_write = [row for row in rows if row["joint_use_recommendation"] != "reject"][: args.write_top_complexes]
    for row in to_write:
        cand_404 = lookup_404[str(row["design_404405"])]
        cand_348 = lookup_348[str(row["design_348"])]
        ref_ca, mob_ca = common_target_ca_points(cand_404, cand_348)
        rotation, translation, _ = kabsch_align(ref_ca, mob_ca)
        out_path = merged_dir / f"{int(row['combo_rank']):02d}_{cand_404.design_name}__{cand_348.design_name}.pdb"
        write_combo_pdb(out_path, cand_404, cand_404, cand_348, rotation, translation)
        row["merged_combo_pdb"] = str(out_path)

    write_csv(output_dir / "all_combo_metrics.csv", rows)
    write_csv(output_dir / "recommended_combos.csv", recommended)
    write_csv(output_dir / "borderline_combos.csv", borderline)
    (output_dir / "screen_summary.txt").write_text(
        summary_text(
            output_dir=output_dir,
            hotspot_map=hotspot_map,
            rows=rows,
            csv_404=csv_404,
            csv_348=csv_348,
            wrote_complexes=len(to_write),
        ),
        encoding="utf-8",
    )

    print(f"404405_candidates={len(candidates_404)}")
    print(f"b348_candidates={len(candidates_348)}")
    print(f"total_combos={len(rows)}")
    print(f"recommended={len(recommended)}")
    print(f"borderline={len(borderline)}")
    print(f"high_priority={len(high_priority)}")
    print(f"medium_priority={len(medium_priority)}")
    print(f"output_dir={output_dir}")


if __name__ == "__main__":
    main()
