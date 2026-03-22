#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(r"D:\Project\cw\ly")
SEED17_SUMMARY = ROOT / "Complexa_B348_secondary_screening_seed17" / "selection_summary.txt"
B348_DRIVER_PID = ROOT / "Proteina-Complexa-runs" / "b348_until_10" / "b348_until_10.pid"
JOINT_WAIT_PID = ROOT / "Joint_peptide_combo_screen" / "joint_screen_wait.pid"
JOINT_OUTPUT = ROOT / "Joint_peptide_combo_screen"
SCRIPT_DIR = Path(__file__).resolve().parent
EVAL_SCRIPT = SCRIPT_DIR / "evaluate_joint_peptide_combos.py"
POLL_SECONDS = 10


def log(msg: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def read_pid(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        return int(path.read_text(encoding="utf-8").strip().splitlines()[0])
    except Exception:
        return None


def terminate_windows_pid(path: Path, label: str) -> None:
    pid = read_pid(path)
    if pid is None:
        log(f"{label}: no pid file or unreadable")
        return
    try:
        os.kill(pid, signal.SIGTERM)
        log(f"{label}: sent SIGTERM to Windows pid {pid}")
    except ProcessLookupError:
        log(f"{label}: pid {pid} not running")
    except Exception as exc:
        log(f"{label}: failed to stop pid {pid}: {exc}")


def wsl_bash(command: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["wsl", "bash", "-lc", command],
        capture_output=True,
        text=True,
        check=False,
    )


def terminate_linux_pidfile(path: Path, label: str) -> None:
    pid = read_pid(path)
    if pid is None:
        log(f"{label}: no linux pid available")
        return
    result = wsl_bash(f"kill {pid} >/dev/null 2>&1 || true")
    if result.returncode == 0:
        log(f"{label}: sent kill to Linux pid {pid}")
    else:
        log(f"{label}: kill returned code {result.returncode}")


def kill_seed23_if_started() -> None:
    script = r"""
ps -eo pid=,args= | grep 'ly_b348_pep8_12_bulk3072_seed23' | grep -v grep | awk '{print $1}' | while read pid; do
  kill "$pid" >/dev/null 2>&1 || true
done
"""
    result = wsl_bash(script)
    if result.returncode == 0:
        log("checked and stopped any seed23 Linux processes if they existed")
    else:
        log(f"seed23 cleanup returned code {result.returncode}")


def as_float(row: dict[str, str], key: str, default: float = -1e9) -> float:
    try:
        return float(row.get(key, ""))
    except Exception:
        return default


def aggregate_b348() -> int:
    screen_dirs = sorted(ROOT.glob("Complexa_B348_secondary_screening*"))
    rows: list[dict[str, str]] = []
    for screen_dir in screen_dirs:
        shortlist = screen_dir / "secondary_shortlist.csv"
        if not shortlist.exists():
            continue
        with shortlist.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                row["screen_dir"] = str(screen_dir)
                rows.append(row)

    unique: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get("pdb_path") or row.get("design_name") or ""
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

    out_dir = ROOT / "Complexa_B348_aggregate_top10"
    final_struct_dir = out_dir / "final_complexes"
    out_dir.mkdir(parents=True, exist_ok=True)
    if final_struct_dir.exists():
        shutil.rmtree(final_struct_dir)
    final_struct_dir.mkdir(parents=True, exist_ok=True)

    top_rows = ranked[:10]
    for idx, row in enumerate(top_rows, start=1):
        src = Path(row["pdb_path"].replace("/mnt/d/", "D:/").replace("/", "\\"))
        dst = final_struct_dir / f"{idx:02d}_{src.name}"
        if src.exists():
            shutil.copy2(src, dst)
            row["final_pdb_path"] = str(dst)
        row["aggregate_rank"] = str(idx)

    if ranked:
        with (out_dir / "all_ranked.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(ranked[0].keys()))
            writer.writeheader()
            writer.writerows(ranked)

    if top_rows:
        with (out_dir / "top10.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(top_rows[0].keys()))
            writer.writeheader()
            writer.writerows(top_rows)

    summary = [
        "B348 aggregate summary",
        f"Loaded rows: {len(rows)}",
        f"Unique candidates: {len(unique)}",
        f"Top N available: {len(top_rows)}",
    ]
    (out_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    log(f"aggregate complete: top_available={len(top_rows)} from_unique={len(unique)}")
    return len(top_rows)


def run_joint_eval() -> int:
    JOINT_OUTPUT.mkdir(parents=True, exist_ok=True)
    cmd = [
        r"D:\software\conda\python.exe",
        str(EVAL_SCRIPT),
        "--no-wait-for-b348",
        "--csv-348",
        str(ROOT / "Complexa_B348_aggregate_top10" / "top10.csv"),
        "--output-dir",
        str(JOINT_OUTPUT),
    ]
    log(f"running joint evaluation: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (JOINT_OUTPUT / "finalize_joint_run.stdout.log").write_text(result.stdout, encoding="utf-8")
    (JOINT_OUTPUT / "finalize_joint_run.stderr.log").write_text(result.stderr, encoding="utf-8")
    log(f"joint evaluation finished with code {result.returncode}")
    return result.returncode


def main() -> int:
    log("waiting for seed17 B348 screening to finish")
    while not SEED17_SUMMARY.exists():
        time.sleep(POLL_SECONDS)

    log(f"detected seed17 completion: {SEED17_SUMMARY}")
    available = aggregate_b348()

    terminate_windows_pid(JOINT_WAIT_PID, "joint_wait")

    if available < 10:
        log(f"available={available} < 10, stopping further B348 expansion and using current set")
        terminate_linux_pidfile(B348_DRIVER_PID, "b348_driver")
        kill_seed23_if_started()
    else:
        log(f"available={available} >= 10, no extra stop action needed beyond finishing current round")
        terminate_linux_pidfile(B348_DRIVER_PID, "b348_driver_post_target")

    return run_joint_eval()


if __name__ == "__main__":
    sys.exit(main())
