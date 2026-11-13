#!/home/qiuyong/.conda/envs/Agent/bin/python3
"""Shared helpers for the LiFSI diffusion workflow."""

from __future__ import annotations

import json
import os
import subprocess
import time
import getpass
from pathlib import Path
from typing import Any

FINAL_STATUSES = {"diffusive", "max_chunks_reached", "failed"}
RECOVERABLE_FAILURE_REASONS = {"resubmit_failed"}


def workflow_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(root: Path | None = None) -> dict:
    base = root or workflow_root()
    return json.loads((base / "workflow_config.json").read_text(encoding="ascii"))


def executable_path(key: str, root: Path | None = None, fallback: str | None = None) -> str:
    config = load_config(root or workflow_root())
    return str(config.get("executables", {}).get(key, fallback or key))


def current_pbs_job_id() -> str:
    job_id = os.environ.get("SLURM_JOB_ID", os.environ.get("PBS_JOBID", ""))
    return job_id.strip()


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    for _ in range(3):
        try:
            text = path.read_text(encoding="ascii")
            if not text.strip():
                time.sleep(0.05)
                continue
            return json.loads(text)
        except json.JSONDecodeError:
            time.sleep(0.05)
    return default


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_text_atomic(path, json.dumps(payload, indent=2))


def write_text_atomic(path: Path, text: str, encoding: str = "ascii") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    tmp_path.write_text(text, encoding=encoding)
    tmp_path.replace(path)


def query_li_fsi_jobs(root: Path | None = None) -> list[dict[str, str]]:
    # A Slurm case must never inspect/cancel similarly named jobs from the
    # legacy PBS campaigns; their short qsub names are not globally unique.
    if os.environ.get("SLURM_JOB_ID"):
        return []
    base = root or workflow_root()
    job_name_map: dict[str, str] = {}
    cases_dir = base / "cases"
    if cases_dir.exists():
        for case_dir in cases_dir.iterdir():
            if not case_dir.is_dir():
                continue
            metadata = read_json(case_dir / "metadata.json", default={}) or {}
            job_name = str(metadata.get("job_name", "")).strip()
            if job_name:
                job_name_map[job_name] = case_dir.name
    qstat_cmd = executable_path("qstat_executable", base, "qstat")
    proc = subprocess.run([qstat_cmd, "-f", "-1"], text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        return []
    jobs: list[dict[str, str]] = []
    current_job_id = ""
    current_name = ""
    current_state = ""
    current_owner = ""
    current_init_work_dir = ""
    current_user = os.environ.get("USER") or getpass.getuser()

    def flush_current() -> None:
        nonlocal current_job_id, current_name, current_state, current_owner, current_init_work_dir
        if not current_job_id or not current_name:
            return
        owner_user = current_owner.split("@", 1)[0] if current_owner else ""
        if owner_user and owner_user != current_user:
            return
        case_name = None
        if current_init_work_dir:
            workdir = current_init_work_dir.strip()
            marker = str(cases_dir) + "/"
            if workdir.startswith(marker):
                case_name = workdir[len(marker):].split("/", 1)[0]
        if case_name is None:
            case_name = job_name_map.get(current_name)
        if case_name is None:
            return
        jobs.append(
            {
                "job_id": current_job_id,
                "job_name": current_name,
                "case_name": case_name,
                "state": current_state,
            }
        )

    for line in proc.stdout.splitlines():
        if line.startswith("Job Id: "):
            flush_current()
            current_job_id = line.split(":", 1)[1].strip()
            current_name = ""
            current_state = ""
            current_owner = ""
            current_init_work_dir = ""
            continue
        stripped = line.strip()
        if stripped.startswith("Job_Name = "):
            current_name = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("job_state = "):
            current_state = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("Job_Owner = "):
            current_owner = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("init_work_dir = "):
            current_init_work_dir = stripped.split("=", 1)[1].strip()
    flush_current()
    return jobs


def active_li_fsi_jobs_by_case(root: Path | None = None) -> dict[str, list[dict[str, str]]]:
    mapping: dict[str, list[dict[str, str]]] = {}
    for job in query_li_fsi_jobs(root):
        if job["state"] not in {"R", "Q", "H", "E"}:
            continue
        mapping.setdefault(job["case_name"], []).append(job)
    return mapping


def cancel_jobs(job_ids: list[str], root: Path | None = None) -> list[str]:
    if not job_ids:
        return []
    base = root or workflow_root()
    qdel_cmd = executable_path("qdel_executable", base, "qdel")
    cancelled: list[str] = []
    for job_id in job_ids:
        proc = subprocess.run([qdel_cmd, job_id], text=True, capture_output=True, check=False)
        if proc.returncode == 0:
            cancelled.append(job_id)
    return cancelled


def cancel_case_jobs(case_name: str, root: Path | None = None, keep_job_ids: set[str] | None = None) -> list[str]:
    if os.environ.get("SLURM_JOB_ID"):
        return []
    keep = {job_id for job_id in (keep_job_ids or set()) if job_id}
    jobs = active_li_fsi_jobs_by_case(root).get(case_name, [])
    doomed = [job["job_id"] for job in jobs if job["job_id"] not in keep]
    return cancel_jobs(doomed, root=root)


def status_path(case_dir: Path) -> Path:
    return case_dir / "analysis" / "status.json"


def completed_chunk_indices(case_dir: Path) -> list[int]:
    state_dir = case_dir / "state"
    indices = []
    for path in state_dir.glob("restart_*.bin"):
        stem = path.stem
        try:
            indices.append(int(stem.split("_")[-1]))
        except ValueError:
            continue
    return sorted(set(indices))


def completed_chunk_count(case_dir: Path) -> int:
    indices = completed_chunk_indices(case_dir)
    return max(indices, default=0)


def finished_chunk_trajectories(case_dir: Path) -> list[Path]:
    traj_dir = case_dir / "traj"
    paths = []
    for chunk_idx in completed_chunk_indices(case_dir):
        traj = traj_dir / f"chunk_{chunk_idx:04d}.dcd"
        if traj.exists():
            paths.append(traj)
    return paths


def is_final_status(payload: dict[str, Any] | None) -> bool:
    if not payload:
        return False
    status = payload.get("status")
    if status != "failed":
        return status in FINAL_STATUSES
    return payload.get("failure_reason") not in RECOVERABLE_FAILURE_REASONS


def effective_status(payload: dict[str, Any] | None) -> str:
    if not payload:
        return "not_started"
    status = payload.get("status", "not_started")
    if status == "failed" and payload.get("failure_reason") in RECOVERABLE_FAILURE_REASONS:
        return "continue"
    return status
