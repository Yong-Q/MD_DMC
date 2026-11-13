#!/home/qiuyong/.conda/envs/Agent/bin/python3
"""Calculate collective LiFSI+DMC conductivity with the Einstein-Helfand relation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.transformations import NoJump

from dmc_concentration_scan.einstein_helfand import charge_moment_msd, conductivity_from_slope, linear_fit


ROOT = Path(__file__).resolve().parents[1]
FIT_WINDOWS_NS = ((1, 4), (2, 5), (3, 6), (4, 7), (5, 8))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--output-dir", default="", help="write analysis into this isolated directory")
    args = parser.parse_args()
    case = Path(args.case).resolve()
    metadata = json.loads((case / "metadata.json").read_text())
    status_path = case / "analysis" / "status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    n_chunks = int(status.get("completed_chunks", 0))
    expected_chunks = int(metadata["expected_chunks"])
    if n_chunks != expected_chunks:
        raise RuntimeError(f"Expected {expected_chunks} completed chunks, found {n_chunks} in {case}")
    trajectories = [case / "traj" / f"chunk_{idx:04d}.dcd" for idx in range(1, n_chunks + 1)]
    missing = [path for path in trajectories if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing or empty completed DCD chunks in {case}: {missing}")

    print(f"START {case.name}: chunks={len(trajectories)} bytes={sum(p.stat().st_size for p in trajectories)}", flush=True)

    universe = mda.Universe(str(case / "system.data"), [str(path) for path in trajectories], format="LAMMPS")
    print(f"TOPOLOGY_LOADED {case.name}: atoms={len(universe.atoms)} frames={len(universe.trajectory)}", flush=True)
    universe.trajectory.add_transformations(NoJump())
    framework_count = int(metadata["framework_atom_count"])
    mobile = universe.atoms[framework_count:]
    charges = np.asarray(mobile.charges, dtype=float)
    if not np.isfinite(charges).all():
        raise RuntimeError("non-finite mobile atom charges")
    net_charge = float(charges.sum())
    if abs(net_charge) > 1e-3:
        raise RuntimeError(f"mobile system is not neutral: net charge={net_charge:.6g} e")

    moment = np.empty((len(universe.trajectory), 3), dtype=float)
    for frame, timestep in enumerate(universe.trajectory):
        moment[frame] = np.einsum("i,ij->j", charges, mobile.positions)
        if (frame + 1) % 1000 == 0:
            print(f"CHARGE_MOMENT {case.name}: {frame + 1}/{len(moment)} frames", flush=True)
    msd = charge_moment_msd(moment)
    frame_dt_ps = float(metadata["frame_dt_ps"])
    lag_ps = np.arange(len(msd), dtype=float) * frame_dt_ps
    required_end_ps = max(end_ns for _, end_ns in FIT_WINDOWS_NS) * 1000
    if lag_ps[-1] < required_end_ps:
        raise RuntimeError(f"Trajectory reaches only {lag_ps[-1]:g} ps; fit windows require {required_end_ps:g} ps")
    volume = float(metadata["box_volume_A3"])
    temperature = float(metadata["temperature_K"])

    analysis = Path(args.output_dir).resolve() if args.output_dir else case / "analysis"
    analysis.mkdir(exist_ok=True)
    curve_path = analysis / "charge_moment_msd.csv"
    pd.DataFrame({"lag_ps": lag_ps, "charge_moment_msd_e2_A2": msd}).to_csv(curve_path, index=False)
    records = []
    for start_ns, end_ns in FIT_WINDOWS_NS:
        mask = (lag_ps >= start_ns * 1000) & (lag_ps <= end_ns * 1000)
        if mask.sum() < 20:
            continue
        slope, intercept, r2 = linear_fit(lag_ps[mask], msd[mask])
        sigma_S_m = conductivity_from_slope(slope, volume, temperature) if slope > 0 else float("nan")
        records.append({
            "case_name": case.name,
            "base_name": metadata["base_name"],
            "dmc_per_lifsi": int(metadata["dmc_per_lifsi"]),
            "completed_chunks": n_chunks,
            "fit_start_ns": start_ns,
            "fit_end_ns": end_ns,
            "slope_e2_A2_ps": slope,
            "intercept_e2_A2": intercept,
            "r2": r2,
            "conductivity_S_m": sigma_S_m,
            "conductivity_mS_cm": sigma_S_m * 10.0,
            "conductivity_signed_mS_cm": slope * conductivity_from_slope(1.0, volume, temperature) * 10.0,
            "net_mobile_charge_e": net_charge,
            "volume_A3": volume,
            "temperature_K": temperature,
            "method": "collective Einstein-Helfand charge displacement; mobile Li+FSI+DMC partial charges; framework fixed",
        })
    if not records:
        raise RuntimeError("trajectory too short for any preregistered 1-ns conductivity fit window")
    result = analysis / "collective_conductivity_windows.csv"
    pd.DataFrame(records).to_csv(result, index=False)
    print(f"DONE {case.name} frames={len(moment)} mobile_atoms={len(mobile)} net_charge={net_charge:.3e} e windows={len(records)}", flush=True)


if __name__ == "__main__":
    main()
