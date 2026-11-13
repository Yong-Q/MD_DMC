#!/home/qiuyong/.conda/envs/Agent/bin/python3
"""Analyze Li diffusion for a chunked LiFSI-in-MOF trajectory."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import MDAnalysis as mda
import matplotlib
import numpy as np
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.transformations import NoJump
from workflow_common import completed_chunk_count, finished_chunk_trajectories, load_config, workflow_root, write_json

sc = None
DiffusionAnalyzer = None
HAVE_KINISI = False
KINISI_IMPORT_ERROR = "disabled for low-memory recovery after out-of-memory"

matplotlib.use("Agg")
from matplotlib import pyplot as plt

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "teal": "#42949E",
    "red_strong": "#B64342",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
}


def apply_publication_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans", "sans-serif"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.size"] = 7.5
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["xtick.major.width"] = 0.8
    plt.rcParams["ytick.major.width"] = 0.8
    plt.rcParams["xtick.minor.width"] = 0.6
    plt.rcParams["ytick.minor.width"] = 0.6
    plt.rcParams["xtick.major.size"] = 3.0
    plt.rcParams["ytick.major.size"] = 3.0
    plt.rcParams["xtick.minor.size"] = 2.0
    plt.rcParams["ytick.minor.size"] = 2.0
    plt.rcParams["legend.frameon"] = False


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.14,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
        color=PALETTE["neutral_black"],
        ha="left",
        va="bottom",
    )


def save_publication_bundle(fig, base: Path) -> None:
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1 or len(values) < window:
        return values
    kernel = np.ones(window) / window
    padded = np.pad(values, (window // 2, window - 1 - window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def parse_thermo_csv(path: Path) -> list[dict[str, float]]:
    rows = []
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                rows.append(
                    {
                        "step": float(row["step"]),
                        "temp_K": float(row["temp_K"]),
                        "pe_kcal_mol": float(row["pe_kcal_mol"]),
                        "ke_kcal_mol": float(row["ke_kcal_mol"]),
                        "etotal_kcal_mol": float(row["etotal_kcal_mol"]),
                    }
                )
            except (KeyError, ValueError):
                continue
    return rows


def load_thermo_rows(case_dir: Path, metadata: dict) -> list[dict[str, float]]:
    rows = []
    timestep_ps = float(metadata["timestep_fs"]) / 1000.0
    for path in sorted((case_dir / "logs").glob("thermo_*.csv")):
        for row in parse_thermo_csv(path):
            row["time_ps"] = row["step"] * timestep_ps
            rows.append(row)
    rows.sort(key=lambda item: (item["time_ps"], item["step"]))
    deduped = []
    seen = set()
    for row in rows:
        key = (row["step"], row["temp_K"], row["pe_kcal_mol"], row["ke_kcal_mol"], row["etotal_kcal_mol"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def window_trend_sigma(times: np.ndarray, values: np.ndarray) -> float | None:
    if len(times) < 6:
        return None
    centered_t = times - times.mean()
    denom = float(np.dot(centered_t, centered_t))
    if denom <= 0.0:
        return None
    centered_v = values - values.mean()
    slope = float(np.dot(centered_t, centered_v) / denom)
    drift = abs(slope) * float(times[-1] - times[0])
    sigma = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    sigma = max(sigma, 1.0e-12)
    return drift / sigma


def assess_thermo_equilibration(rows: list[dict[str, float]], config: dict, target_temperature_K: float) -> dict:
    payload = {"thermo_equilibrated": False}
    if not rows:
        return payload
    times = np.array([row["time_ps"] for row in rows], dtype=float)
    temps = np.array([row["temp_K"] for row in rows], dtype=float)
    pes = np.array([row["pe_kcal_mol"] for row in rows], dtype=float)
    window_ps = float(config["analysis"]["thermo_window_ps"])
    mask = times >= max(times[-1] - window_ps, times[0])
    if mask.sum() < 6:
        return payload
    tw = times[mask]
    tempw = temps[mask]
    pew = pes[mask]
    temp_mean = float(np.mean(tempw))
    temp_std = float(np.std(tempw, ddof=1)) if len(tempw) > 1 else 0.0
    pe_mean = float(np.mean(pew))
    pe_std = float(np.std(pew, ddof=1)) if len(pew) > 1 else 0.0
    temp_drift_sigma = window_trend_sigma(tw, tempw)
    pe_drift_sigma = window_trend_sigma(tw, pew)
    temp_ok = abs(temp_mean - target_temperature_K) <= float(config["analysis"]["temperature_tolerance_K"])
    temp_trend_ok = temp_drift_sigma is not None and temp_drift_sigma <= float(config["analysis"]["max_temp_drift_sigma"])
    pe_trend_ok = pe_drift_sigma is not None and pe_drift_sigma <= float(config["analysis"]["max_pe_drift_sigma"])
    payload.update(
        {
            "thermo_window_ps": window_ps,
            "thermo_window_points": int(mask.sum()),
            "thermo_temp_mean_K": temp_mean,
            "thermo_temp_std_K": temp_std,
            "thermo_pe_mean_kcal_mol": pe_mean,
            "thermo_pe_std_kcal_mol": pe_std,
            "thermo_temp_drift_sigma": temp_drift_sigma,
            "thermo_pe_drift_sigma": pe_drift_sigma,
            "thermo_equilibrated": bool(temp_ok and temp_trend_ok and pe_trend_ok),
        }
    )
    return payload


def apparent_diffusion_cm2_s(lag_ps: np.ndarray, msd_A2: np.ndarray) -> np.ndarray:
    values = np.full_like(msd_A2, np.nan, dtype=float)
    mask = lag_ps > 0.0
    values[mask] = msd_A2[mask] * 1.0e-4 / (6.0 * lag_ps[mask])
    return values


def assess_diffusion_plateau(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    start_ps: float,
    config: dict,
    end_ps: float | None = None,
) -> dict:
    payload = {"d_plateau_ok": False}
    d_app = apparent_diffusion_cm2_s(lag_ps, msd_A2)
    mask = (lag_ps >= start_ps) & np.isfinite(d_app) & (d_app > 0.0)
    if end_ps is not None:
        mask &= lag_ps <= end_ps
    if mask.sum() < 12:
        return payload
    tail = d_app[mask]
    tail = tail[max(0, len(tail) // 2) :]
    if len(tail) < 9:
        return payload
    median = float(np.median(tail))
    if median <= 0.0:
        return payload
    q10, q90 = np.quantile(tail, [0.1, 0.9])
    rel_spread = float((q90 - q10) / median)
    block_count = int(config["analysis"]["d_block_count"])
    block_values = [block for block in np.array_split(tail, block_count) if len(block) > 0]
    block_medians = np.array([np.median(block) for block in block_values], dtype=float)
    rel_block_dev = float(np.max(np.abs(block_medians - median)) / median) if len(block_medians) else np.inf
    payload.update(
        {
            "d_app_cm2_s_median": median,
            "d_plateau_rel_spread": rel_spread,
            "d_block_medians_cm2_s": [float(value) for value in block_medians],
            "d_block_rel_deviation": rel_block_dev,
            "d_plateau_ok": bool(
                rel_spread <= float(config["analysis"]["d_plateau_rel_spread_max"])
                and rel_block_dev <= float(config["analysis"]["d_block_rel_deviation_max"])
            ),
        }
    )
    return payload


def find_diffusive_start(lag_ps: np.ndarray, msd_A2: np.ndarray, config: dict) -> tuple[float | None, float | None, int]:
    """Legacy helper retained for diagnostics.

    The actual fitting path now prefers an explicit diffusive window selection
    so that the fit window and the displayed local slope refer to the same
    time span.
    """
    max_lag_ps = reliable_lag_limit_ps(lag_ps, config)
    mask = (lag_ps > 0.0) & (msd_A2 > 0.0)
    if max_lag_ps is not None:
        mask &= lag_ps <= max_lag_ps
    if mask.sum() < 20:
        return None, None, 0
    log_t = np.log(lag_ps[mask])
    log_msd = np.log(msd_A2[mask])
    slope = moving_average(np.gradient(log_msd, log_t), 7)
    valid = (
        (lag_ps[mask] >= float(config["analysis"]["min_start_ps"]))
        & (slope >= float(config["analysis"]["min_loglog_slope"]))
        & (slope <= float(config["analysis"]["max_loglog_slope"]))
    )
    best = []
    current = []
    for idx, ok in enumerate(valid):
        if ok:
            current.append(idx)
        else:
            if len(current) > len(best):
                best = current
            current = []
    if len(current) > len(best):
        best = current
    if len(best) < int(config["analysis"]["min_slope_window_frames"]):
        return None, float(slope[-1]), len(best)
    start = float(lag_ps[mask][best[0]])
    tail = float(np.mean(slope[max(0, best[-1] - 4) : best[-1] + 1]))
    return start, tail, len(best)


def linear_diffusion_fit_window(x: np.ndarray, y: np.ndarray) -> tuple[float | None, float | None]:
    if len(x) < 8:
        return None, None
    try:
        coeffs, cov = np.polyfit(x, y, deg=1, cov=True)
    except (np.linalg.LinAlgError, ValueError):
        return None, None
    slope_A2_per_ps = float(coeffs[0])
    slope_std = float(np.sqrt(max(cov[0, 0], 0.0))) if cov.size else 0.0
    if slope_A2_per_ps <= 0.0:
        return None, None
    mean = slope_A2_per_ps * 1.0e-4 / 6.0
    std = slope_std * 1.0e-4 / 6.0
    return mean, std


def linear_fit_stats_window(x: np.ndarray, y: np.ndarray) -> tuple[float | None, float | None, float | None]:
    if len(x) < 8:
        return None, None, None
    try:
        coeffs, cov = np.polyfit(x, y, deg=1, cov=True)
    except (np.linalg.LinAlgError, ValueError):
        return None, None, None
    slope_A2_per_ps = float(coeffs[0])
    slope_std = float(np.sqrt(max(cov[0, 0], 0.0))) if cov.size else 0.0
    if slope_A2_per_ps <= 0.0:
        return None, None, None
    y_hat = np.polyval(coeffs, x)
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 0.0 else max(0.0, 1.0 - ss_res / ss_tot)
    mean = slope_A2_per_ps * 1.0e-4 / 6.0
    std = slope_std * 1.0e-4 / 6.0
    return mean, std, r2


def reliable_lag_limit_ps(lag_ps: np.ndarray, config: dict) -> float | None:
    positive = lag_ps[lag_ps > 0.0]
    n_points = len(positive)
    if n_points < 8:
        return None
    analysis_cfg = config["analysis"]
    min_origins = int(analysis_cfg.get("fit_min_origins", 500))
    min_fraction = float(analysis_cfg.get("fit_min_origins_fraction", 0.10))
    required_origins = max(min_origins, int(np.ceil(n_points * min_fraction)))
    max_index = n_points - required_origins
    if max_index < 1:
        return float(positive[-1])
    max_index = min(max_index, n_points - 1)
    return float(positive[max_index])


def linear_diffusion_fit(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    start_ps: float,
    end_ps: float | None = None,
) -> tuple[float | None, float | None]:
    mask = mask_window(lag_ps, msd_A2, start_ps, end_ps)
    if mask.sum() < 8:
        return None, None
    return linear_diffusion_fit_window(lag_ps[mask], msd_A2[mask])


def run_kinisi_diffusion(universe: mda.Universe, li_type: int, frame_dt_ps: float, start_ps: float) -> tuple[float | None, float | None, str | None]:
    if not HAVE_KINISI:
        return None, None, None
    analyzer = DiffusionAnalyzer.from_universe(
        trajectory=universe,
        specie=str(li_type),
        time_step=frame_dt_ps * sc.Unit("ps"),
        step_skip=1 * sc.Unit("dimensionless"),
        progress=False,
    )
    analyzer.diffusion(start_ps * sc.Unit("ps"), progress=False)
    return (
        float(sc.mean(analyzer.D).value),
        float(sc.std(analyzer.D, ddof=1).value),
        str(analyzer.D.unit),
    )


def mask_window(lag_ps: np.ndarray, msd_A2: np.ndarray, start_ps: float, end_ps: float | None = None) -> np.ndarray:
    mask = np.isfinite(lag_ps) & np.isfinite(msd_A2) & (lag_ps >= start_ps) & (lag_ps > 0.0)
    if end_ps is not None:
        mask &= lag_ps <= end_ps
    return mask


def robust_window_diffusion_fit(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    start_ps: float,
    end_ps: float | None = None,
    *,
    min_points: int,
    allow_zero_floor: bool = False,
    zero_method: str = "window_zero_floor",
) -> tuple[float | None, float | None, float | None, float | None, str | None]:
    mask = mask_window(lag_ps, msd_A2, start_ps, end_ps)
    if mask.sum() < min_points:
        return None, None, None, None, None

    x_window = lag_ps[mask]
    y_window = msd_A2[mask]
    fit_start = float(x_window[0])
    fit_end = float(x_window[-1])

    mean, std = linear_diffusion_fit_window(x_window, y_window)
    if mean is not None:
        return mean, std, fit_start, fit_end, "tail_linear"

    y_envelope = np.maximum.accumulate(y_window)
    mean, std = linear_diffusion_fit_window(x_window, y_envelope)
    if mean is not None:
        return mean, std, fit_start, fit_end, "tail_monotone_envelope"

    dx = np.diff(x_window)
    dy = np.diff(y_envelope)
    slope_samples = dy[dx > 0.0] / dx[dx > 0.0]
    positive_slopes = slope_samples[slope_samples > 0.0]
    if len(positive_slopes) >= 3:
        slope = float(np.median(positive_slopes))
        slope_std = float(np.std(positive_slopes, ddof=1)) if len(positive_slopes) > 1 else 0.0
        mean = slope * 1.0e-4 / 6.0
        std = slope_std * 1.0e-4 / 6.0
        return mean, std, fit_start, fit_end, "tail_increment_median"

    if allow_zero_floor:
        return 0.0, 0.0, fit_start, fit_end, zero_method

    return None, None, fit_start, fit_end, None


def tail_window(lag_ps: np.ndarray, msd_A2: np.ndarray, config: dict) -> tuple[np.ndarray | None, np.ndarray | None, float | None]:
    mask = np.isfinite(lag_ps) & np.isfinite(msd_A2) & (lag_ps > 0.0)
    max_lag_ps = reliable_lag_limit_ps(lag_ps, config)
    if max_lag_ps is not None:
        mask &= lag_ps <= max_lag_ps
    if mask.sum() < int(config["analysis"]["tail_fit_min_points"]):
        return None, None, None
    x = lag_ps[mask]
    y = msd_A2[mask]
    start_idx = int(max(0, np.floor((1.0 - float(config["analysis"]["tail_fit_fraction"])) * len(x))))
    x_tail = x[start_idx:]
    y_tail = y[start_idx:]
    if len(x_tail) < int(config["analysis"]["tail_fit_min_points"]):
        return None, None, None
    return x_tail, y_tail, float(x_tail[0])


def tail_linear_diffusion_fit(lag_ps: np.ndarray, msd_A2: np.ndarray, config: dict) -> tuple[float | None, float | None, float | None]:
    x_tail, y_tail, start_ps = tail_window(lag_ps, msd_A2, config)
    if x_tail is None or y_tail is None or start_ps is None:
        return None, None, None
    mean, std = linear_diffusion_fit_window(x_tail, y_tail)
    return mean, std, start_ps


def robust_tail_diffusion_fit(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    config: dict,
    *,
    allow_zero_floor: bool = False,
    zero_method: str = "tail_zero_floor",
) -> tuple[float | None, float | None, float | None, str | None]:
    x_tail, y_tail, start_ps = tail_window(lag_ps, msd_A2, config)
    if x_tail is None or y_tail is None or start_ps is None:
        return None, None, None, None
    mean, std, fit_start, _fit_end, method = robust_window_diffusion_fit(
        x_tail,
        y_tail,
        float(x_tail[0]),
        min_points=int(config["analysis"]["tail_fit_min_points"]),
        allow_zero_floor=allow_zero_floor,
        zero_method=zero_method,
    )
    return mean, std, fit_start, method


def window_slope_metrics(
    lag_ps: np.ndarray, loglog_slope: np.ndarray, start_ps: float, end_ps: float | None = None
) -> tuple[float | None, float | None]:
    mask = (lag_ps >= start_ps) & np.isfinite(loglog_slope)
    if end_ps is not None:
        mask &= lag_ps <= end_ps
    if mask.sum() < 8:
        return None, None
    slopes = loglog_slope[mask]
    median = float(np.median(slopes))
    q25, q75 = np.quantile(slopes, [0.25, 0.75])
    return median, float(q75 - q25)


def candidate_start_times(lag_ps: np.ndarray, config: dict) -> list[float]:
    max_lag_ps = reliable_lag_limit_ps(lag_ps, config)
    positive = lag_ps[lag_ps > 0.0]
    if max_lag_ps is not None:
        positive = positive[positive <= max_lag_ps]
    if len(positive) < 8:
        return []
    total_ps = float(positive[-1])
    analysis_cfg = config["analysis"]
    min_window_ps = float(analysis_cfg.get("dynamic_window_min_ps", 1000.0))
    step_ps = float(analysis_cfg.get("dynamic_window_candidate_step_ps", 250.0))
    start_fraction = float(analysis_cfg.get("dynamic_window_search_start_fraction", 0.2))
    end_fraction = float(analysis_cfg.get("dynamic_window_search_end_fraction", 0.8))
    search_start = max(float(analysis_cfg["min_start_ps"]), total_ps * start_fraction)
    latest_start = max(search_start, total_ps - min_window_ps)
    search_end = min(total_ps * end_fraction, latest_start)
    starts: list[float] = []
    current = search_start
    while current <= search_end + 1.0e-9:
        if total_ps - current >= min_window_ps:
            starts.append(float(round(current, 6)))
        current += step_ps
    if total_ps - latest_start >= min_window_ps:
        starts.append(float(round(latest_start, 6)))
    return sorted(set(starts))


def candidate_end_times(lag_ps: np.ndarray, start_ps: float, config: dict) -> list[float]:
    max_lag_ps = reliable_lag_limit_ps(lag_ps, config)
    positive = lag_ps[lag_ps > 0.0]
    if max_lag_ps is not None:
        positive = positive[positive <= max_lag_ps]
    if len(positive) < 8:
        return []
    total_ps = float(positive[-1])
    analysis_cfg = config["analysis"]
    min_duration_ps = max(
        float(analysis_cfg.get("dynamic_subwindow_min_ps", 400.0)),
        total_ps * float(analysis_cfg.get("dynamic_subwindow_min_fraction", 0.0)),
    )
    step_ps = float(analysis_cfg.get("dynamic_window_candidate_step_ps", 250.0))
    first_end = start_ps + min_duration_ps
    if first_end > total_ps:
        return []
    ends: list[float] = []
    current = first_end
    while current <= total_ps + 1.0e-9:
        ends.append(float(round(current, 6)))
        current += step_ps
    ends.append(float(round(total_ps, 6)))
    return sorted(set(ends))


def select_dynamic_fit_window(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    loglog_slope: np.ndarray,
    config: dict,
) -> dict | None:
    analysis_cfg = config["analysis"]
    min_slope = float(analysis_cfg["min_loglog_slope"])
    max_slope = float(analysis_cfg["max_loglog_slope"])
    max_slope_iqr = float(analysis_cfg.get("dynamic_max_slope_iqr", 0.20))
    max_abs_slope_deviation = float(analysis_cfg.get("dynamic_max_abs_slope_deviation", 0.05))
    require_plateau = bool(analysis_cfg.get("dynamic_require_plateau", True))
    min_window_ps = float(
        analysis_cfg.get(
            "dynamic_accept_min_window_ps",
            max(
                analysis_cfg.get("dynamic_subwindow_min_ps", 400.0),
                2.0 * analysis_cfg.get("dynamic_window_candidate_step_ps", 250.0),
            ),
        )
    )
    candidates = []
    for start_ps in candidate_start_times(lag_ps, config):
        for end_ps in candidate_end_times(lag_ps, start_ps, config):
            mean, std, fit_start, fit_end, method = robust_window_diffusion_fit(
                lag_ps,
                msd_A2,
                start_ps,
                end_ps=end_ps,
                min_points=int(config["analysis"]["tail_fit_min_points"]),
                allow_zero_floor=False,
            )
            if mean is None or fit_start is None or fit_end is None or method is None or mean <= 0.0:
                continue
            plateau = assess_diffusion_plateau(lag_ps, msd_A2, fit_start, config, end_ps=fit_end)
            slope_median, slope_iqr = window_slope_metrics(lag_ps, loglog_slope, fit_start, end_ps=fit_end)
            rel_std = float(std / mean) if std is not None and mean > 0.0 else 1.0
            rel_spread = float(plateau.get("d_plateau_rel_spread", 2.0))
            rel_block_dev = float(plateau.get("d_block_rel_deviation", 2.0))
            slope_penalty = abs((slope_median if slope_median is not None else 0.0) - 1.0) if slope_median is not None else 1.5
            slope_iqr_penalty = slope_iqr if slope_iqr is not None else 1.5
            plateau_penalty = 0.0 if plateau.get("d_plateau_ok", False) else 1.5
            window_span_ps = float(fit_end - fit_start)
            slope_ok = (
                slope_median is not None
                and min_slope <= float(slope_median) <= max_slope
                and abs(float(slope_median) - 1.0) <= max_abs_slope_deviation
                and slope_iqr is not None
                and float(slope_iqr) <= max_slope_iqr
            )
            plateau_ok = bool(plateau.get("d_plateau_ok", False))
            if not slope_ok:
                continue
            if require_plateau and not plateau_ok:
                continue
            if window_span_ps < min_window_ps:
                continue
            method_penalty = {
                "tail_linear": 0.0,
                "tail_monotone_envelope": 0.15,
                "tail_increment_median": 0.3,
            }.get(method, 0.5)
            score = (
                2.5 * plateau_penalty
                + 2.0 * slope_penalty
                + 0.5 * min(slope_iqr_penalty, 2.0)
                + 0.5 * min(rel_spread, 2.0)
                + 0.5 * min(rel_block_dev, 2.0)
                + 0.35 * min(rel_std, 1.0)
                + method_penalty
            )
            candidates.append(
                {
                    "start_ps": fit_start,
                    "end_ps": fit_end,
                    "mean": mean,
                    "std": 0.0 if std is None else std,
                    "method": method,
                    "score": score,
                    "plateau_ok": bool(plateau.get("d_plateau_ok", False)),
                    "rel_spread": rel_spread,
                    "rel_block_dev": rel_block_dev,
                    "slope_median": slope_median,
                    "slope_iqr": slope_iqr,
                }
            )
    if not candidates:
        return None
    candidates.sort(
        key=lambda item: (
            item["score"],
            0 if item["plateau_ok"] else 1,
            abs((item["slope_median"] if item["slope_median"] is not None else 0.0) - 1.0),
            -(item["end_ps"] - item["start_ps"]),
            -item["start_ps"],
        )
    )
    return candidates[0]


def select_diffusive_window(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    loglog_slope: np.ndarray,
    config: dict,
) -> dict | None:
    analysis_cfg = config["analysis"]
    min_slope = float(analysis_cfg["min_loglog_slope"])
    max_slope = float(analysis_cfg["max_loglog_slope"])
    max_slope_iqr = float(analysis_cfg.get("diffusive_max_slope_iqr", 0.22))
    max_abs_slope_deviation = float(analysis_cfg.get("diffusive_max_abs_slope_deviation", 0.10))
    min_window_ps = float(
        analysis_cfg.get(
            "diffusive_min_window_ps",
            max(
                analysis_cfg.get("dynamic_subwindow_min_ps", 400.0),
                2.0 * analysis_cfg.get("dynamic_window_candidate_step_ps", 250.0),
            ),
        )
    )
    require_plateau = bool(analysis_cfg.get("diffusive_require_plateau", True))
    candidates = []
    for start_ps in candidate_start_times(lag_ps, config):
        for end_ps in candidate_end_times(lag_ps, start_ps, config):
            mask = mask_window(lag_ps, msd_A2, start_ps, end_ps)
            if mask.sum() < int(config["analysis"]["tail_fit_min_points"]):
                continue
            x_window = lag_ps[mask]
            y_window = msd_A2[mask]
            fit_start = float(x_window[0])
            fit_end = float(x_window[-1])
            window_span_ps = fit_end - fit_start
            if window_span_ps < min_window_ps:
                continue
            mean, std, r2 = linear_fit_stats_window(x_window, y_window)
            if mean is None or std is None or r2 is None:
                continue
            slope_median, slope_iqr = window_slope_metrics(lag_ps, loglog_slope, fit_start, end_ps=fit_end)
            if slope_median is None or slope_iqr is None:
                continue
            slope_ok = (
                min_slope <= float(slope_median) <= max_slope
                and abs(float(slope_median) - 1.0) <= max_abs_slope_deviation
                and float(slope_iqr) <= max_slope_iqr
            )
            if not slope_ok:
                continue
            plateau = assess_diffusion_plateau(lag_ps, msd_A2, fit_start, config, end_ps=fit_end)
            plateau_ok = bool(plateau.get("d_plateau_ok", False))
            if require_plateau and not plateau_ok:
                continue
            rel_std = float(std / mean) if mean > 0.0 else 1.0
            rel_spread = float(plateau.get("d_plateau_rel_spread", 2.0))
            rel_block_dev = float(plateau.get("d_block_rel_deviation", 2.0))
            score = (
                4.0 * abs(float(slope_median) - 1.0)
                + 1.0 * min(float(slope_iqr), 1.0)
                + 0.5 * min(rel_std, 1.0)
                + 0.35 * min(rel_spread, 2.0)
                + 0.35 * min(rel_block_dev, 2.0)
                + 0.2 * (1.0 - r2)
            )
            candidates.append(
                {
                    "start_ps": fit_start,
                    "end_ps": fit_end,
                    "mean": mean,
                    "std": std,
                    "r2": r2,
                    "score": score,
                    "plateau_ok": plateau_ok,
                    "rel_spread": rel_spread,
                    "rel_block_dev": rel_block_dev,
                    "slope_median": slope_median,
                    "slope_iqr": slope_iqr,
                }
            )
    if not candidates:
        return None
    candidates.sort(
        key=lambda item: (
            item["score"],
            item["start_ps"],
            -(item["end_ps"] - item["start_ps"]),
            -item["r2"],
        )
    )
    return candidates[0]


def select_dynamic_linear_window(
    lag_ps: np.ndarray,
    msd_A2: np.ndarray,
    loglog_slope: np.ndarray,
    config: dict,
) -> dict | None:
    analysis_cfg = config["analysis"]
    min_r2 = float(analysis_cfg.get("dynamic_linear_min_r2", 0.995))
    min_window_ps = float(
        analysis_cfg.get(
            "dynamic_accept_min_window_ps",
            max(
                analysis_cfg.get("dynamic_subwindow_min_ps", 400.0),
                2.0 * analysis_cfg.get("dynamic_window_candidate_step_ps", 250.0),
            ),
        )
    )
    require_plateau = bool(analysis_cfg.get("dynamic_require_plateau", True))
    candidates = []
    for start_ps in candidate_start_times(lag_ps, config):
        for end_ps in candidate_end_times(lag_ps, start_ps, config):
            mask = mask_window(lag_ps, msd_A2, start_ps, end_ps)
            if mask.sum() < int(config["analysis"]["tail_fit_min_points"]):
                continue
            x_window = lag_ps[mask]
            y_window = msd_A2[mask]
            fit_start = float(x_window[0])
            fit_end = float(x_window[-1])
            window_span_ps = fit_end - fit_start
            if window_span_ps < min_window_ps:
                continue
            mean, std, r2 = linear_fit_stats_window(x_window, y_window)
            if mean is None or std is None or r2 is None or r2 < min_r2:
                continue
            plateau = assess_diffusion_plateau(lag_ps, msd_A2, fit_start, config, end_ps=fit_end)
            plateau_ok = bool(plateau.get("d_plateau_ok", False))
            if require_plateau and not plateau_ok:
                continue
            slope_median, slope_iqr = window_slope_metrics(lag_ps, loglog_slope, fit_start, end_ps=fit_end)
            rel_std = float(std / mean) if mean > 0.0 else 1.0
            rel_spread = float(plateau.get("d_plateau_rel_spread", 2.0))
            rel_block_dev = float(plateau.get("d_block_rel_deviation", 2.0))
            score = (
                4.0 * (1.0 - r2)
                + 0.5 * min(rel_std, 1.0)
                + 0.5 * min(rel_spread, 2.0)
                + 0.5 * min(rel_block_dev, 2.0)
            )
            candidates.append(
                {
                    "start_ps": fit_start,
                    "end_ps": fit_end,
                    "mean": mean,
                    "std": std,
                    "method": "linear_window",
                    "score": score,
                    "r2": r2,
                    "plateau_ok": plateau_ok,
                    "rel_spread": rel_spread,
                    "rel_block_dev": rel_block_dev,
                    "slope_median": slope_median,
                    "slope_iqr": slope_iqr,
                }
            )
    if not candidates:
        return None
    candidates.sort(
        key=lambda item: (
            item["score"],
            -item["r2"],
            0 if item["plateau_ok"] else 1,
            -(item["end_ps"] - item["start_ps"]),
            item["start_ps"],
        )
    )
    return candidates[0]


def topology_path(case_dir: Path) -> Path | None:
    state_dir = case_dir / "state"
    data_files = sorted(state_dir.glob("data_*.data"))
    if data_files:
        return data_files[-1]
    eq_path = state_dir / "equilibrated_npt.data"
    if eq_path.exists():
        return eq_path
    system_path = case_dir / "system.data"
    if system_path.exists():
        return system_path
    return None


def run_einstein_msd(universe: mda.Universe, selection: str) -> tuple[np.ndarray, bool]:
    msd = EinsteinMSD(universe, select=selection, msd_type="xyz", fft=True)
    try:
        msd.run()
        return np.asarray(msd.results.timeseries, dtype=float), True
    except ImportError as exc:
        if "tidynamics" not in str(exc):
            raise
        msd = EinsteinMSD(universe, select=selection, msd_type="xyz", fft=False)
        msd.run()
        return np.asarray(msd.results.timeseries, dtype=float), False


def readable_trajectories(top_path: Path, trajs: list[Path]) -> tuple[list[Path], dict[str, str]]:
    valid: list[Path] = []
    skipped: dict[str, str] = {}
    for path in trajs:
        try:
            probe = mda.Universe(str(top_path), [str(path)], format="LAMMPS")
            probe.trajectory.add_transformations(NoJump())
            frame_count = len(probe.trajectory)
            if frame_count <= 0:
                raise ValueError("empty trajectory")
            _ = probe.trajectory[0]
            _ = probe.trajectory[-1]
            valid.append(path)
        except Exception as exc:
            skipped[str(path)] = str(exc)
    return valid, skipped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    args = parser.parse_args()

    config = load_config(workflow_root())
    case_dir = Path(args.case).resolve()
    metadata = json.loads((case_dir / "metadata.json").read_text(encoding="ascii"))
    completed_chunks = completed_chunk_count(case_dir)
    trajs = finished_chunk_trajectories(case_dir)
    out_dir = case_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    msd_csv = out_dir / "msd.csv"
    li_msd_csv = out_dir / "li_msd.csv"
    fsi_msd_csv = out_dir / "fsi_msd.csv"
    msd_base = out_dir / "msd"
    msd_png = msd_base.with_suffix(".png")
    msd_svg = msd_base.with_suffix(".svg")
    msd_pdf = msd_base.with_suffix(".pdf")
    msd_tiff = msd_base.with_suffix(".tiff")

    if not trajs:
        payload = {"status": "no_trajectory", "completed_chunks": completed_chunks, "max_chunks": metadata.get("max_chunks", 0)}
        write_json(out_dir / "status.json", payload)
        return

    top_path = topology_path(case_dir)
    if top_path is None:
        payload = {"status": "missing_topology", "completed_chunks": completed_chunks, "max_chunks": metadata.get("max_chunks", 0)}
        write_json(out_dir / "status.json", payload)
        return

    valid_trajs, skipped_trajs = readable_trajectories(top_path, trajs)
    if not valid_trajs:
        payload = {
            "status": "continue",
            "selection": f"type {metadata.get('li_type', 0)}",
            "completed_chunks": completed_chunks,
            "max_chunks": metadata["max_chunks"],
            "topology_path": str(top_path),
            "trajectory_count": 0,
            "trajectory_candidate_count": len(trajs),
            "skipped_trajectories": skipped_trajs,
        }
        write_json(out_dir / "status.json", payload)
        return

    universe = mda.Universe(str(top_path), [str(path) for path in valid_trajs], format="LAMMPS")
    universe.trajectory.add_transformations(NoJump())
    guest_types = metadata.get("guest_types", {})
    li_type = int(guest_types.get("li", metadata.get("li_type")))
    fsi_anchor_type = int(guest_types.get("fsi_anchor", metadata.get("fsi_anchor_type", guest_types.get("fsi_all", [0])[0])))
    li_selection = f"type {li_type}"
    fsi_selection = f"type {fsi_anchor_type}"
    msd_values, fft_mode = run_einstein_msd(universe, li_selection)
    fsi_msd_values, fsi_fft_mode = run_einstein_msd(universe, fsi_selection)

    lag_ps = np.arange(len(msd_values), dtype=float) * float(metadata["frame_dt_ps"])
    legacy_start_ps, tail_slope, slope_window_frames = find_diffusive_start(lag_ps, msd_values, config)
    total_ns = float(lag_ps[-1] / 1000.0) if len(lag_ps) else 0.0
    thermo_rows = load_thermo_rows(case_dir, metadata)
    thermo_payload = assess_thermo_equilibration(thermo_rows, config, float(metadata["temperature_K"]))

    loglog_slope = np.full_like(msd_values, np.nan, dtype=float)
    positive_mask = (lag_ps > 0.0) & (msd_values > 0.0)
    if positive_mask.sum() >= 3:
        loglog_slope[positive_mask] = moving_average(
            np.gradient(np.log(msd_values[positive_mask]), np.log(lag_ps[positive_mask])),
            7,
        )

    with msd_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lag_ps", "li_msd_A2", "fsi_anchor_msd_A2", "li_loglog_slope"])
        for lag, li_val, fsi_val, slope in zip(lag_ps, msd_values, fsi_msd_values, loglog_slope):
            writer.writerow([f"{lag:.8f}", f"{li_val:.8f}", f"{fsi_val:.8f}", "" if np.isnan(slope) else f"{slope:.8f}"])
    with li_msd_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lag_ps", "msd_A2"])
        for lag, val in zip(lag_ps, msd_values):
            writer.writerow([f"{lag:.8f}", f"{val:.8f}"])
    with fsi_msd_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lag_ps", "msd_A2"])
        for lag, val in zip(lag_ps, fsi_msd_values):
            writer.writerow([f"{lag:.8f}", f"{val:.8f}"])

    payload = {
        "status": "continue",
        "selection": li_selection,
        "li_selection": li_selection,
        "fsi_selection": fsi_selection,
        "topology_path": str(top_path),
        "trajectory_count": len(valid_trajs),
        "trajectory_candidate_count": len(trajs),
        "skipped_trajectories": skipped_trajs,
        "completed_chunks": completed_chunks,
        "max_chunks": metadata["max_chunks"],
        "total_ns": total_ns,
        "msd_fft_mode": fft_mode,
        "fsi_msd_fft_mode": fsi_fft_mode,
        "tail_loglog_slope": tail_slope,
        "slope_window_frames": slope_window_frames,
        "diffusive_start_ps": legacy_start_ps,
        "msd_csv": str(msd_csv),
        "li_msd_csv": str(li_msd_csv),
        "fsi_msd_csv": str(fsi_msd_csv),
        "msd_png": str(msd_png),
        "msd_svg": str(msd_svg),
        "msd_pdf": str(msd_pdf),
        "msd_tiff": str(msd_tiff),
    }
    payload.update(thermo_payload)
    diffusive_window = select_diffusive_window(lag_ps, msd_values, loglog_slope, config)
    start_ps = diffusive_window["start_ps"] if diffusive_window is not None else None
    fit_start_ps = start_ps
    fit_end_ps = diffusive_window["end_ps"] if diffusive_window is not None else (float(lag_ps[-1]) if len(lag_ps) else start_ps)
    fit_window_source = "diffusive_window" if diffusive_window is not None else ""
    fit_window_is_fallback = diffusive_window is None
    if diffusive_window is not None:
        payload["diffusive_window_slope_median"] = diffusive_window["slope_median"]
        payload["diffusive_window_slope_iqr"] = diffusive_window["slope_iqr"]
        payload["diffusive_window_r2"] = diffusive_window["r2"]
    final_kinisi_allowed = completed_chunks >= int(metadata["max_chunks"])

    if diffusive_window is not None and start_ps is not None:
        linear_mean, linear_std = linear_diffusion_fit(lag_ps, msd_values, start_ps, end_ps=fit_end_ps)
        if linear_mean is not None:
            payload["linear_diffusion_cm2_s_mean"] = linear_mean
            payload["linear_diffusion_cm2_s_std"] = 0.0 if linear_std is None else linear_std
            payload["linear_fit_start_ps"] = start_ps
            payload["linear_fit_end_ps"] = fit_end_ps
        payload.update(assess_diffusion_plateau(lag_ps, msd_values, start_ps, config))
        diffusive_candidate = (
            total_ns >= float(config["analysis"]["min_diffusive_ns"])
            and payload.get("linear_diffusion_cm2_s_mean", 0.0) > 0.0
            and payload.get("thermo_equilibrated", False)
            and payload.get("d_plateau_ok", False)
        )
        if diffusive_candidate:
            payload["diffusion_criteria_met"] = True
            # This concentration scan requires the same full production length
            # for every ratio, even if self-diffusion converges earlier.
            if completed_chunks >= int(metadata["max_chunks"]):
                payload["status"] = "diffusive"
                final_kinisi_allowed = True
        if final_kinisi_allowed:
            if HAVE_KINISI:
                mean, std, unit = run_kinisi_diffusion(universe, li_type, float(metadata["frame_dt_ps"]), start_ps)
                if mean is not None:
                    payload["kinisi_diffusion_cm2_s_mean"] = mean
                    payload["kinisi_diffusion_cm2_s_std"] = 0.0 if std is None else std
                    payload["kinisi_unit"] = unit
                    payload["diffusion_method"] = "kinisi"
            else:
                payload["kinisi_available"] = False
                payload["kinisi_import_error"] = KINISI_IMPORT_ERROR
                if linear_mean is not None:
                    payload["diffusion_method"] = "linear_msd_recovery_no_kinisi"
        else:
            payload["kinisi_deferred_until_final"] = True
    else:
        zero_method = "tail_zero_floor" if completed_chunks >= int(metadata["max_chunks"]) else "tail_zero_candidate"
        dynamic_window = select_dynamic_fit_window(lag_ps, msd_values, loglog_slope, config)
        if dynamic_window is not None:
            payload["linear_diffusion_cm2_s_mean"] = dynamic_window["mean"]
            payload["linear_diffusion_cm2_s_std"] = dynamic_window["std"]
            payload["linear_fit_start_ps"] = dynamic_window["start_ps"]
            payload["linear_fit_end_ps"] = dynamic_window["end_ps"]
            payload["fallback_diffusion_cm2_s_mean"] = dynamic_window["mean"]
            payload["fallback_diffusion_cm2_s_std"] = dynamic_window["std"]
            payload["fallback_fit_start_ps"] = dynamic_window["start_ps"]
            payload["fallback_fit_end_ps"] = dynamic_window["end_ps"]
            payload["fallback_quality"] = dynamic_window["method"]
            payload["fallback_window_score"] = dynamic_window["score"]
            payload["fallback_window_slope_median"] = dynamic_window["slope_median"]
            payload["fallback_window_slope_iqr"] = dynamic_window["slope_iqr"]
            payload["fallback_window_plateau_ok"] = dynamic_window["plateau_ok"]
            fit_start_ps = dynamic_window["start_ps"]
            fit_end_ps = dynamic_window["end_ps"]
            fit_window_source = f"fallback_{dynamic_window['method']}"
            fit_window_is_fallback = True
            payload["diffusion_method"] = fit_window_source
        else:
            linear_window = select_dynamic_linear_window(lag_ps, msd_values, loglog_slope, config)
            if linear_window is not None:
                payload["linear_diffusion_cm2_s_mean"] = linear_window["mean"]
                payload["linear_diffusion_cm2_s_std"] = linear_window["std"]
                payload["linear_fit_start_ps"] = linear_window["start_ps"]
                payload["linear_fit_end_ps"] = linear_window["end_ps"]
                payload["fallback_diffusion_cm2_s_mean"] = linear_window["mean"]
                payload["fallback_diffusion_cm2_s_std"] = linear_window["std"]
                payload["fallback_fit_start_ps"] = linear_window["start_ps"]
                payload["fallback_fit_end_ps"] = linear_window["end_ps"]
                payload["fallback_quality"] = linear_window["method"]
                payload["fallback_window_score"] = linear_window["score"]
                payload["fallback_window_r2"] = linear_window["r2"]
                payload["fallback_window_slope_median"] = linear_window["slope_median"]
                payload["fallback_window_slope_iqr"] = linear_window["slope_iqr"]
                payload["fallback_window_plateau_ok"] = linear_window["plateau_ok"]
                fit_start_ps = linear_window["start_ps"]
                fit_end_ps = linear_window["end_ps"]
                fit_window_source = f"fallback_{linear_window['method']}"
                fit_window_is_fallback = True
                payload["diffusion_method"] = fit_window_source
            else:
                fallback_mean, fallback_std, fallback_start_ps, fallback_quality = robust_tail_diffusion_fit(
                    lag_ps,
                    msd_values,
                    config,
                    allow_zero_floor=True,
                    zero_method=zero_method,
                )
                if fallback_mean is not None:
                    payload["linear_diffusion_cm2_s_mean"] = fallback_mean
                    payload["linear_diffusion_cm2_s_std"] = 0.0 if fallback_std is None else fallback_std
                    payload["linear_fit_start_ps"] = fallback_start_ps
                    payload["linear_fit_end_ps"] = float(lag_ps[-1]) if len(lag_ps) else fallback_start_ps
                    payload["fallback_diffusion_cm2_s_mean"] = fallback_mean
                    payload["fallback_diffusion_cm2_s_std"] = 0.0 if fallback_std is None else fallback_std
                    payload["fallback_fit_start_ps"] = fallback_start_ps
                    payload["fallback_fit_end_ps"] = float(lag_ps[-1]) if len(lag_ps) else fallback_start_ps
                    payload["fallback_quality"] = fallback_quality or "tail_fallback"
                    fit_start_ps = fallback_start_ps
                    fit_end_ps = float(lag_ps[-1]) if len(lag_ps) else fallback_start_ps
                    fit_window_source = f"fallback_{payload['fallback_quality']}"
                    fit_window_is_fallback = True
                    payload["diffusion_method"] = fit_window_source

    fsi_fit_start = fit_start_ps if fit_start_ps is not None else max(float(lag_ps[-1]) * 0.5, 0.0)
    fsi_fit_end = fit_end_ps if fit_end_ps is not None else (float(lag_ps[-1]) if len(lag_ps) else fsi_fit_start)
    fsi_mean, fsi_std = linear_diffusion_fit(lag_ps, fsi_msd_values, fsi_fit_start, end_ps=fsi_fit_end) if len(lag_ps) else (None, None)
    if fsi_mean is not None:
        payload["fsi_diffusion_cm2_s_mean"] = fsi_mean
        payload["fsi_diffusion_cm2_s_std"] = 0.0 if fsi_std is None else fsi_std

    if payload["status"] == "continue" and completed_chunks >= int(metadata["max_chunks"]):
        payload["status"] = "max_chunks_reached"
    if fit_start_ps is not None:
        payload["fit_window_start_ps"] = fit_start_ps
        payload["fit_window_end_ps"] = fit_end_ps if fit_end_ps is not None else (float(lag_ps[-1]) if len(lag_ps) else fit_start_ps)
        payload["fit_window_source"] = fit_window_source
        payload["fit_window_is_fallback"] = fit_window_is_fallback

    # Persist analysis results before figure export so summary collection does
    # not get blocked by slow or hanging image writers on a particular case.
    write_json(out_dir / "status.json", payload)

    apply_publication_style()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(6.9, 2.9),
        gridspec_kw={"width_ratios": [1.2, 1.0]},
        constrained_layout=True,
    )

    axes[0].plot(lag_ps, msd_values, lw=1.6, color=PALETTE["blue_main"], label="Li")
    axes[0].plot(lag_ps, fsi_msd_values, lw=1.2, color=PALETTE["red_strong"], label="FSI anchor")
    plot_fit_end_ps = fit_end_ps if fit_end_ps is not None else (float(lag_ps[-1]) if len(lag_ps) else fit_start_ps)
    if fit_start_ps is not None:
        axes[0].axvspan(fit_start_ps, plot_fit_end_ps, color=PALETTE["neutral_light"], alpha=0.35, lw=0.0)
        axes[0].axvline(fit_start_ps, color=PALETTE["red_strong"], ls="--", lw=1.0)
    axes[0].set_xlabel("Lag time / ps")
    axes[0].set_ylabel(r"MSD / $\mathrm{\AA^2}$")
    axes[0].set_title(case_dir.name, fontsize=8.0, pad=4.0)
    axes[0].legend(loc="upper left", fontsize=6.2)

    axes[1].plot(lag_ps[positive_mask], msd_values[positive_mask], lw=1.2, color=PALETTE["blue_main"])
    fsi_positive_mask = (lag_ps > 0.0) & (fsi_msd_values > 0.0)
    axes[1].plot(lag_ps[fsi_positive_mask], fsi_msd_values[fsi_positive_mask], lw=1.0, color=PALETTE["red_strong"], alpha=0.9)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes_slope = axes[1].twinx()
    slope_mask = positive_mask & (~np.isnan(loglog_slope)) & (loglog_slope > 0.0)
    axes_slope.plot(
        lag_ps[slope_mask],
        loglog_slope[slope_mask],
        color=PALETTE["teal"],
        alpha=0.95,
        lw=1.0,
    )
    axes_slope.axhline(float(config["analysis"]["min_loglog_slope"]), color=PALETTE["neutral_mid"], ls=":", lw=0.8)
    axes_slope.axhline(float(config["analysis"]["max_loglog_slope"]), color=PALETTE["neutral_mid"], ls=":", lw=0.8)
    if fit_start_ps is not None:
        axes[1].axvspan(fit_start_ps, plot_fit_end_ps, color=PALETTE["neutral_light"], alpha=0.25, lw=0.0)
        axes[1].axvline(fit_start_ps, color=PALETTE["red_strong"], ls="--", lw=1.0)
        axes_slope.axvline(fit_start_ps, color=PALETTE["red_strong"], ls="--", lw=1.0)
    axes[1].set_xlabel("Lag time / ps")
    axes[1].set_ylabel(r"Li MSD / $\mathrm{\AA^2}$")
    axes_slope.set_ylabel("log-log slope")
    finite_positive_slopes = loglog_slope[np.isfinite(loglog_slope) & (loglog_slope > 0.0)]
    if finite_positive_slopes.size:
        slope_upper = float(np.quantile(finite_positive_slopes, 0.99))
        slope_upper = max(1.6, min(2.5, slope_upper * 1.1))
    else:
        slope_upper = 1.6
    axes_slope.set_ylim(0.0, slope_upper)

    summary_lines = [
        f"status: {payload['status']}",
        f"sampled: {total_ns:.2f} ns",
        f"LiFSI: {metadata['actual_pairs']} pairs",
        f"thermo eq: {payload.get('thermo_equilibrated', False)}",
    ]
    if fit_start_ps is not None:
        summary_lines.append(f"fit window: {fit_start_ps:.0f}-{plot_fit_end_ps:.0f} ps")
        summary_lines.append(f"fit source: {payload.get('fit_window_source', payload.get('diffusion_method', ''))}")
    if start_ps is not None:
        summary_lines.append(f"D-app plateau: {payload.get('d_plateau_ok', False)}")
    if "kinisi_diffusion_cm2_s_mean" in payload:
        summary_lines.append(f"D(Li) = {payload['kinisi_diffusion_cm2_s_mean']:.2e}")
    if "fsi_diffusion_cm2_s_mean" in payload:
        summary_lines.append(f"D(FSI) = {payload['fsi_diffusion_cm2_s_mean']:.2e}")
    axes[0].text(
        0.03,
        0.97,
        "\n".join(summary_lines),
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        color=PALETTE["neutral_black"],
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": PALETTE["neutral_light"], "lw": 0.6},
    )

    axes[1].text(
        0.03,
        0.96,
        "MSD",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        color=PALETTE["blue_main"],
    )
    axes_slope.text(
        0.97,
        0.08,
        "slope",
        transform=axes_slope.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=PALETTE["teal"],
    )

    add_panel_label(axes[0], "a")
    add_panel_label(axes[1], "b")
    save_publication_bundle(fig, msd_base)
    plt.close(fig)


if __name__ == "__main__":
    main()
