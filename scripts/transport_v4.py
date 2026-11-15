#!/usr/bin/env python3
"""Adaptive full-trajectory lag-window screening (v4).

Window start and width are both scanned over a continuous grid (no fixed
200/400/750/1000 ps widths, no fixed 3 ns cap).  A candidate must satisfy
alpha in [0.9,1.1], positive slope, linear R2 >= 0.98 and a half-window
slope ratio in [0.5,2]; the candidate closest to alpha=1 is selected.
"""
import argparse, json, traceback
from pathlib import Path
import numpy as np
from transport_v2 import read_positions, curves, metrics, passes, temperature_diagnostic, kinisi_check
from robust_transport import conductivity_from_slope


def select_adaptive(t, y, lag_frac, max_lag, wmin, wmax, wstep, sstep):
    limit = min(float(max_lag), lag_frac * float(t[-1]))
    widths = np.arange(float(wmin), min(float(wmax), limit) + 1e-9, float(wstep))
    candidates = []
    for width in widths:
        for start in np.arange(50.0, limit - width + 1e-9, float(sstep)):
            win = (float(start), float(start + width))
            m = metrics(t, y, win)
            if not passes(m):
                continue
            mid = (win[0] + win[1]) / 2.0
            left = metrics(t, y, (win[0], mid))
            right = metrics(t, y, (mid, win[1]))
            if left.get("status") != "ok" or right.get("status") != "ok":
                continue
            if left.get("slope", 0.0) <= 0.0:
                continue
            ratio = right["slope"] / left["slope"]
            if not (0.5 <= ratio <= 2.0):
                continue
            candidates.append({
                "status": "candidate",
                "window_ps": [win[0], win[1]],
                "fit": m,
                "half_window_slope_ratio": ratio,
            })
    if not candidates:
        return {"status": "no_candidate", "n_candidates": 0, "lag_limit_ps": limit}
    best = min(candidates, key=lambda r: (abs(r["fit"]["alpha"] - 1.0), -(r["window_ps"][1] - r["window_ps"][0]), -r["fit"]["r2"]))
    best["n_candidates"] = len(candidates)
    best["lag_limit_ps"] = limit
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--kinisi", action="store_true")
    p.add_argument("--lag-frac", type=float, default=0.8)
    p.add_argument("--max-lag-ps", type=float, default=12000.0)
    p.add_argument("--width-min", type=float, default=500.0)
    p.add_argument("--width-max", type=float, default=6000.0)
    p.add_argument("--width-step", type=float, default=250.0)
    p.add_argument("--start-step", type=float, default=100.0)
    a = p.parse_args()
    case = a.case.resolve()
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    try:
        meta = json.loads((case / "metadata.json").read_text())
        li, fsi, moment, dimensions, prov = read_positions(case, meta)
        dt = prov["dt_ps"]
        t = np.arange(len(li)) * dt
        data = {}
        results = {}
        volume = float(meta["box_volume_A3"])
        temp = float(meta["temperature_K"])
        for name, x, kind in [("Li", li, "self"), ("FSI_COM", fsi, "self"), ("collective_charge", moment, "collective")]:
            y = curves(x, kind)
            data[name] = y
            r = select_adaptive(t, y, a.lag_frac, a.max_lag_ps, a.width_min, a.width_max, a.width_step, a.start_step)
            if r["status"] == "candidate":
                slope = r["fit"]["slope"]
                if name == "collective_charge":
                    r["conductivity_mS_cm"] = 10 * conductivity_from_slope(slope, volume, temp)
                else:
                    r["diffusion_cm2_s"] = slope * 1e-4 / 6
            results[name] = r
        if all(results[x]["status"] == "candidate" for x in ("Li", "FSI_COM")):
            ql = prov["li_charge_e"]
            qf = prov["fsi_charge_e"]
            factor = 10 * int(meta["actual_pairs"]) * 1.602176634e-19 ** 2 * 1e-8 / (6 * volume * 1e-30 * 1.380649e-23 * temp)
            results["nernst_einstein_mS_cm"] = factor * (ql ** 2 * results["Li"]["fit"]["slope"] + qf ** 2 * results["FSI_COM"]["fit"]["slope"])
        if a.kinisi and results["Li"]["status"] == "candidate":
            results["Li"]["kinisi_crosscheck"] = kinisi_check(li, dimensions, meta, results["Li"]["window_ps"])
        elif a.kinisi:
            results["Li"]["kinisi_crosscheck"] = {"status": "not_run_no_candidate"}
        results.update(
            case_name=meta["case_name"], base_name=meta["base_name"], dmc_per_lifsi=meta["dmc_per_lifsi"],
            replica_index=meta["replica_index"], velocity_seed=meta["velocity_seed"],
            temperature_K=temp, volume_A3=volume, provenance=prov,
            temperature_diagnostic=temperature_diagnostic(case, meta),
            method_version="v4 adaptive full-trajectory lag windows",
            scan={"lag_frac": a.lag_frac, "max_lag_ps": a.max_lag_ps, "width_min_ps": a.width_min,
                  "width_max_ps": a.width_max, "width_step_ps": a.width_step, "start_step_ps": a.start_step},
            warning="candidate is exploratory, not proof of convergence or equilibrium",
        )
        import pandas as pd
        pd.DataFrame({"lag_ps": t, **{k + "_msd": v for k, v in data.items()}}).to_csv(out / "curves.csv", index=False)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
        for i, (name, y) in enumerate(data.items()):
            r = results[name]
            mask = (t > 0) & (y > 0)
            for ax in axes[i]:
                ax.plot(t[mask], y[mask])
                ax.set_xlabel("Lag time (ps)")
                ax.set_ylabel(name + " MSD")
                ax.grid(alpha=.2)
                if "window_ps" in r:
                    ax.axvspan(*r["window_ps"], color="orange", alpha=.25)
            axes[i, 0].set_xscale("log"); axes[i, 0].set_yscale("log")
            axes[i, 0].set_title(name + " | " + r["status"] + (f" | alpha={r['fit']['alpha']:.3f}" if r["status"] == "candidate" else ""))
            if r["status"] == "candidate":
                fit = r["fit"]
                m2 = (t >= r["window_ps"][0]) & (t <= r["window_ps"][1])
                axes[i, 1].plot(t[m2], fit["slope"] * t[m2] + fit["intercept"], color="crimson")
                v = (f"sigma={r['conductivity_mS_cm']:.3g} mS/cm" if name == "collective_charge" else f"D={r['diffusion_cm2_s']:.3g} cm2/s")
                axes[i, 1].set_title(v)
            else:
                axes[i, 1].set_title("No screened estimate")
        fig.savefig(out / "diagnostics.png", dpi=160)
        fig.savefig(out / "diagnostics.pdf")
        plt.close(fig)
        tmp = out / "results.json.tmp"
        tmp.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
        tmp.replace(out / "results.json")
        print(json.dumps({k: results[k]["status"] for k in data}), flush=True)
    except Exception:
        (out / "failure.json").write_text(json.dumps({"traceback": traceback.format_exc()}, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
