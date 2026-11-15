#!/usr/bin/env python3
"""Conservative full-trajectory lag-window screening (v5).

Predeclared lag windows, first eligible window, adjacent-window checks.

"""
import argparse, json, traceback
from pathlib import Path
import numpy as np
from transport_v2 import read_positions, curves, metrics, passes, temperature_diagnostic, kinisi_check
from robust_transport import conductivity_from_slope


def select_stable(t, y):
    t=np.asarray(t,float);y=np.asarray(y,float)
    if len(t)!=len(y) or len(t)<2 or not np.isfinite(t).all() or not np.isfinite(y).all() or np.any(np.diff(t)<=0):return {'status':'invalid_curve'}
    limit=.5*float(t[-1]); candidates=[];n_tested=0
    # Fixed ordering: longest predefined width then earliest start, not the best-looking alpha.
    for width in (2000.,1000.,500.):
        shift=100.
        for start in np.arange(150.,limit-width-shift+1,100.):
            win=(float(start),float(start+width));n_tested+=1;m=metrics(t,y,win)
            if not passes(m):continue
            mid=(win[0]+win[1])/2;left=metrics(t,y,(win[0],mid));right=metrics(t,y,(mid,win[1]))
            if left.get('status')!='ok' or right.get('status')!='ok' or left['slope']<=0:continue
            ratio=right['slope']/left['slope']
            if not (2/3<=ratio<=1.5):continue
            before=metrics(t,y,(win[0]-shift,win[1]-shift));after=metrics(t,y,(win[0]+shift,win[1]+shift))
            if before.get('status')!='ok' or after.get('status')!='ok':continue
            sensitivity=max(abs(before['slope']/m['slope']-1),abs(after['slope']/m['slope']-1))
            if sensitivity>.35:continue
            candidates.append({'status':'candidate','window_ps':list(win),'fit':m,'half_window_slope_ratio':ratio,'shift_slope_max_relative_difference':sensitivity})
    if not candidates:return {'status':'no_candidate','n_candidates':0,'n_tested':n_tested,'lag_limit_ps':limit}
    best=candidates[0];best.update(n_candidates=len(candidates),n_tested=n_tested,lag_limit_ps=limit)
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--kinisi", action="store_true")
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
            r = select_stable(t, y)
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
            method_version="v5 unified conservative full-trajectory lag screening",
            scan={"lag_fraction": 0.5, "widths_ps": [2000,1000,500], "start_step_ps": 100, "shift_ps": 100, "max_shift_slope_difference": 0.35},
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
