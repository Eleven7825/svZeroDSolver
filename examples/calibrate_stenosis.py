"""Calibrate the aortic-band `stenosis_coefficient` in the Eberth 0D model
against the Table-1 pulsatility-index (PI) targets.

PI is velocity-based in the paper, PI = (v_max - v_min)/v_mean; for a fixed
cross-section that equals the flow-based PI = (Q_max - Q_min)/Q_mean, which is
what we compute from the carotid segment flow.

Usage:
    uv run python examples/calibrate_stenosis.py [config.json]
Default config: examples/eberth_aortic_band_chronic.json (the chronic /
remodeled-geometry state that the Table-1 PI values correspond to).
"""
import sys, json, copy
import numpy as np
from scipy.optimize import minimize_scalar
import pysvzerod

REPO = "/home/shiyi/projects/svZeroDSolver"
CFG = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/examples/eberth_aortic_band_chronic.json"

# Eberth Table 1 targets (velocity PI)
TARGET = {"rcca_b": 3.11, "lcca_b": 1.65}
TARGET_RATIO = TARGET["rcca_b"] / TARGET["lcca_b"]   # 1.885 : the pulsatility SPLIT
MAP_MAX = 130.0                                       # physiological guard [mmHg]
BAND_VESSEL_NAME = "aortic_band"

# NOTE on the objective: a single stenosis_coefficient cannot match BOTH absolute
# PIs - increasing S raises the upstream (RCCA-B) PI via wave reflection while
# lowering the downstream (LCCA-B) PI, and forcing the absolute values drives MAP
# non-physiological. The band's physical role is the RELATIVE split, so we
# calibrate S to the PI RATIO (with a MAP guard). The absolute-magnitude offset is
# set by the ejection waveform / total compliance (separate free knobs).

with open(CFG) as f:
    base = json.load(f)

def band_idx(cfg):
    for i, v in enumerate(cfg["vessels"]):
        if v["vessel_name"] == BAND_VESSEL_NAME:
            return i
    raise KeyError(BAND_VESSEL_NAME)

def flow_PI(df, seg):
    d = df[df.name == seg]
    q = d.flow_in
    return (q.max() - q.min()) / q.mean()

def simulate_with_S(S):
    cfg = copy.deepcopy(base)
    cfg["vessels"][band_idx(cfg)]["zero_d_element_values"]["stenosis_coefficient"] = float(S)
    df = pysvzerod.simulate(cfg)
    return {seg: flow_PI(df, seg) for seg in TARGET}, df

def map_upstream(df):
    return df[df.name == "rcca_b"].pressure_in.mean()

def objective(S):
    pis, df = simulate_with_S(S)
    ratio = pis["rcca_b"] / pis["lcca_b"]
    err = (ratio - TARGET_RATIO) ** 2
    if map_upstream(df) > MAP_MAX:           # penalize non-physiological MAP
        err += (map_upstream(df) - MAP_MAX) ** 2
    return err

if __name__ == "__main__":
    # 1) coarse scan for context
    print("Scan of stenosis_coefficient S vs flow-PI:")
    print(f"  {'S':>8}  {'PI_RCCA-B':>10}  {'PI_LCCA-B':>10}  {'ratio':>6}  {'objective':>10}")
    for S in [0, 10, 25, 50, 100, 200, 400, 800, 1600]:
        pis, _ = simulate_with_S(S)
        r = pis['rcca_b'] / pis['lcca_b'] if pis['lcca_b'] else float('nan')
        print(f"  {S:>8.0f}  {pis['rcca_b']:>10.2f}  {pis['lcca_b']:>10.2f}  {r:>6.2f}  {objective(S):>10.4f}")

    # 2) bounded 1-D minimization on the PI-ratio objective
    res = minimize_scalar(objective, bounds=(0.0, 1000.0), method="bounded",
                          options={"xatol": 1e-2})
    Sopt = res.x
    pis, df = simulate_with_S(Sopt)
    ratio = pis["rcca_b"] / pis["lcca_b"]
    print(f"\nBest-fit stenosis_coefficient S* = {Sopt:.1f}  (objective {res.fun:.4f})")
    print(f"PI ratio: target {TARGET_RATIO:.3f}  ->  model {ratio:.3f}")
    print(f"Absolute PI (offset by baseline pulsatility): "
          f"RCCA-B {pis['rcca_b']:.2f} (target 3.11), LCCA-B {pis['lcca_b']:.2f} (target 1.65)")
    for seg, lab in [("rcca_b","RCCA-B"),("lcca_b","LCCA-B")]:
        d = df[df.name == seg]; p = d.pressure_in; q = d.flow_in
        print(f"   {lab}: MAP={p.mean():.1f} mmHg  pulseP={p.max()-p.min():.1f}  Qmean={q.mean():.4f} ml/s")
    print(f"\n==> set {BAND_VESSEL_NAME}.stenosis_coefficient = {Sopt:.1f}")

    # 3) write the calibrated value back into the config in place
    if "--write" in sys.argv:
        base["vessels"][band_idx(base)]["zero_d_element_values"]["stenosis_coefficient"] = round(Sopt, 1)
        base.setdefault("description", {})["calibrated_stenosis_coefficient"] = {
            "value": round(Sopt, 1), "objective": "PI ratio 3.11/1.65 = 1.885 (MAP-guarded)",
            "model_ratio": round(ratio, 3),
        }
        with open(CFG, "w") as f:
            json.dump(base, f, indent=2)
        print(f"WROTE calibrated S into {CFG}")
