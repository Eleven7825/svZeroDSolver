"""Create the REMODELED-geometry variant of the Eberth aortic-band model.
Starts from the acute model (identical carotids, band as only difference) and
gives each carotid SEGMENT the Table-1 remodeled geometry, deriving R, L, C
from the note's Eq. (1)-(3):
    R = 8 mu l /(pi r^4)      ->  R  ~ 1/r^4
    L = rho l /(pi r^2)       ->  L  ~ 1/r^2
    C = 3 pi r^3 l /(2 E h)   ->  C  ~ r^3 / (E h)
Expressed as factors relative to the baseline CCA (E assumed constant across
vessels; see caveat). This isolates the effect of the measured wall remodeling."""
import json, copy

BASE = "/home/shiyi/projects/svZeroDSolver/examples/eberth_aortic_band_pulsatile.json"
OUT  = "/home/shiyi/projects/svZeroDSolver/examples/eberth_aortic_band_remodeled.json"

# Table 1 of Eberth et al. (2009), values at 100 mmHg
#            inner diameter 2ri [um],  wall thickness h [um]
geom = {
    "CCA":    {"d": 496.0, "h": 24.8},   # baseline (reference)
    "RCCA_B": {"d": 633.0, "h": 88.1},   # right, upstream of band  (remodeled)
    "LCCA_B": {"d": 482.0, "h": 41.6},   # left,  downstream of band (remodeled)
}
r = {k: v["d"] / 2.0 for k, v in geom.items()}
h = {k: v["h"] for k, v in geom.items()}

def factors(seg):
    rr = r[seg] / r["CCA"]
    hh = h[seg] / h["CCA"]
    return {
        "R": (1.0 / rr) ** 4,          # R ~ 1/r^4
        "L": (1.0 / rr) ** 2,          # L ~ 1/r^2
        "C": (rr ** 3) / hh,           # C ~ r^3 / h   (E assumed constant)
    }

fR = factors("RCCA_B")
fL = factors("LCCA_B")

# Illustrative baseline carotid-SEGMENT wall properties (same scale as the acute
# model's conduit R=1.0). These represent a normal CCA segment; the factors above
# turn them into the remodeled RCCA-B / LCCA-B segments.
R0, C0, L0 = 1.0, 2.0e-4, 0.02

def seg_vals(f):
    return {
        "R_poiseuille": round(R0 * f["R"], 5),
        "C":            round(C0 * f["C"], 8),
        "L":            round(L0 * f["L"], 5),
    }

with open(BASE) as fp:
    model = json.load(fp)

# vessel_id 2 = rcca_b (RCCA-B, node A), vessel_id 4 = lcca_b (LCCA-B, node B)
for v in model["vessels"]:
    if v["vessel_name"] == "rcca_b":
        v["zero_d_element_values"].update(seg_vals(fR))
    elif v["vessel_name"] == "lcca_b":
        v["zero_d_element_values"].update(seg_vals(fL))

model["description"] = {
    "model": "Eberth (2009) aortic banding - REMODELED-geometry variant (chronic 5-8 week state)",
    "reference": "Eberth et al., J. Hypertens. 27:2010-2021 (2009). Geometry from Table 1 (values at 100 mmHg).",
    "difference_from_acute_model": "The RCCA-B and LCCA-B carotid SEGMENTS now carry the measured "
        "remodeled wall geometry (R, C, L from Eq. 1-3), so the two carotids differ by band position "
        "PLUS wall remodeling. The acute model (eberth_aortic_band_pulsatile.json) instead keeps the "
        "carotids identical (band as the only difference).",
    "timing_note": "Eberth measured hemodynamics AND geometry 5 or 8 weeks after banding, chosen to give "
        "a 'nearly steady state adaptive response' - so Table-1 PI (RCCA-B 3.11, LCCA-B 1.65) and the "
        "geometry below are the CHRONIC, post-remodeling values, not the acute post-op state.",
    "Table1_geometry_um": {
        "inner_diameter_2ri": {"CCA": 496, "RCCA_B": 633, "LCCA_B": 482},
        "wall_thickness_h":   {"CCA": 24.8, "RCCA_B": 88.1, "LCCA_B": 41.6},
    },
    "scaling_from_Eq_1_3": {
        "formulas": "R~1/r^4, L~1/r^2, C~r^3/h  (relative to baseline CCA; E assumed constant)",
        "RCCA_B_factors": {k: round(val, 4) for k, val in fR.items()},
        "LCCA_B_factors": {k: round(val, 4) for k, val in fL.items()},
        "baseline_segment": {"R_poiseuille": R0, "C": C0, "L": L0},
        "RCCA_B_segment": seg_vals(fR),
        "LCCA_B_segment": seg_vals(fL),
        "interpretation": "RCCA-B is wider -> lower R and L; both carotids are wall-thickened -> lower "
            "compliance C than baseline. The stiffer, thicker RCCA-B wall alters local pulse transmission.",
    },
    "units": "pressure [mmHg], flow [ml/s], time [s]",
    "caveat": "Illustrative: E (incremental wall stiffness) is assumed equal across vessels because the "
        "paper reports circumferential STRESS, not modulus; segment length/mu/rho set the baseline scale. "
        "The RCR terminal beds (distal head microvasculature) are left identical - remodeling here is the "
        "carotid ARTERY wall only. To match Table-1 mean-flow split (Qmean RCCA-B 0.022 > LCCA-B 0.012), "
        "additionally tune the RCR Rd. See eberth_aortic_band_pulsatile.md.",
}

with open(OUT, "w") as fp:
    json.dump(model, fp, indent=2)

print("wrote", OUT)
print("RCCA-B factors:", {k: round(v, 4) for k, v in fR.items()}, "->", seg_vals(fR))
print("LCCA-B factors:", {k: round(v, 4) for k, v in fL.items()}, "->", seg_vals(fL))
