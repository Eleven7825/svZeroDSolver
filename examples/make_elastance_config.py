"""EXPERIMENTAL: Eberth aortic-band model driven by a CONTRACTILE
time-varying-elastance left ventricle instead of a prescribed-flow inlet.

Topology:
  preload (LA pressure) -> mitral ValveTanh -> ChamberElastanceInductor
  (ventricle) -> aortic ValveTanh -> ascending aorta -> node A -> band -> node B
  -> descending aorta ; carotid RCR beds off nodes A and B.

Reuses the physiological vessel / bed / band constants from make_configs.py; the
ONLY new part is the heart. Goal: test whether heart-afterload coupling amplifies
the banded upstream-carotid pulsatility (the 1.16 -> 3.11 gap the prescribed-flow
model cannot reach).

LV parameters are placeholders until the literature lookup returns; units are the
solver's mmHg / ml / s  (volumes in ml, elastance in mmHg/ml).
"""
import json, copy
import make_configs as mc   # committed builder: reuse aorta/carotid/band/bed values

OUTDIR = "/home/shiyi/projects/svZeroDSolver/examples"

# ---- LV time-varying elastance (mouse literature; ml, mmHg, s) --------------
# Ees~8-10 mmHg/uL (Pacher 2008 Nat Protoc; Sci Rep 2019 sham 8.4); EDV~45,
# ESV~18, SV~27 uL, ESP~105, LVEDP~5, V0~5 uL, filling~5 mmHg. Same heart is
# used for all states; banding raises afterload via the band alone.
LV = dict(
    Emax=9000.0,      # mmHg/ml  (= 9 mmHg/uL)   end-systolic elastance (Ees)
    Emin=150.0,       # mmHg/ml  (= 0.15 mmHg/uL) diastolic elastance -> LVEDP~5
    Vrd=0.008,        # ml  rest diastolic volume
    Vrs=0.005,        # ml  rest systolic volume (ESPVR intercept V0)
    Impedance=5.0e-4, # inductor (Pacher-scale)
)
V0_EDV   = 0.045      # ml  initial chamber volume (~EDV, 45 uL)
PRELOAD  = 10.0       # mmHg  filling pressure (raised so the stiff mouse LV fills
                      #       in the short diastole -> EDV~52, SV~33 uL, CO~0.2)
MITRAL_RMIN = 0.05    # low mitral resistance so filling is not rate-limited
T        = mc.T       # cardiac period (from make_configs, HR 6.09 Hz)
T_ACTIVE = 0.010      # s  delay before contraction (filling)
T_TWITCH = 0.070      # s  twitch duration (~systole at mouse HR)

def rv(x, n=6): return round(x, n)

def build(mode):
    # carotid geometry + segment between nodes A/B, exactly as make_configs
    if mode == "prebanding":
        rkey = lkey = "CCA_baseline"; seg = "arch"
    elif mode == "acute":
        rkey = lkey = "CCA_baseline"; seg = "band"
    elif mode == "chronic":
        rkey, lkey = "RCCA_B", "LCCA_B"; seg = "band"
    else:
        raise ValueError(mode)
    rg, lg = mc.carotid_geom(rkey), mc.carotid_geom(lkey)
    R_rc, C_rc, L_rc = mc.vessel_RLC(rg, mc.E_car)
    R_lc, C_lc, L_lc = mc.vessel_RLC(lg, mc.E_car)

    if seg == "arch":
        seg3_name = "aortic_arch"
        seg3_vals = {"R_poiseuille": rv(mc.R_arch,5), "C": rv(mc.C_arch,9), "L": rv(mc.L_arch,6)}
        seg3_len = mc.arch["l"]*10
    else:
        seg3_name = "aortic_band"
        seg3_vals = {"R_poiseuille": rv(mc.R_band,5), "L": rv(mc.L_band,6),
                     "stenosis_coefficient": mc.S_CALIB}
        seg3_len = mc.Ls*10

    model = {
      "description": {
        "model": f"Eberth aortic banding - {mode} - ELASTANCE HEART (experimental)",
        "units": "pressure [mmHg], flow [ml/s], volume [ml], time [s]",
        "heart": "ChamberElastanceInductor ventricle + mitral/aortic ValveTanh; "
                 "prescribed-flow inlet replaced by a contractile pump.",
        "LV_params": LV, "preload_mmHg": PRELOAD,
      },
      "simulation_parameters": {
        "number_of_cardiac_cycles": 30,
        "number_of_time_pts_per_cardiac_cycle": 201,
        "output_all_cycles": False,
        "steady_initial": False,
        "cardiac_period": T,
        "absolute_tolerance": 1e-9
      },
      "boundary_conditions": [
        {"bc_name": "PRELOAD", "bc_type": "PRESSURE",
         "bc_values": {"P": [PRELOAD, PRELOAD], "t": [0.0, T]}},
        {"bc_name": "RCR_RIGHT", "bc_type": "RCR", "bc_values": dict(mc.rcr_head)},
        {"bc_name": "RCR_LEFT",  "bc_type": "RCR", "bc_values": dict(mc.rcr_head)},
        {"bc_name": "RCR_SYS",   "bc_type": "RCR", "bc_values": dict(mc.rcr_sys)}
      ],
      "chambers": [
        {"type": "ChamberElastanceInductor", "name": "ventricle",
         "values": {"Emax": LV["Emax"], "Emin": LV["Emin"], "Vrd": LV["Vrd"],
                    "Vrs": LV["Vrs"], "Impedance": LV["Impedance"]},
         "activation_function": {"type": "half_cosine",
                                 "t_active": T_ACTIVE, "t_twitch": T_TWITCH}}
      ],
      "valves": [
        {"type": "ValveTanh", "name": "mitral_valve",
         "params": {"Rmax": 1.0e5, "Rmin": MITRAL_RMIN, "Steepness": 1.0,
                    "upstream_block": "PRELOAD", "downstream_block": "ventricle"}},
        {"type": "ValveTanh", "name": "aortic_valve",
         "params": {"Rmax": 1.0e5, "Rmin": 1.0, "Steepness": 1.0,
                    "upstream_block": "ventricle", "downstream_block": "ascending_aorta"}}
      ],
      "junctions": [
        {"junction_name": "node_A", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [0], "outlet_vessels": [1, 2]},
        {"junction_name": "node_B", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [2], "outlet_vessels": [3, 4]}
      ],
      "vessels": [
        {"vessel_id": 0, "vessel_name": "ascending_aorta", "vessel_length": mc.asc["l"]*10,
         "zero_d_element_type": "BloodVessel",
         "zero_d_element_values": {"R_poiseuille": rv(mc.R_asc,5), "C": rv(mc.C_asc,9), "L": rv(mc.L_asc,6)}},
        {"vessel_id": 1, "vessel_name": "rcca_b", "vessel_length": rg["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_RIGHT"},
         "zero_d_element_values": {"R_poiseuille": rv(R_rc,5), "C": rv(C_rc,10), "L": rv(L_rc,6)}},
        {"vessel_id": 2, "vessel_name": seg3_name, "vessel_length": seg3_len,
         "zero_d_element_type": "BloodVessel", "zero_d_element_values": seg3_vals},
        {"vessel_id": 3, "vessel_name": "lcca_b", "vessel_length": lg["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_LEFT"},
         "zero_d_element_values": {"R_poiseuille": rv(R_lc,5), "C": rv(C_lc,10), "L": rv(L_lc,6)}},
        {"vessel_id": 4, "vessel_name": "descending_aorta", "vessel_length": mc.des["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_SYS"},
         "zero_d_element_values": {"R_poiseuille": rv(mc.R_des,5), "C": rv(mc.C_des,9), "L": rv(mc.L_des,6)}}
      ],
      "initial_condition": {"Vc:ventricle": V0_EDV}
    }
    return model

if __name__ == "__main__":
    for mode in ("prebanding", "acute", "chronic"):
        out = f"{OUTDIR}/eberth_elastance_{mode}.json"
        with open(out, "w") as f:
            json.dump(build(mode), f, indent=2)
        print("wrote", out)
    print(f"LV: {LV}  preload {PRELOAD} mmHg  T {T:.4f}s  t_active {T_ACTIVE} t_twitch {T_TWITCH}")
