"""Build the physiologically-parameterized Eberth aortic-band 0D configs.

Emits BOTH model states, which differ ONLY in carotid-wall geometry:
  * acute   -> both carotids use BASELINE CCA geometry (un-remodeled). Right
               after banding; the band is the only difference between them.
  * chronic -> carotids use the Table-1 REMODELED geometry (RCCA-B wider +
               ~3.5x thicker wall; LCCA-B thickened). The 5-8 week state Eberth
               measured, and the one the Table-1 PI targets correspond to.
Everything else (aorta, RCR beds, inflow, blood, band) is identical, so the two
files are a controlled comparison: band-alone vs band + wall remodeling.

Every vessel R, L, C is COMPUTED from researched mouse geometry + wall stiffness
via Eq. (1)-(3):
    R = 8 mu l / (pi r^4)   ;   L = rho l / (pi r^2)   ;   C = 3 pi r^3 l /(2 E h)
in CGS, then converted to the solver's mmHg / ml / s units.

The band spacer is the same physical object in both states, so the calibrated
stenosis_coefficient (S = 41.1, from calibrate_stenosis.py on the chronic
geometry / PI ratio) is used for BOTH. Parameters are tagged FIXED / DERIVED /
CALIBRATE in description.parameter_tags.

Sources: CO/HR/MAP - Janssen 2002, Constantinides 2011, Mills 2000; blood mu/rho
- Windberger 2003, Aslanidou 2016; aorta dims - Casteleyn 2010, Guo&Kassab 2002;
carotid dims - Eberth 2009 Table 1, Ferruzzi 2013; E - Ferruzzi 2013 / Wagenseil
2019; flow split - Feintuch 2007, Trachet 2009.
"""
import json, math

OUTDIR = "/home/shiyi/projects/svZeroDSolver/examples"
S_CALIB = 41.1   # calibrated band stenosis_coefficient (PI ratio, chronic geometry @MAP)

# ---- unit conversion --------------------------------------------------------
MMHG = 1333.22          # 1 mmHg in dyn/cm^2
def R_cgs(mu, l, r):   return 8.0 * mu * l / (math.pi * r**4)
def L_cgs(rho, l, r):  return rho * l / (math.pi * r**2)
def C_cgs(r, l, E, h): return 3.0 * math.pi * r**3 * l / (2.0 * E * h)
def toR(x): return x / MMHG
def toL(x): return x / MMHG
def toC(x): return x * MMHG

# ---- FIXED physiological constants (literature / Eberth 2009) ---------------
rho = 1.06            # g/cm^3    (Aslanidou 2016)
mu  = 0.035           # Poise = 3.5 cP  (Windberger 2003, high-shear)
E_ao  = 1.0e7         # dyn/cm^2 = 1.0 MPa  aorta  (Wagenseil 2019 / Bersi)
E_car = 1.5e7         # dyn/cm^2 = 1.5 MPa  carotid (Ferruzzi 2013, systolic)

HR   = 6.09           # Hz  (Eberth 2009, banded mice)
T    = 1.0 / HR       # s   cardiac period
CO   = 0.20           # ml/s = 12 ml/min  (anesthetized mouse; Aslanidou 2016)
Tsf  = 0.40           # systolic ejection fraction of the cycle (assumed)
MAP  = 90.0           # mmHg target (anesthetized; Constantinides 2011)

# ---- geometry [cm] : (radius, length, wall thickness) -----------------------
asc = dict(r=0.070, l=0.25, h=0.0040)     # ascending aorta, ID~1.4 mm
des = dict(r=0.045, l=0.80, h=0.0040)     # descending aorta, ID~0.9 mm
# carotid geometry sets from Eberth Table 1 (ID [um] @ MAP, wall [um]); length
# ~7 mm. Table 1 reports inner diameter at TWO conditions -- "@ MAP" (each
# vessel's own physiological pressure) and "@ 100 mmHg" (a common ex vivo
# reference pressure for cross-group comparison). We use "@ MAP" since this
# model runs each vessel at its own distinct physiological MAP, not a shared
# 100 mmHg -- the "@100mmHg" values (496/633/482) over-distend each vessel by
# an amount that grows with how far its own MAP sits below 100 (worst for
# LCCA-B: MAP~77, ~18% diameter inflation vs ~2.5% for CCA at MAP~92).
CAROTID = {
    #                       ID (@ MAP)   wall
    "CCA_baseline": dict(d=484.0, h=24.8),   # un-remodeled (acute, both sides)
    "RCCA_B":       dict(d=591.0, h=88.1),   # remodeled right (chronic)
    "LCCA_B":       dict(d=410.0, h=41.6),   # remodeled left  (chronic)
}
def carotid_geom(key):
    g = CAROTID[key]
    return dict(r=g["d"]*1e-4/2, l=0.70, h=g["h"]*1e-4)

# band (aortic arch stenosis): throat 406 um, local aortic lumen D0 875 um
r_throat = 406e-4 / 2
r0_band  = 875e-4 / 2
Ls       = 0.05
A0 = math.pi * r0_band**2
As = math.pi * r_throat**2
Kt = 1.5

def vessel_RLC(g, E):
    return (toR(R_cgs(mu, g["l"], g["r"])),
            toC(C_cgs(g["r"], g["l"], E, g["h"])),
            toL(L_cgs(rho, g["l"], g["r"])))

# mode-independent segments
R_asc, C_asc0, L_asc = vessel_RLC(asc, E_ao)
R_des, C_des0, L_des = vessel_RLC(des, E_ao)

# aortic buffering compliance (Ca, Cb): measured central aortic compliance
# (Aslanidou 2016) - the per-segment geometric C underestimates it ~20x.
C_aorta_central = 2.67e-4
C_asc = round(0.64 * C_aorta_central, 8)   # Ca @ node A
C_des = round(0.36 * C_aorta_central, 8)   # Cb @ node B

# band viscous R, inertial L (from throat); geometry prior for S
R_band = toR(R_cgs(mu, Ls, r_throat))
L_band = toL(rho * Ls / As)
S_prior = (Kt * rho / (2 * A0**2) * (A0/As - 1)**2) / MMHG

# normal aortic-arch segment between the two carotid takeoffs, for the
# PRE-BANDING state (no stenosis). Uses the arch lumen (D0=875 um), not the throat.
arch = dict(r=r0_band, l=0.10, h=0.0040)          # arch segment ~1 mm long
R_arch, C_arch, L_arch = vessel_RLC(arch, E_ao)

# RCR terminal beds (distal microvasculature)
def rcr(Q, tau, Pd=0.0):
    Rtot = MAP / Q
    return dict(Rp=round(0.10*Rtot, 2), C=round(tau/(0.90*Rtot), 8), Rd=round(0.90*Rtot, 2), Pd=Pd)
rcr_head = rcr(Q=0.017, tau=0.10)
rcr_sys  = rcr(Q=CO - 2*0.017, tau=0.45)

# inflow ejection waveform (half-sine during systole)
Ts = Tsf * T
Qpeak = CO / ((2.0/math.pi) * (Ts/T))
N = 81
t = [round(i*T/(N-1), 6) for i in range(N)]
Q = [round(Qpeak*math.sin(math.pi*ti/Ts), 6) if ti <= Ts else 0.0 for ti in t]

def rv(x, n=6): return round(x, n)

def build(mode):
    """mode in {'prebanding','acute','chronic'} -> (config dict, RCCA-B RLC, LCCA-B RLC)."""
    if mode == "prebanding":
        rkey, lkey = "CCA_baseline", "CCA_baseline"
        state = "PRE-BANDING (healthy baseline, no band)"
        note = ("Healthy control BEFORE surgery: carotids are baseline CCA and the "
                "segment between the two takeoffs is a NORMAL aortic arch (no "
                "stenosis, S=0). Both carotids are symmetric; this is the reference "
                "against which the acute band effect is measured. Table-1 baseline "
                "CCA PI is 1.16.")
    elif mode == "acute":
        rkey, lkey = "CCA_baseline", "CCA_baseline"
        state = "ACUTE (right after banding, before wall remodeling)"
        note = ("Carotids use BASELINE CCA geometry (un-remodeled) - acutely both "
                "carotids are still normal, so the band is the ONLY inter-carotid "
                "difference. Stenosis reuses the calibrated chronic value (same "
                "physical band); the resulting PI is a PREDICTION (no acute PI in paper).")
    elif mode == "chronic":
        rkey, lkey = "RCCA_B", "LCCA_B"
        state = "CHRONIC (5-8 weeks post-band, after wall remodeling)"
        note = ("Carotids use the Table-1 REMODELED geometry. This is the state "
                "Eberth measured; the band stenosis_coefficient was calibrated on "
                "this geometry to the Table-1 PI ratio (3.11/1.65).")
    else:
        raise ValueError(mode)

    rg, lg = carotid_geom(rkey), carotid_geom(lkey)
    R_rc, C_rc, L_rc = vessel_RLC(rg, E_car)
    R_lc, C_lc, L_lc = vessel_RLC(lg, E_car)

    # segment between node A and node B: normal arch (pre-banding) or band (post)
    if mode == "prebanding":
        seg3_name = "aortic_arch"
        seg3_vals = {"R_poiseuille": rv(R_arch, 5), "C": rv(C_arch, 9), "L": rv(L_arch, 6)}
        seg3_len = arch["l"] * 10
        tag_seg_RL = {"value": {"R": rv(R_arch, 3), "L": rv(L_arch, 4)}, "tag": "DERIVED", "src": "normal aortic-arch segment (no band)"}
        tag_seg_S  = {"value": 0.0, "tag": "FIXED", "src": "no band before surgery"}
    else:
        seg3_name = "aortic_band"
        seg3_vals = {"R_poiseuille": rv(R_band, 5), "L": rv(L_band, 6), "stenosis_coefficient": S_CALIB}
        seg3_len = Ls * 10
        tag_seg_RL = {"value": {"R": rv(R_band, 3), "L": rv(L_band, 4)}, "tag": "DERIVED", "src": "throat geometry (Young-Tsai viscous+inertial)"}
        tag_seg_S  = {"value": S_CALIB, "prior": rv(S_prior, 1), "tag": "CALIBRATE", "src": f"geometry prior Kt={Kt}; calibrated on chronic geometry to PI ratio, reused (same physical band)"}

    model = {
      "description": {
        "model": f"Eberth (2009) aortic banding - {state}, physiologically parameterized",
        "units": "pressure [mmHg], flow [ml/s], time [s]",
        "condition": f"anesthetized banded mouse; HR={HR} Hz (T={T:.4f} s), CO={CO} ml/s, MAP~{MAP} mmHg",
        "state_note": note,
        "parameter_tags": {
          "_legend": "FIXED = literature/paper measured; DERIVED = computed from geometry+E via Eq(1-3); CALIBRATE = tuned to Table-1 targets",
          "blood_density_rho_g_cm3":   {"value": rho, "tag": "FIXED", "src": "Aslanidou 2016"},
          "blood_viscosity_mu_Poise":  {"value": mu,  "tag": "FIXED", "src": "Windberger 2003 (3.5 cP high-shear)"},
          "E_aorta_MPa":   {"value": E_ao/1e7,  "tag": "FIXED", "src": "Wagenseil 2019 / Bersi"},
          "E_carotid_MPa": {"value": E_car/1e7, "tag": "FIXED", "src": "Ferruzzi 2013 (systolic linearized)"},
          "heart_rate_Hz": {"value": HR, "tag": "FIXED", "src": "Eberth 2009 (banded)"},
          "cardiac_output_ml_s": {"value": CO, "tag": "FIXED", "src": "Aslanidou 2016 (anesthetized ~12 ml/min)"},
          "MAP_target_mmHg": {"value": MAP, "tag": "FIXED", "src": "Constantinides 2011"},
          "systolic_fraction": {"value": Tsf, "tag": "CALIBRATE", "src": "assumed ejection shape"},
          "carotid_geometry_ID_wall_um": {"value": {"RCCA": [CAROTID[rkey]["d"], CAROTID[rkey]["h"]], "LCCA": [CAROTID[lkey]["d"], CAROTID[lkey]["h"]]}, "tag": "FIXED", "src": f"Eberth 2009 Table 1 ({mode}: {'baseline CCA both sides' if mode=='acute' else 'remodeled RCCA-B/LCCA-B'})"},
          "aorta_geometry": {"value": {"asc_ID_mm": 1.4, "des_ID_mm": 0.9, "wall_um": 40}, "tag": "FIXED", "src": "Casteleyn 2010, Guo&Kassab 2002"},
          "segment_lengths_cm": {"value": {"asc": asc["l"], "des": des["l"], "carotid": rg["l"], "band": Ls}, "tag": "FIXED*", "src": "carotid/band lengths are literature-guided ESTIMATES"},
          "vessel_R_L_C": {"value": "all vessel R_poiseuille/L (and carotid C) below", "tag": "DERIVED", "src": "Eq (1-3) from geometry+E"},
          "aortic_buffer_compliance_Ca_Cb_ml_mmHg": {"value": {"Ca_asc": C_asc, "Cb_des": C_des, "total": C_aorta_central}, "tag": "FIXED", "src": "Aslanidou 2016 measured central aortic compliance (overrides ~20x-smaller per-segment geometric estimate)"},
          "AB_segment_R_L": tag_seg_RL,
          "AB_segment_stenosis_coefficient": tag_seg_S,
          "RCR_beds": {"value": "Rp/C/Rd", "tag": "DERIVED", "src": "MAP/Q_target split 10/90, tau=Rd*C; Q split Feintuch 2007/Trachet 2009"},
          "valve_params": {"value": "Rmax/Rmin/Steepness", "tag": "FIXED", "src": "assumed diode (not in paper)"},
          "distal_pressure_Pd": {"value": 0.0, "tag": "FIXED", "src": "venous reference assumed 0"}
        },
        "PI_targets_Table1": {"RCCA_B": 3.11, "LCCA_B": 1.65, "CCA_baseline": 1.16},
      },
      "simulation_parameters": {
        "number_of_cardiac_cycles": 25,
        "number_of_time_pts_per_cardiac_cycle": 201,
        "output_all_cycles": False, "steady_initial": True, "absolute_tolerance": 1e-9
      },
      "boundary_conditions": [
        {"bc_name": "INFLOW", "bc_type": "FLOW", "bc_values": {"Q": Q, "t": t}},
        {"bc_name": "RCR_RIGHT", "bc_type": "RCR", "bc_values": dict(rcr_head)},
        {"bc_name": "RCR_LEFT",  "bc_type": "RCR", "bc_values": dict(rcr_head)},
        {"bc_name": "RCR_SYS",   "bc_type": "RCR", "bc_values": dict(rcr_sys)}
      ],
      "valves": [
        {"type": "ValveTanh", "name": "aortic_valve",
         "params": {"Rmax": 1.0e5, "Rmin": max(R_asc, 1.0), "Steepness": 1.0,
                    "upstream_block": "aortic_root", "downstream_block": "ascending_aorta"}}
      ],
      "junctions": [
        {"junction_name": "node_A", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [1], "outlet_vessels": [2, 3]},
        {"junction_name": "node_B", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [3], "outlet_vessels": [4, 5]}
      ],
      "vessels": [
        {"vessel_id": 0, "vessel_name": "aortic_root", "vessel_length": asc["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"inlet": "INFLOW"},
         "zero_d_element_values": {"R_poiseuille": rv(R_asc*0.3,5), "L": rv(L_asc*0.3,6)}},
        {"vessel_id": 1, "vessel_name": "ascending_aorta", "vessel_length": asc["l"]*10,
         "zero_d_element_type": "BloodVessel",
         "zero_d_element_values": {"R_poiseuille": rv(R_asc,5), "C": rv(C_asc,9), "L": rv(L_asc,6)}},
        {"vessel_id": 2, "vessel_name": "rcca_b", "vessel_length": rg["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_RIGHT"},
         "zero_d_element_values": {"R_poiseuille": rv(R_rc,5), "C": rv(C_rc,10), "L": rv(L_rc,6)}},
        {"vessel_id": 3, "vessel_name": seg3_name, "vessel_length": seg3_len,
         "zero_d_element_type": "BloodVessel",
         "zero_d_element_values": seg3_vals},
        {"vessel_id": 4, "vessel_name": "lcca_b", "vessel_length": lg["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_LEFT"},
         "zero_d_element_values": {"R_poiseuille": rv(R_lc,5), "C": rv(C_lc,10), "L": rv(L_lc,6)}},
        {"vessel_id": 5, "vessel_name": "descending_aorta", "vessel_length": des["l"]*10,
         "zero_d_element_type": "BloodVessel", "boundary_conditions": {"outlet": "RCR_SYS"},
         "zero_d_element_values": {"R_poiseuille": rv(R_des,5), "C": rv(C_des,9), "L": rv(L_des,6)}}
      ]
    }
    return model, (R_rc, C_rc, L_rc), (R_lc, C_lc, L_lc)

if __name__ == "__main__":
    for mode in ("prebanding", "acute", "chronic"):
        model, rc, lc = build(mode)
        out = f"{OUTDIR}/eberth_aortic_band_{mode}.json"
        with open(out, "w") as f:
            json.dump(model, f, indent=2)
        print(f"wrote {out}")
        print(f"   RCCA-B carotid: R={rc[0]:.2f} C={rc[1]:.3e} L={rc[2]:.4f}")
        print(f"   LCCA-B carotid: R={lc[0]:.2f} C={lc[1]:.3e} L={lc[2]:.4f}")
    print(f"\nBoth use stenosis_coefficient = {S_CALIB} (geometry prior {S_prior:.0f}), "
          f"band R={R_band:.2f}, L={L_band:.4f}")
