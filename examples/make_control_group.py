"""Plan C, step 1: build the Eberth CONTROL group on its OWN terms.

Rationale: control (CCA) and banded (RCCA-B/LCCA-B) are two separate animal
groups, so they need not share R/C/L. Here the control is built from the
control group's own measured quantities:
  * HR = 7.17 Hz (control; banded is 6.09)          [Tier 2, measured]
  * carotid geometry = baseline CCA 496/24.8 um       [Tier 2, measured]
  * terminal beds PINNED to the measured CCA Qbar=0.016 ml/s and MAP~92 mmHg
  * no band (normal arch)
  * aortic compliance = the free/uncertain knob; we test the measured value
    and also what value would be needed to force the flow-PI to 1.16.

Everything else (blood, topology, ejection shape, aortic/carotid E) is a shared
species invariant (Tier 1), taken from make_configs.
"""
import json, copy, math, pysvzerod
import make_configs as mc

HR=7.17; T=1.0/HR; CO=0.20; MAP=92.0; Qcar=0.016; Tsf=0.40
C_AORTA_MEAS=2.67e-4   # measured central aortic compliance (Aslanidou 2016)
OUT="/home/shiyi/projects/svZeroDSolver/examples/eberth_control_group.json"

def rcr(Q,tau):
    R=MAP/Q
    return dict(Rp=round(0.1*R,2), C=round(tau/(0.9*R),8), Rd=round(0.9*R,2), Pd=0.0)

def inflow():
    Ts=Tsf*T; Qpk=CO/((2/math.pi)*(Ts/T)); N=81
    t=[round(i*T/(N-1),6) for i in range(N)]
    Q=[round(Qpk*math.sin(math.pi*ti/Ts),6) if ti<=Ts else 0.0 for ti in t]
    return {"Q":Q,"t":t}

def build(Caorta):
    c=json.load(open("/home/shiyi/projects/svZeroDSolver/examples/eberth_aortic_band_prebanding.json"))
    for b in c["boundary_conditions"]:
        if b["bc_name"]=="INFLOW": b["bc_values"]=inflow()
        if b["bc_name"] in ("RCR_RIGHT","RCR_LEFT"): b["bc_values"]=rcr(Qcar,0.10)
        if b["bc_name"]=="RCR_SYS": b["bc_values"]=rcr(CO-2*Qcar,0.45)
    for v in c["vessels"]:
        if v["vessel_name"]=="ascending_aorta":  v["zero_d_element_values"]["C"]=round(0.64*Caorta,9)
        if v["vessel_name"]=="descending_aorta": v["zero_d_element_values"]["C"]=round(0.36*Caorta,9)
    return c

def stats(cfg):
    df=pysvzerod.simulate(cfg); d=df[df.name=="rcca_b"]; q=d.flow_in; p=d.pressure_in
    return dict(PI=(q.max()-q.min())/q.mean(), Qm=q.mean(), MAP=p.mean(), PP=p.max()-p.min())

if __name__=="__main__":
    s=stats(build(C_AORTA_MEAS))
    print(f"control @ measured aortic C={C_AORTA_MEAS:.2e}: "
          f"flow-PI={s['PI']:.2f}  Qbar={s['Qm']:.4f}  MAP={s['MAP']:.0f}  PP={s['PP']:.0f}")
    print("  (paper CCA: PI 1.16, Qbar 0.016, MAP ~92, PP 42.5)")
    # what aortic C would force flow-PI to 1.16?
    lo,hi=1e-4,3e-3
    for _ in range(34):
        mid=0.5*(lo+hi)
        if stats(build(mid))["PI"]>1.16: lo=mid
        else: hi=mid
    s2=stats(build(mid))
    print(f"to force flow-PI=1.16 needs aortic C={mid:.2e} (x{mid/C_AORTA_MEAS:.1f} measured) -> PP={s2['PP']:.0f} (measured 42.5)")
    # save the physiological-compliance control config (matches pressure PP; flow-PI over-predicted)
    cfg=build(C_AORTA_MEAS)
    cfg["description"]={"model":"Eberth CONTROL group, own terms (HR 7.17, beds pinned to CCA Qbar 0.016 & MAP 92, baseline carotids, no band, measured aortic compliance)",
                        "result":f"Qbar {s['Qm']:.4f}, MAP {s['MAP']:.0f}, PP {s['PP']:.0f} (match measured); flow-PI {s['PI']:.2f} over-predicts velocity PI 1.16",
                        "see":"two_group_calibration.md"}
    json.dump(cfg,open(OUT,"w"),indent=2)
    print("wrote",OUT)
