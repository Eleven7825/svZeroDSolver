"""Plan C, step 2: build the Eberth BANDED group on its OWN terms.

Own-terms inputs (Tier 2, measured for the banded group):
  * HR = 6.09 Hz
  * remodeled carotids (RCCA-B 591/88.1, LCCA-B 410/41.6 @ MAP)
  * band present
  * terminal-bed resistances PINNED to each carotid's measured Q̄ and MAP
    (RCCA-B: 0.022 ml/s @ MAP 86.5; LCCA-B: 0.012 @ 76.8; systemic @ 76.8)
Calibrated knobs (Tier 3):
  * band stenosis_coefficient S -> the A->B mean pressure drop (MAP_A - MAP_B ~ 10)
  * RCCA-B bed C -> PI 3.11 ; LCCA-B bed C -> PI 1.65   (flow-admittance knobs)
Aortic compliance kept at the measured value. Template: the chronic config
(already HR 6.09, remodeled carotids, band).
"""
import json, copy, pysvzerod

MAP_A, MAP_B = 86.5, 76.8       # node A (RCCA-B) and node B (LCCA-B) MAP, Table 1
Q_R, Q_L, CO = 0.022, 0.012, 0.20
OUT="/home/shiyi/projects/svZeroDSolver/examples/eberth_banded_group.json"
chr=json.load(open("/home/shiyi/projects/svZeroDSolver/examples/eberth_aortic_band_chronic.json"))

R_R=MAP_A/Q_R; R_L=MAP_B/Q_L; R_S=MAP_B/(CO-Q_R-Q_L)   # pinned bed resistances

def build(S, C_R, C_L):
    c=copy.deepcopy(chr)
    for b in c["boundary_conditions"]:
        v=b["bc_values"]
        if b["bc_name"]=="RCR_RIGHT": v["Rp"]=round(0.1*R_R,2); v["Rd"]=round(0.9*R_R,2); v["C"]=round(C_R,10)
        if b["bc_name"]=="RCR_LEFT":  v["Rp"]=round(0.1*R_L,2); v["Rd"]=round(0.9*R_L,2); v["C"]=round(C_L,10)
        if b["bc_name"]=="RCR_SYS":   v["Rp"]=round(0.1*R_S,2); v["Rd"]=round(0.9*R_S,2)
    for vv in c["vessels"]:
        if vv["vessel_name"]=="aortic_band": vv["zero_d_element_values"]["stenosis_coefficient"]=round(S,3)
    return c

def sim(cfg):
    df=pysvzerod.simulate(cfg)
    def st(s):
        d=df[df.name==s]; q=d.flow_in; p=d.pressure_in
        return dict(Q=q.mean(), PI=(q.max()-q.min())/q.mean(), MAP=p.mean(), PP=p.max()-p.min())
    return st("rcca_b"), st("lcca_b")

def bisect(f, target, lo, hi, incr, n=28):
    """find x in [lo,hi] with f(x)=target; incr=True if f increases with x."""
    for _ in range(n):
        mid=0.5*(lo+hi); v=f(mid)
        if (v>target)==incr: hi=mid
        else: lo=mid
    return mid

if __name__=="__main__":
    C_R0=1.06e-6; C_L0=0.99e-6   # nominal remodeled-carotid bed C (starting point)
    # 1) band S -> A-B MAP drop ~ (MAP_A - MAP_B)
    dropT=MAP_A-MAP_B
    S=bisect(lambda S: (lambda r,l:r["MAP"]-l["MAP"])(*sim(build(S,C_R0,C_L0))), dropT, 20, 600, incr=True)
    # 2) RCCA-B bed C -> PI 3.11 (higher C -> higher PI)
    C_R=bisect(lambda C: sim(build(S,C,C_L0))[0]["PI"], 3.11, 1e-7, 2e-4, incr=True)
    # 3) LCCA-B bed C -> PI 1.65
    C_L=bisect(lambda C: sim(build(S,C_R,C))[1]["PI"], 1.65, 1e-8, 2e-4, incr=True)
    r,l=sim(build(S,C_R,C_L))
    print(f"calibrated: band S={S:.1f}, RCCA bed C={C_R:.2e}, LCCA bed C={C_L:.2e}")
    print(f"{'':10}{'PI':>6}{'Q':>8}{'MAP':>6}{'PP':>6}   (targets)")
    print(f"RCCA-B    {r['PI']:>6.2f}{r['Q']:>8.4f}{r['MAP']:>6.0f}{r['PP']:>6.0f}   (3.11 / 0.022 / 86 / 56)")
    print(f"LCCA-B    {l['PI']:>6.2f}{l['Q']:>8.4f}{l['MAP']:>6.0f}{l['PP']:>6.0f}   (1.65 / 0.012 / 77 / 27)")
    print(f"band A-B MAP drop = {r['MAP']-l['MAP']:.1f} (target {dropT:.1f})")
    cfg=build(S,C_R,C_L)
    cfg["description"]={"model":"Eberth BANDED group, own terms (HR 6.09; remodeled carotids; beds pinned to measured RCCA-B/LCCA-B Q & MAP; band S -> A-B MAP drop; per-carotid bed C -> PIs)",
                        "calibrated":{"band_S":round(S,3),"RCCA_bed_C":round(C_R,10),"LCCA_bed_C":round(C_L,10)},
                        "result":f"RCCA-B PI {r['PI']:.2f}/Q {r['Q']:.4f}; LCCA-B PI {l['PI']:.2f}/Q {l['Q']:.4f}",
                        "see":"two_group_calibration.md"}
    json.dump(cfg,open(OUT,"w"),indent=2); print("wrote",OUT)
