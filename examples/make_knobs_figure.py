"""Settled calibration knobs: terminal-bed resistance R and compliance C for the
three carotid states (CCA control, RCCA-B, LCCA-B). Same bar style as the
two-group comparison figure.
"""
import json
import numpy as np
import matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

REPO="/home/shiyi/projects/svZeroDSolver"; OUT=f"{REPO}/examples/figures"
MUTED, BLUE, GREEN, INK, INK2, GRID = "#8a8a86", "#2a78d6", "#008300", "#0b0b0b", "#52514e", "#e6e6e2"
mpl.rcParams.update({"figure.dpi":200,"savefig.dpi":200,"savefig.bbox":"tight","font.size":10.5,
    "font.family":"sans-serif","axes.edgecolor":MUTED,"axes.labelcolor":INK,"text.color":INK,
    "xtick.color":INK2,"ytick.color":INK2,"axes.linewidth":0.8,"axes.grid":True,"grid.color":GRID,
    "grid.linewidth":0.8,"axes.axisbelow":True,"axes.spines.top":False,"axes.spines.right":False})

def bed(cfg_path, bc_name):
    c=json.load(open(cfg_path))
    v=[b["bc_values"] for b in c["boundary_conditions"] if b["bc_name"]==bc_name][0]
    return v["Rp"]+v["Rd"], v["C"]
ctrl=f"{REPO}/examples/eberth_control_group.json"; band=f"{REPO}/examples/eberth_banded_group.json"
R_cca,C_cca = bed(ctrl,"RCR_RIGHT")     # control: both carotids = CCA
R_rb, C_rb  = bed(band,"RCR_RIGHT")     # banded right  = RCCA-B
R_lb, C_lb  = bed(band,"RCR_LEFT")      # banded left   = LCCA-B

labels=["CCA\n(control)","RCCA-B","LCCA-B"]; colors=[MUTED,BLUE,GREEN]
fig,(axR,axC)=plt.subplots(1,2,figsize=(9.0,4.0),constrained_layout=True)
x=np.arange(3)
for ax,vals,title,unit,fmt in [
    (axR,[R_cca,R_rb,R_lb],"terminal-bed resistance $R$","mmHg·s/ml","{:.0f}"),
    (axC,[C_cca,C_rb,C_lb],"terminal-bed compliance $C$","ml/mmHg","{:.2e}")]:
    ax.bar(x,vals,0.62,color=colors)
    for xi,v in zip(x,vals):
        ax.text(xi,v,fmt.format(v),ha="center",va="bottom",fontsize=9,color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=9.5)
    ax.set_title(f"{title}  [{unit}]",fontsize=11,color=INK); ax.margins(y=0.16)
fig.suptitle("Settled calibration knobs: carotid terminal-bed $R$ and $C$",fontsize=12.5,color=INK)
fig.savefig(f"{OUT}/fig_knobs.png"); plt.close(fig)
print("wrote",f"{OUT}/fig_knobs.png")
print(f"  CCA   : R={R_cca:.0f}  C={C_cca:.2e}")
print(f"  RCCA-B: R={R_rb:.0f}  C={C_rb:.2e}")
print(f"  LCCA-B: R={R_lb:.0f}  C={C_lb:.2e}")
