"""Settled calibration knobs (bed R and bed C, excluding the band S and the
fixed/literature aortic C_a,C_b): control group in grey, banded group in blue
(two columns for the per-carotid knobs: RCCA-B and LCCA-B). Same bar style as
the two-group figure. The two resistance quantities (carotid bed and systemic
bed R_p+R_d) share one panel since both are pinned directly from experimental
(Qbar, MAP) via Ohm's law -- not searched/fit.
"""
import json
import numpy as np
import matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

REPO="/home/shiyi/projects/svZeroDSolver"; OUT=f"{REPO}/examples/figures"
GREY, BLUE, INK, INK2, GRID = "#8a8a86", "#2a78d6", "#0b0b0b", "#52514e", "#e6e6e2"
mpl.rcParams.update({"figure.dpi":200,"savefig.dpi":200,"savefig.bbox":"tight","font.size":10,
    "font.family":"sans-serif","axes.edgecolor":GREY,"axes.labelcolor":INK,"text.color":INK,
    "xtick.color":INK2,"ytick.color":INK2,"axes.linewidth":0.8,"axes.grid":True,"grid.color":GRID,
    "grid.linewidth":0.8,"axes.axisbelow":True,"axes.spines.top":False,"axes.spines.right":False})

def rcr(cfg, name):
    v=[b["bc_values"] for b in cfg["boundary_conditions"] if b["bc_name"]==name][0]
    return v["Rp"]+v["Rd"], v["C"]
ctrl=json.load(open(f"{REPO}/examples/eberth_control_group.json"))
band=json.load(open(f"{REPO}/examples/eberth_banded_group.json"))
Rc,Cc=rcr(ctrl,"RCR_RIGHT")                 # CCA (control, symmetric)
Rr,Cr=rcr(band,"RCR_RIGHT"); Rl,Cl=rcr(band,"RCR_LEFT")   # RCCA-B, LCCA-B
Rsys_c=rcr(ctrl,"RCR_SYS")[0]; Rsys_b=rcr(band,"RCR_SYS")[0]

fig,axes=plt.subplots(1,2,figsize=(8.0,3.9),constrained_layout=True)

# panel 1: all R_p+R_d resistances together (carotid bed + systemic bed) --
# a visual gap separates the two clusters; both are pinned directly from
# experimental (Qbar, MAP), not fit.
ax = axes[0]
x = [0,1,2,3.5,4.5]
vals = [Rc,Rr,Rl,Rsys_c,Rsys_b]
cols = [GREY,BLUE,BLUE,GREY,BLUE]
xl = ["CCA","RCCA-B","LCCA-B","sys\n(control)","sys\n(banded)"]
ax.bar(x, vals, 0.62, color=cols)
for xi,v in zip(x,vals): ax.text(xi,v,f"{v:.0f}",ha="center",va="bottom",fontsize=8.5,color=INK2)
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=8.5)
ax.set_title(r"bed $R_p\!+\!R_d$ (carotid + systemic)  [mmHg·s/ml]",fontsize=10,color=INK)
ax.margins(y=0.17)

# panel 2: bed C (the only real knob among these two panels)
ax = axes[1]
x = np.arange(3)
vals = [Cc,Cr,Cl]
ax.bar(x, vals, 0.62, color=[GREY,BLUE,BLUE])
for xi,v in zip(x,vals): ax.text(xi,v,f"{v:.1e}",ha="center",va="bottom",fontsize=8.5,color=INK2)
ax.set_xticks(x); ax.set_xticklabels(["CCA","RCCA-B","LCCA-B"],fontsize=9)
ax.set_title(r"bed $C$  [ml/mmHg]",fontsize=10.5,color=INK); ax.margins(y=0.17)

# legend proxies
import matplotlib.patches as mp
h=[mp.Patch(color=GREY,label="control group"),mp.Patch(color=BLUE,label="banded group")]
fig.legend(handles=h,loc="lower center",ncol=2,fontsize=10,bbox_to_anchor=(0.5,-0.14))
fig.suptitle(r"Settled calibration knobs: control vs banded"+"\n"+
             r"(left: $R_p\!+\!R_d$, derived directly from experimental $\bar Q$,MAP -- not fit; right: bed $C$, the actual knob)",
             fontsize=11,color=INK)
fig.savefig(f"{OUT}/fig_knobs.png"); plt.close(fig)
print("wrote",f"{OUT}/fig_knobs.png")
print(f"  bed R : CCA {Rc:.0f} | RCCA-B {Rr:.0f} | LCCA-B {Rl:.0f}")
print(f"  bed C : CCA {Cc:.2e} | RCCA-B {Cr:.2e} | LCCA-B {Cl:.2e}")
print(f"  sys R : control {Rsys_c:.0f} | banded {Rsys_b:.0f}")
