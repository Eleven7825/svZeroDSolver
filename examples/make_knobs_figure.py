"""Settled calibration knobs (bed R and bed C, excluding the band S and the
fixed/literature aortic C_a,C_b): control group in grey, banded group in blue
(two columns for the per-carotid knobs: RCCA-B and LCCA-B). Same bar style as
the two-group figure.
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

# each panel: (title, unit, values, colors, xlabels, fmt)
# titles use the same symbols as the circuit diagram (slide 4): R_p+R_d and C
# for the RCR_RIGHT/RCR_LEFT/RCR_SYS beds.
panels=[
 (r"bed $R_p\!+\!R_d$","mmHg·s/ml",[Rc,Rr,Rl],[GREY,BLUE,BLUE],["CCA","RCCA-B","LCCA-B"],"{:.0f}"),
 (r"bed $C$","ml/mmHg",[Cc,Cr,Cl],[GREY,BLUE,BLUE],["CCA","RCCA-B","LCCA-B"],"{:.1e}"),
 (r"systemic $R_p^S\!+\!R_d^S$","mmHg·s/ml",[Rsys_c,Rsys_b],[GREY,BLUE],["control","banded"],"{:.0f}"),
]
fig,axes=plt.subplots(1,3,figsize=(10.5,3.9),constrained_layout=True)
for ax,(title,unit,vals,cols,xl,fmt) in zip(axes,panels):
    x=np.arange(len(vals))
    ax.bar(x,vals,0.62,color=cols)
    for xi,v in zip(x,vals): ax.text(xi,v,fmt.format(v),ha="center",va="bottom",fontsize=8.5,color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=9)
    ax.set_title(f"{title}  [{unit}]",fontsize=10.5,color=INK); ax.margins(y=0.17)
# legend proxies
import matplotlib.patches as mp
h=[mp.Patch(color=GREY,label="control group"),mp.Patch(color=BLUE,label="banded group")]
fig.legend(handles=h,loc="lower center",ncol=2,fontsize=10,bbox_to_anchor=(0.5,-0.14))
fig.suptitle(r"Settled calibration knobs ($R_p\!+\!R_d$ and bed $C$; band $S$ and fixed aortic $C_a,C_b$ excluded): control vs banded",fontsize=12,color=INK)
fig.savefig(f"{OUT}/fig_knobs.png"); plt.close(fig)
print("wrote",f"{OUT}/fig_knobs.png")
print(f"  bed R : CCA {Rc:.0f} | RCCA-B {Rr:.0f} | LCCA-B {Rl:.0f}")
print(f"  bed C : CCA {Cc:.2e} | RCCA-B {Cr:.2e} | LCCA-B {Cl:.2e}")
print(f"  sys R : control {Rsys_c:.0f} | banded {Rsys_b:.0f}")
