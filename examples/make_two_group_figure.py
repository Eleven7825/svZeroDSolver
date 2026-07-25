"""Two-group calibration result vs Eberth Table 1.
Compares the model (control + banded configs) against the paper for the three
carotid states (CCA control, RCCA-B, LCCA-B) on PI, mean flow, MAP, pulse pressure.
"""
import pysvzerod
import numpy as np
import matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

REPO="/home/shiyi/projects/svZeroDSolver"; OUT=f"{REPO}/examples/figures"
BLUE, MUTED, INK, INK2, GRID = "#2a78d6", "#8a8a86", "#0b0b0b", "#52514e", "#e6e6e2"
mpl.rcParams.update({"figure.dpi":200,"savefig.dpi":200,"savefig.bbox":"tight","font.size":10.5,
    "font.family":"sans-serif","axes.edgecolor":MUTED,"axes.labelcolor":INK,"text.color":INK,
    "xtick.color":INK2,"ytick.color":INK2,"axes.linewidth":0.8,"axes.grid":True,"grid.color":GRID,
    "grid.linewidth":0.8,"axes.axisbelow":True,"legend.frameon":False,
    "axes.spines.top":False,"axes.spines.right":False})

def vstats(cfg, seg):
    df=pysvzerod.simulate(cfg); d=df[df.name==seg]; q=d.flow_in; p=d.pressure_in
    return dict(PI=(q.max()-q.min())/q.mean(), Q=q.mean(), MAP=p.mean(), PP=p.max()-p.min())

ctrl=vstats(f"{REPO}/examples/eberth_control_group.json","rcca_b")
rb  =vstats(f"{REPO}/examples/eberth_banded_group.json","rcca_b")
lb  =vstats(f"{REPO}/examples/eberth_banded_group.json","lcca_b")
model={"CCA":ctrl,"RCCA-B":rb,"LCCA-B":lb}
paper={"CCA":dict(PI=1.16,Q=0.016,MAP=91.7,PP=42.5),
       "RCCA-B":dict(PI=3.11,Q=0.022,MAP=86.5,PP=56.3),
       "LCCA-B":dict(PI=1.65,Q=0.012,MAP=76.8,PP=26.7)}

labels=["CCA\n(control)","RCCA-B","LCCA-B"]; keys=["CCA","RCCA-B","LCCA-B"]
panels=[("PI","pulsatility index",""),("Q","mean flow","ml/s"),
        ("MAP","MAP","mmHg"),("PP","pulse pressure","mmHg")]
fig,axes=plt.subplots(1,4,figsize=(13.5,3.9),constrained_layout=True)
x=np.arange(3); w=0.38
for ax,(q,title,unit) in zip(axes,panels):
    pv=[paper[k][q] for k in keys]; mv=[model[k][q] for k in keys]
    b1=ax.bar(x-w/2,pv,w,color=MUTED,label="Eberth Table 1")
    b2=ax.bar(x+w/2,mv,w,color=BLUE,label="model (two-group)")
    for rects,vals,c in [(b1,pv,INK2),(b2,mv,BLUE)]:
        for r,v in zip(rects,vals):
            ax.text(r.get_x()+r.get_width()/2,v,(f"{v:.3f}" if q=="Q" else f"{v:.0f}" if q in("MAP","PP") else f"{v:.2f}"),
                    ha="center",va="bottom",fontsize=8,color=c)
    ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=9)
    ax.set_title(f"{title}"+(f"  [{unit}]" if unit else ""),fontsize=11,color=INK)
    ax.margins(y=0.15)
h,l=axes[0].get_legend_handles_labels()
fig.legend(h,l,loc="lower center",ncol=2,fontsize=9.5,bbox_to_anchor=(0.5,-0.06))
fig.suptitle("Two-group calibration vs Eberth Table 1 (control CCA + banded RCCA-B/LCCA-B)",fontsize=12.5,color=INK)
fig.savefig(f"{OUT}/fig_two_group.png"); plt.close(fig)
print("wrote",f"{OUT}/fig_two_group.png")
for k in keys: print(f"  {k:8} model PI {model[k]['PI']:.2f} Q {model[k]['Q']:.4f} MAP {model[k]['MAP']:.0f} PP {model[k]['PP']:.0f}")
