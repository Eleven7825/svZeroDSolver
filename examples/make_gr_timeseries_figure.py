"""Flow, pressure, and wall shear stress time series at RCCA-B/LCCA-B (banded
group) and CCA (control group) over each config's own converged cardiac cycle
-- the inputs a future G&R model needs. WSS is not a solver output: computed
via Poiseuille tau_w = 4*mu*Q(t)/(pi*r^3) using each vessel's fixed lumen
radius (Eberth Table 1: chronic RCCA-B/LCCA-B for banded, baseline CCA for
control).
"""
import numpy as np
import pysvzerod
import matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

REPO="/home/shiyi/projects/svZeroDSolver"; OUT=f"{REPO}/examples/figures"
BLUE, GREEN, GREY, INK, INK2, GRID = "#2a78d6", "#008300", "#8a8a86", "#0b0b0b", "#52514e", "#e6e6e2"
mpl.rcParams.update({"figure.dpi":200,"savefig.dpi":200,"savefig.bbox":"tight","font.size":10,
    "font.family":"sans-serif","axes.edgecolor":INK2,"axes.labelcolor":INK,"text.color":INK,
    "xtick.color":INK2,"ytick.color":INK2,"axes.linewidth":0.8,"axes.grid":True,"grid.color":GRID,
    "grid.linewidth":0.8,"axes.axisbelow":True,"legend.frameon":False,
    "axes.spines.top":False,"axes.spines.right":False})

MU = 0.035                              # poise (3.5 cP, Windberger 2003)
R_RCCA = 633.0e-4 / 2                    # cm, chronic RCCA-B lumen radius (Eberth Table 1)
R_LCCA = 482.0e-4 / 2                    # cm, chronic LCCA-B lumen radius
R_CCA  = 496.0e-4 / 2                    # cm, baseline CCA lumen radius (control)

def load(cfg, vessel, r_lumen):
    df = pysvzerod.simulate(f"{REPO}/examples/{cfg}")
    d = df[df.name == vessel].reset_index(drop=True)
    t = d.time.values - d.time.values[0]
    wss = 4 * MU * d.flow_in.values / (np.pi * r_lumen**3)   # dyn/cm^2
    return t, d.flow_in.values, d.pressure_in.values, wss

t_r, q_r, p_r, w_r = load("eberth_banded_group.json", "rcca_b", R_RCCA)
t_l, q_l, p_l, w_l = load("eberth_banded_group.json", "lcca_b", R_LCCA)
t_c, q_c, p_c, w_c = load("eberth_control_group.json", "rcca_b", R_CCA)

series = [("RCCA-B", BLUE, t_r, q_r, p_r, w_r),
          ("LCCA-B", GREEN, t_l, q_l, p_l, w_l),
          ("CCA (control)", GREY, t_c, q_c, p_c, w_c)]
panels = [(r"flow $Q(t)$", "ml/s", 1), (r"pressure $P(t)$", "mmHg", 2), (r"wall shear stress $\tau_w(t)$", r"dyn/cm$^2$", 3)]
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.7), constrained_layout=True)
for ax, (title, unit, idx) in zip(axes, panels):
    for name, col, t, q, p, w in series:
        y = (q, p, w)[idx-1]
        ax.plot(t, y, color=col, lw=2, label=name)
    ax.set_title(f"{title}  [{unit}]", fontsize=10.5, color=INK)
    ax.set_xlabel("time [s]", fontsize=9)
    ax.margins(x=0.02)
h, lab = axes[0].get_legend_handles_labels()
fig.legend(h, lab, loc="lower center", ncol=3, fontsize=9.5, bbox_to_anchor=(0.5, -0.08))
fig.suptitle("RCCA-B / LCCA-B (banded) vs. CCA (control) time series, each over its own converged cycle",
             fontsize=12.5, color=INK)
fig.savefig(f"{OUT}/fig_gr_timeseries.png")
plt.close(fig)
print("wrote", f"{OUT}/fig_gr_timeseries.png")
for name, col, t, q, p, w in series:
    print(f"  {name:15} Q mean {q.mean():.4f} min/max {q.min():.4f}/{q.max():.4f} ml/s"
          f" | P mean {p.mean():.0f} min/max {p.min():.0f}/{p.max():.0f} mmHg"
          f" | WSS mean {w.mean():.1f} min/max {w.min():.1f}/{w.max():.1f} dyn/cm^2")
