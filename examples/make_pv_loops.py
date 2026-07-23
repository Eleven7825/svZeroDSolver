"""Plot LV pressure-volume loops from the experimental elastance-heart model,
pre-op vs acute band, showing the afterload response (loop shifts up, SV falls).

Uses the committed eberth_elastance_{prebanding,acute}.json configs.
"""
import json, pysvzerod
import numpy as np
import matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

REPO = "/home/shiyi/projects/svZeroDSolver"
OUT = f"{REPO}/examples/figures"

BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#e6e6e2"
mpl.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 11, "font.family": "sans-serif", "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.axisbelow": True, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
})

def loop(mode):
    cfg = json.load(open(f"{REPO}/examples/eberth_elastance_{mode}.json"))
    cfg["simulation_parameters"]["output_variable_based"] = True
    df = pysvzerod.simulate(cfg)
    V = df[df.name == "Vc:ventricle"].y.to_numpy() * 1000.0        # ml -> uL
    P = df[df.name == "pressure:ventricle:aortic_valve"].y.to_numpy()
    return V, P

fig, ax = plt.subplots(figsize=(5.8, 4.8))
for mode, color in [("prebanding", BLUE), ("acute", ORANGE)]:
    V, P = loop(mode)
    EDV, ESV = V.max(), V.min(); SV = EDV - ESV; EF = 100 * SV / EDV; ESP = P.max()
    Vc, Pc = np.append(V, V[0]), np.append(P, P[0])   # close the loop
    ax.plot(Vc, Pc, color=color, lw=2.0,
            label=f"{'pre-op' if mode=='prebanding' else 'acute band'}: "
                  f"SV {SV:.0f} µL, EF {EF:.0f}%, ESP {ESP:.0f} mmHg")

ax.set_xlabel("LV volume  [µL]")
ax.set_ylabel("LV pressure  [mmHg]")
ax.set_title("LV pressure–volume loops (elastance-heart model)\n"
             "banding shifts the loop up-and-right: higher pressure, smaller stroke volume",
             fontsize=11, color=INK)
ax.legend(loc="lower center", fontsize=9)
fig.savefig(f"{OUT}/fig_pv_loops.png"); plt.close(fig)
print(f"wrote {OUT}/fig_pv_loops.png")
