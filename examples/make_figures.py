"""Generate illustrative figures for the Eberth aortic-banding 0D model.
Applies the dataviz method: colorblind-safe categorical palette (validated
blue/green/... slots), thin 2px marks, recessive grid, legend + direct labels,
single y-axis per panel."""
import json, copy, math
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pysvzerod

REPO = "/home/shiyi/projects/svZeroDSolver"
# Calibrated, physiologically-parameterized model (stenosis_coefficient = 53.8).
CFG = f"{REPO}/examples/eberth_aortic_band_physio.json"
OUT = f"{REPO}/examples/figures"
import os; os.makedirs(OUT, exist_ok=True)

# ---- validated categorical palette (dataviz references/palette.md, light) ----
BLUE, GREEN, MAGENTA, YELLOW = "#2a78d6", "#008300", "#e87ba4", "#eda100"
ORANGE, RED = "#eb6834", "#e34948"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a86"
GRID = "#e6e6e2"

mpl.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 11, "font.family": "sans-serif", "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.linewidth": 0.8, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False,
})
LW = 2.0

with open(CFG) as f:
    base = json.load(f)

def last_cycle(df, seg):
    d = df[df.name == seg].copy()
    # solver writes only the last cycle here; keep as-is, sort by time
    return d.sort_values("time")

def run(cfg):
    return pysvzerod.simulate(cfg)

df = run(base)
segs = {"ascending_aorta": "Ascending aorta", "rcca_b": "RCCA-B (node A, upstream)",
        "aortic_band": "Band", "lcca_b": "LCCA-B (node B, downstream)",
        "descending_aorta": "Descending aorta"}

# =====================================================================
# FIG 1 — Pressure waveforms: node A (RCCA-B) vs node B (LCCA-B)
# =====================================================================
fig, ax = plt.subplots(figsize=(6.2, 4.0))
a = last_cycle(df, "rcca_b"); b = last_cycle(df, "lcca_b")
t0 = a.time.min()
ax.plot(a.time - t0, a.pressure_in, color=BLUE, lw=LW, label="RCCA-B (upstream)")
ax.plot(b.time - t0, b.pressure_in, color=GREEN, lw=LW, label="LCCA-B (downstream)")
# annotate the pulse-pressure gap at systolic peak instead of colliding end-labels
ipk = a.pressure_in.values.argmax(); xpk = (a.time - t0).values[ipk]
ax.annotate("", xy=(xpk, a.pressure_in.max()), xytext=(xpk, b.pressure_in.max()),
            arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.0))
ax.text(xpk + 0.004, (a.pressure_in.max() + b.pressure_in.max())/2,
        f"pulse\n{a.pressure_in.max()-a.pressure_in.min():.0f} vs {b.pressure_in.max()-b.pressure_in.min():.0f} mmHg",
        color=INK2, fontsize=9, va="center")
ax.set_xlabel("time within cardiac cycle  [s]")
ax.set_ylabel("pressure  [mmHg]")
ax.set_title("Pulse is damped downstream of the band\n(same mean, smaller swing at LCCA-B)", fontsize=11.5, color=INK)
ax.legend(loc="upper right")
fig.savefig(f"{OUT}/fig1_pressure_waveforms.png"); plt.close(fig)

# =====================================================================
# FIG 2 — Flow waveforms: right vs left carotid (+ inflow reference)
# =====================================================================
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.axhline(0, color=MUTED, lw=0.8)
ax.plot(a.time - t0, a.flow_in, color=BLUE, lw=LW, label="RCCA-B (right carotid)")
ax.plot(b.time - t0, b.flow_in, color=GREEN, lw=LW, label="LCCA-B (left carotid)")
ax.set_xlabel("time within cardiac cycle  [s]")
ax.set_ylabel("flow  [ml/s]")
ax.set_title("Carotid inflow waveforms", fontsize=11.5, color=INK)
ax.legend(loc="upper right")
fig.savefig(f"{OUT}/fig2_carotid_flow_waveforms.png"); plt.close(fig)

# =====================================================================
# FIG 3 — Band characteristic: ΔP vs Q over the cycle (nonlinear loop)
# =====================================================================
fig, ax = plt.subplots(figsize=(5.6, 4.4))
bd = last_cycle(df, "aortic_band")
dP = bd.pressure_in.values - bd.pressure_out.values
Q = bd.flow_in.values
ax.plot(Q, dP, color=BLUE, lw=LW)
# quadratic reference (S|Q|Q) using the config's stenosis_coefficient
S = base["vessels"][3]["zero_d_element_values"]["stenosis_coefficient"]
R = base["vessels"][3]["zero_d_element_values"]["R_poiseuille"]
qq = np.linspace(Q.min(), Q.max(), 100)
ax.plot(qq, (R + S*np.abs(qq))*qq, color=MUTED, lw=1.4, ls="--", label=r"$(R+S|Q|)\,Q$  (static)")
ax.text(Q[np.argmax(dP)], dP.max(), "  cycle loop\n  (inertance L)", color=BLUE, fontsize=9.5, va="top")
ax.set_xlabel("flow through band  Q  [ml/s]")
ax.set_ylabel(r"$\Delta P_\mathrm{band}$  [mmHg]")
ax.set_title("Nonlinear band characteristic (Young–Tsai)", fontsize=11.5, color=INK)
ax.legend(loc="upper left")
fig.savefig(f"{OUT}/fig3_band_deltaP_Q.png"); plt.close(fig)

# =====================================================================
# FIG 4 — Pulsatility index by location: model vs Eberth Table 1
# =====================================================================
def PI(seg):
    d = last_cycle(df, seg); q = d.flow_in
    return (q.max() - q.min()) / q.mean()
model = {"RCCA-B": PI("rcca_b"), "LCCA-B": PI("lcca_b")}
target = {"CCA (baseline)": 1.16, "RCCA-B": 3.11, "LCCA-B": 1.65}
labels = ["CCA (baseline)", "RCCA-B", "LCCA-B"]
tvals = [target[k] for k in labels]
mvals = [np.nan, model["RCCA-B"], model["LCCA-B"]]
x = np.arange(len(labels)); w = 0.38
fig, ax = plt.subplots(figsize=(6.2, 4.0))
b1 = ax.bar(x - w/2, tvals, w, color=MUTED, label="Eberth Table 1 (target)")
b2 = ax.bar(x + w/2, mvals, w, color=BLUE, label="this model")
for rect, v in zip(b1, tvals):
    ax.text(rect.get_x()+rect.get_width()/2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9, color=INK2)
for rect, v in zip(b2, mvals):
    if not math.isnan(v):
        ax.text(rect.get_x()+rect.get_width()/2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9, color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("pulsatility index  PI = (Q$_{max}$−Q$_{min}$)/Q$_{mean}$")
ax.set_title("Pulsatility by location — calibrated model vs paper Table 1\n(stenosis calibrated to the PI ratio; both PIs within ~7%)", fontsize=11.5, color=INK)
ax.legend(loc="upper left")
fig.savefig(f"{OUT}/fig4_PI_validation.png"); plt.close(fig)

# =====================================================================
# FIG 5 — Pressure envelope along the aortic path
# =====================================================================
path = ["ascending_aorta", "aortic_band", "descending_aorta"]
xpos = [0, 1, 2]; xlab = ["ascending\naorta", "band\n(node A→B)", "descending\naorta"]
sysP, diaP, meanP = [], [], []
for seg in path:
    d = last_cycle(df, seg)
    sysP.append(d.pressure_in.max()); diaP.append(d.pressure_in.min()); meanP.append(d.pressure_in.mean())
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.fill_between(xpos, diaP, sysP, color=BLUE, alpha=0.16, label="systolic–diastolic range")
ax.plot(xpos, meanP, color=BLUE, lw=LW, marker="o", ms=7, label="mean (MAP)")
ax.plot(xpos, sysP, color=INK2, lw=1.2, ls=":"); ax.plot(xpos, diaP, color=INK2, lw=1.2, ls=":")
for xp, s, d_ in zip(xpos, sysP, diaP):
    ax.annotate("", xy=(xp, s), xytext=(xp, d_), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.0))
    ax.text(xp+0.04, (s+d_)/2, f"{s-d_:.0f}", color=INK2, fontsize=9, va="center")
ax.set_xticks(xpos); ax.set_xticklabels(xlab)
ax.set_ylabel("pressure  [mmHg]")
ax.set_title("Pressure along the aortic path\n(pulse pressure [mmHg] collapses across the band; MAP ~flat)", fontsize=11.5, color=INK)
ax.legend(loc="lower left")
fig.savefig(f"{OUT}/fig5_pressure_envelope.png"); plt.close(fig)

# =====================================================================
# FIG 6 — Sensitivity sweep: PI vs band severity (stenosis_coefficient)
# =====================================================================
Svals = [0, 20, 40, 80, 160, 320, 640]
piR, piL = [], []
for s in Svals:
    cfg = copy.deepcopy(base)
    cfg["vessels"][3]["zero_d_element_values"]["stenosis_coefficient"] = float(s)
    d = run(cfg)
    def pi(seg):
        dd = d[d.name == seg]; q = dd.flow_in
        return (q.max()-q.min())/q.mean()
    piR.append(pi("rcca_b")); piL.append(pi("lcca_b"))
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.plot(Svals, piR, color=BLUE, lw=LW, marker="o", ms=6, label="RCCA-B (upstream)")
ax.plot(Svals, piL, color=GREEN, lw=LW, marker="s", ms=6, label="LCCA-B (downstream)")
S_cal = base["vessels"][3]["zero_d_element_values"].get("stenosis_coefficient", 0)
ax.axvline(S_cal, color=MUTED, lw=1.0, ls="--")
ax.text(S_cal, ax.get_ylim()[1], f" calibrated S={S_cal:g}", color=MUTED, fontsize=9, va="top")
ax.set_xlabel("band severity — stenosis_coefficient  S")
ax.set_ylabel("pulsatility index  PI")
ax.set_title("Sensitivity: band severity sets the pulsatility split", fontsize=11.5, color=INK)
ax.legend(loc="center right")
fig.savefig(f"{OUT}/fig6_severity_sweep.png"); plt.close(fig)

print("PI model:", {k: round(v,2) for k,v in model.items()})
print("pulse pressures [sys,dia,mean] along path:")
for seg,s,d_,m in zip(path,sysP,diaP,meanP):
    print(f"  {seg:18s} sys={s:.1f} dia={d_:.1f} mean={m:.1f} pp={s-d_:.1f}")
print("sweep S:", Svals)
print("  PI_R:", [round(v,2) for v in piR])
print("  PI_L:", [round(v,2) for v in piL])
print("wrote figures to", OUT)
