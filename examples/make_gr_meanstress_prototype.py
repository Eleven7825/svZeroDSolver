"""G&R prototype, step 1: does a MEAN-PRESSURE-ONLY stimulus explain Eberth's
remodeled carotid geometry?

This is the fast go/no-go test in the staged G&R plan: every one of the four
teaching-repo theories in Growth-Remodeling/src/gr (kinematic growth, full CMM,
homogenized CMM, equilibrated CMM) -- and svGrowth's own default gains -- grows
tissue off a SINGLE stimulus, mean circumferential (intramural) Cauchy stress,
sigma_bar = P*a/h via the Laplace law. We use the fastest of the four (the
equilibrated CMM, an O(1) algebraic solve -- Latorre & Humphrey 2018) to ask: if
that is literally the only signal cells sense, does feeding it our calibrated
per-vessel MAP (relative to baseline CCA) reproduce Table 1's measured chronic
geometry?

We do NOT re-fit the generic teaching-grade constituent parameters (elastin/
collagen/smc stiffnesses) to mouse carotid tissue -- that would take a separate
mechanical-testing calibration this repo doesn't have. So we compare FOLD-CHANGE
(relative to the CCA baseline), which the equilibrated-CMM solution gives for
free and does not depend on those absolute material constants -- only on the
load factor (our MAP ratio) and the model's dimensionless gains/turnover
constants (which are of-order-1 teaching values, not mouse-specific, but the
qualitative test -- does mean pressure alone predict thickening or thinning --
does not hinge on their exact value).

Requires the Growth-Remodeling repo checked out one level up (../../Growth-Remodeling).
"""
import sys, os
sys.path.insert(0, os.path.expanduser("~/projects/Growth-Remodeling/src"))
from gr.parameters import Model, Insult
from gr.geometry import artery
from gr import equilibrated_cmm as eq

import make_control_group as mcg   # for the calibrated CCA MAP (92 mmHg)
import make_banded_group as mbg    # for the calibrated RCCA-B/LCCA-B MAP (86.5/76.8 mmHg)
import make_configs as mc          # for Table-1 "@ MAP" geometry (the ground truth)

MAP_CCA = mcg.MAP                  # 92.0 mmHg
MAP_RCCA, MAP_LCCA = mbg.MAP_A, mbg.MAP_B   # 86.5, 76.8 mmHg

model = Model()   # generic teaching-grade constituents (elastin/collagen/smc); see module docstring
geom = artery(model)


def solve(map_group):
    """Sustained equilibrated-CMM solve, load factor = this group's MAP / CCA MAP."""
    r = eq.solve(geom, Insult(pressure_factor=map_group / MAP_CCA))
    return r.lam, (r.mass / r.lam if r.exists else float("nan")), r.exists


def table1_fold(name):
    d0, h0 = mc.CAROTID["CCA_baseline"]["d"], mc.CAROTID["CCA_baseline"]["h"]
    d, h = mc.CAROTID[name]["d"], mc.CAROTID[name]["h"]
    return d / d0, h / h0


if __name__ == "__main__":
    print("Mean-pressure-only stimulus (equilibrated CMM) vs. Eberth Table 1 (@ MAP)\n")
    print(f"{'group':10}{'MAP factor':>11}{'':>3}{'model: radius fold':>20}{'model: wall fold':>18}"
          f"{'':>3}{'Table1: diam fold':>19}{'Table1: wall fold':>19}")
    rows = [("CCA (ref)", MAP_CCA, "CCA_baseline"),
            ("RCCA-B", MAP_RCCA, "RCCA_B"),
            ("LCCA-B", MAP_LCCA, "LCCA_B")]
    for name, MAP, key in rows:
        lam, thick_fold, exists = solve(MAP)
        d_fold, h_fold = table1_fold(key)
        print(f"{name:10}{MAP/MAP_CCA:>11.3f}{'':>3}{lam:>20.4f}{thick_fold:>18.4f}"
              f"{'':>3}{d_fold:>19.3f}{h_fold:>19.3f}")

    print("""
Reading this: the model column is what a mean-pressure-only stimulus predicts
(radius fold = lambda*, wall fold = mass/lambda* from the equilibrated CMM);
the Table1 columns are what Eberth actually measured. Both banded carotids have
LOWER MAP than baseline CCA, so a mean-pressure-only stimulus predicts mild
THINNING in both (wall fold < 1) and near-unchanged radius. Table 1 shows the
opposite in both cases -- pronounced THICKENING (RCCA-B wall fold 3.55x, LCCA-B
1.68x) and diverging radius (RCCA-B dilates, LCCA-B constricts). Mean pressure
(hence mean intramural stress, the standard stimulus in every gr/svGrowth
theory) cannot be the sole driver of this remodeling -- consistent with
Eberth's own finding (wall thickness correlates with pulse pressure, r*=0.632,
but NOT mean pressure, r*=0.020 NS). This motivates adding a
flow/WSS/pulsatility channel (svGrowth's existing wss gain, or a new PI-based
one) as the next step.""")
