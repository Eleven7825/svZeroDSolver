"""G&R prototype, step 1: does a MEAN-PRESSURE-ONLY stimulus, with a REAL mouse
carotid constitutive law, explain Eberth's remodeled carotid geometry?

Every one of the four teaching-repo theories in Growth-Remodeling/src/gr
(kinematic growth, full CMM, homogenized CMM, equilibrated CMM) -- and svGrowth's
own default gains -- grows tissue off a SINGLE stimulus, mean circumferential
(intramural) Cauchy stress, sigma_bar = P*a/h via the Laplace law. We use the
fastest of the four (the equilibrated CMM, an O(1) algebraic solve -- Latorre &
Humphrey 2018) to ask: if that is literally the only signal cells sense, does
feeding it our calibrated per-vessel MAP (relative to baseline CCA) reproduce
Table 1's measured chronic geometry?

-------------------------------------------------------------------------------
Where the numbers come from (replacing gr's generic "large elastic artery"
teaching defaults, which its own docstring says are "not a fit to any one
animal")
-------------------------------------------------------------------------------
* Constitutive law (elastin NeoHookean c; circumferential collagen+SMC Fung
  c1,c2) -- Bersi, Ferruzzi, Eberth JF, Gleason, Humphrey (2014), "Consistent
  Biomechanical Phenotyping of Common Carotid Arteries from Seven Genetic,
  Pharmacological, and Surgical Mouse Models," Ann Biomed Eng 42(6):1207-1223,
  Table 2, the "Aortic Banding, control" fit -- the SAME aortic-banding mouse
  model as Eberth (2009) (J.F. Eberth is a co-author on both), four-fiber-family
  model, pooled control group (n=16):
      elastin c        = 8.126 kPa
      circumferential collagen+SMC  c12 = 4.782 kPa,  c22 = 0.041
  We use the CONTROL fit as the tissue's fixed material law for ALL THREE rows
  below (CCA / RCCA-B / LCCA-B) -- i.e. we do NOT swap in that paper's
  "35-56 day post-banding" fit (c12=9.920, c22=11.579) as a different material
  law for the banded vessels. Doing so would hand the model the answer: a
  constrained-mixture theory's whole premise is that each constituent's
  intrinsic stress-stretch behavior is FIXED, and only mass/composition/
  natural-configuration evolve under a stimulus. If we want to test "does mean
  pressure alone, acting on a fixed material law, predict the remodeling,"
  the post-band fit cannot be an input -- it can only be a (very informative)
  independent check, see the printed comparison at the bottom.
* Reference geometry (R = mid-wall radius, P_h = MAP) -- our own calibrated CCA
  baseline (make_configs.CAROTID["CCA_baseline"], make_control_group.MAP).
* Mass fractions, deposition stretches, SMC fiber parameters, turnover
  rate/gain -- NEITHER Bersi 2014 NOR Eberth 2009 report these (checked:
  no mass-fraction or pre-stretch table in either paper; SMC is lumped into
  the "circumferential collagen+SMC" family, not separately identified). These
  stay at gr's generic teaching-grade defaults -- flagged here as an
  assumption, not a literature value.

Requires the Growth-Remodeling repo checked out one level up (../../Growth-Remodeling).
"""
import sys, os
sys.path.insert(0, os.path.expanduser("~/projects/Growth-Remodeling/src"))
from gr.parameters import Model, Insult, Constituent, default_constituents
from gr.mechanics import NeoHookean, FungFiber
from gr.geometry import artery
from gr import equilibrated_cmm as eq

import make_control_group as mcg   # for the calibrated CCA MAP (92 mmHg)
import make_banded_group as mbg    # for the calibrated RCCA-B/LCCA-B MAP (86.5/76.8 mmHg)
import make_configs as mc          # for Table-1 "@ MAP" geometry (the ground truth)

MMHG_TO_KPA = 0.133322
MAP_CCA = mcg.MAP                           # 92.0 mmHg
MAP_RCCA, MAP_LCCA = mbg.MAP_A, mbg.MAP_B   # 86.5, 76.8 mmHg

# ---- mouse-carotid reference geometry: CCA baseline mid-wall radius @ MAP ----
d0, h0 = mc.CAROTID["CCA_baseline"]["d"], mc.CAROTID["CCA_baseline"]["h"]   # um
R_MID_MM = (d0 / 2.0 + h0 / 2.0) * 1e-3     # lumen radius + half wall -> mid-wall, mm
P_H_KPA = MAP_CCA * MMHG_TO_KPA

# ---- Bersi/Ferruzzi/Eberth/Gleason/Humphrey (2014) Table 2, control fit -----
BERSI_CONTROL = dict(elastin_c=8.126, collagen_c1=4.782, collagen_c2=0.041)
BERSI_BANDING_35_56D = dict(elastin_c=4.025, collagen_c1=9.920, collagen_c2=11.579)   # reported only, NOT used as model input -- see docstring

_generic = {c.name: c for c in default_constituents()}   # for G, k_d, gain, smc law (unfit, flagged)
model = Model(
    R=R_MID_MM,
    P_h=P_H_KPA,
    constituents=[
        Constituent("elastin", phi0=_generic["elastin"].phi0, G=_generic["elastin"].G,
                    law=NeoHookean(c=BERSI_CONTROL["elastin_c"]),
                    k_d=_generic["elastin"].k_d, gain=_generic["elastin"].gain, degradable=True),
        Constituent("collagen", phi0=_generic["collagen"].phi0, G=_generic["collagen"].G,
                    law=FungFiber(c1=BERSI_CONTROL["collagen_c1"], c2=BERSI_CONTROL["collagen_c2"]),
                    k_d=_generic["collagen"].k_d, gain=_generic["collagen"].gain),
        Constituent("smc", phi0=_generic["smc"].phi0, G=_generic["smc"].G,
                    law=_generic["smc"].law,   # unchanged generic Fung -- not separately reported by either paper
                    k_d=_generic["smc"].k_d, gain=_generic["smc"].gain),
    ],
)
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
    print("Mean-pressure-only stimulus (equilibrated CMM, real mouse-carotid constitutive law)")
    print("vs. Eberth Table 1 (@ MAP)\n")
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

    print(f"""
Reading this: the model column is what a mean-pressure-only stimulus, acting on
a FIXED real mouse-carotid material law (Bersi/Ferruzzi/Eberth/Gleason/Humphrey
2014 control fit), predicts (radius fold = lambda*, wall fold = mass/lambda*
from the equilibrated CMM); the Table1 columns are what Eberth actually
measured. Both banded carotids have LOWER MAP than baseline CCA, so a
mean-pressure-only stimulus predicts mild THINNING in both (wall fold < 1) and
near-unchanged radius. Table 1 shows the opposite in both cases -- pronounced
THICKENING (RCCA-B wall fold 3.55x, LCCA-B 1.68x) and diverging radius
(RCCA-B dilates, LCCA-B constricts).

A second, independent line of evidence points the same way. The equilibrated
CMM (like every theory in this package) assumes each constituent's intrinsic
stress-stretch behavior is FIXED over time -- only mass and natural
configuration evolve. But the SAME source paper's own "35-56 day post-banding"
fit is a genuinely different material:
    circumferential collagen+SMC (c12, c22):  control {BERSI_CONTROL['collagen_c1']:.3f}, {BERSI_CONTROL['collagen_c2']:.3f}  ->  banded {BERSI_BANDING_35_56D['collagen_c1']:.3f}, {BERSI_BANDING_35_56D['collagen_c2']:.3f}
    elastin c (kPa):                          control {BERSI_CONTROL['elastin_c']:.3f}  ->  banded {BERSI_BANDING_35_56D['elastin_c']:.3f}
A ~280x jump in the stiffening exponent c22 is not a fixed-law CMM predicting
"more of the same material" -- it is new collagen being laid down differently
(Eberth 2009 itself reports the new RCCA-B collagen fibers as thinner and less
ordered; fiber angle shifts 29.4 deg -> 37.1 deg). No fixed-material-law,
mean-stress-only theory in this package can produce that by construction.

Together: mean pressure (hence mean intramural stress) cannot be the sole
driver of this remodeling -- consistent with Eberth's own finding (wall
thickness correlates with pulse pressure, r*=0.632, but NOT mean pressure,
r*=0.020 NS). This motivates adding a flow/WSS/pulsatility channel (svGrowth's
existing wss gain, or a new PI-based one) as the next step.""")
