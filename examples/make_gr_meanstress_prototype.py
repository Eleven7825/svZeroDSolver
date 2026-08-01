"""G&R prototype, step 1: does a MEAN-PRESSURE-ONLY stimulus, with the closest
available REAL mouse-carotid parameters, explain Eberth's remodeled carotid
geometry?

Every one of the four teaching-repo theories in Growth-Remodeling/src/gr
(kinematic growth, full CMM, homogenized CMM, equilibrated CMM) -- and svGrowth's
own default gains -- grows tissue off a SINGLE stimulus, mean circumferential
(intramural) Cauchy stress, sigma_bar = P*a/h via the Laplace law. We use the
fastest of the four (the equilibrated CMM, an O(1) algebraic solve -- Latorre &
Humphrey 2018) to ask: if that is literally the only signal cells sense, does
feeding it our calibrated per-vessel MAP (relative to baseline CCA) reproduce
Table 1's measured chronic geometry?

-------------------------------------------------------------------------------
Where every parameter comes from (full table + citations also in gr_extension.md)
-------------------------------------------------------------------------------
gr's own generic "large elastic artery" teaching defaults are, by its own
docstring, "not a fit to any one animal." We replace every field we can with a
real mouse common carotid artery (CCA) source, and flag the (few) fields no
paper reports.

* Mass fractions (phi0) and the turnover fiber's deposition stretch (G) --
  Bellini, Ferruzzi, Roccabianca, Di Martino, Humphrey (2014), "A
  Microstructurally Motivated Model of Arterial Wall Mechanics with
  Mechanobiological Implications," Ann Biomed Eng 42(3):488-502. Table 1/text:
  10-week C57BL/6 wild-type mouse CCA, homeostatic MAP ~93 mmHg -- elastin
  phi_e=0.249, media (smooth muscle + circumferential collagen) phi_m=0.279,
  total collagen phi_c=0.458. Deposition stretches (Table/text): circumferential
  collagen and smooth muscle SHARE one value in their model, G_hc2=G_hm in
  [1.07,1.09] -- i.e. this literature does not treat circumferential collagen
  and SMC as mechanically separable, which is exactly why we don't try to
  either (see below).
* Elastin and circumferential-family (collagen+SMC) STIFFNESS -- Bersi,
  Ferruzzi, Eberth JF, Gleason, Humphrey (2014), "Consistent Biomechanical
  Phenotyping of Common Carotid Arteries from Seven Genetic, Pharmacological,
  and Surgical Mouse Models," Ann Biomed Eng 42(6):1207-1223, Table 2,
  "Aortic Banding, control" fit (n=16) -- the SAME aortic-banding mouse model as
  Eberth 2009 (J.F. Eberth is a co-author on both), same isotropic-neo-Hookean
  + Fung-exponential functional forms gr uses:
      elastin c = 8.126 kPa
      circumferential collagen+SMC family: c1(=c12) = 4.782 kPa, c2(=c22) = 0.041
  We deliberately use ONLY this "control" fit as the fixed material law for ALL
  THREE rows (CCA / RCCA-B / LCCA-B), never the same paper's "35-56 day
  post-banding" fit (c1=9.920, c2=11.579). A constrained-mixture theory's
  premise is that each constituent's intrinsic stress-stretch behavior is
  FIXED; only mass/composition/natural-configuration evolve under a stimulus.
  Swapping in the post-band fit as a different material law for the banded
  rows would hand the model the answer instead of testing whether mean
  pressure predicts it. The post-band fit is used only as an independent check
  (printed at the bottom), never as a solve input.
* Why only 2 constituents, not gr's usual 3 (elastin/collagen/smc): BOTH
  mouse-CCA papers above model circumferential collagen and smooth muscle as
  ONE mechanical family (not separable) -- Bellini's own notation sets their
  deposition stretches literally equal, and Bersi's four-fiber-family fit
  reports one (c1,c2) pair for "circumferential collagen+SMC" together, not
  two. Forcing gr's separate elastin/collagen/smc structure here would mean
  inventing a collagen-vs-SMC split with no cited source and double-counting
  the real family's stiffness across two constituents. So we build a 2-
  constituent Model (elastin; one turnover fiber representing circumferential
  collagen+SMC combined) -- gr's math (parameters.py, equilibrated_cmm.py)
  supports any number of turnover constituents with no code change.
* Reference geometry (R = mid-wall radius, P_h = MAP) -- our OWN calibrated CCA
  baseline (make_configs.CAROTID["CCA_baseline"], make_control_group.MAP), not
  Bellini's own reported mouse cohort (r_i=299 um, h=26 um, MAP~93 mmHg) --
  those mice are a different cohort/strain/age than Eberth's, and the whole
  point of this project is Eberth's specific animals.
* NOT available from either paper, kept at gr's generic default (flagged): the
  turnover fiber's degradation rate k_d and mechanosensitivity gain. These
  don't actually enter the equilibrated-CMM's algebraic equilibrium equation at
  all (turnover rate only matters for TRANSIENT dynamics -- homogenized/full
  CMM, not this fast equilibrium test) -- so this gap doesn't affect the result
  below, only future transient-model steps. Elastin's deposition stretch G is
  also kept at gr's generic default: Bellini's elastin deposition stretch is an
  anisotropic (circumferential/axial/radial) triad for a structurally
  different, anisotropic elastin law, not compatible with gr's simpler
  isotropic single-scalar elastin formulation.

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

# ---- mouse-carotid reference geometry: OUR CCA baseline mid-wall radius @ MAP
d0, h0 = mc.CAROTID["CCA_baseline"]["d"], mc.CAROTID["CCA_baseline"]["h"]   # um
R_MID_MM = (d0 / 2.0 + h0 / 2.0) * 1e-3     # lumen radius + half wall -> mid-wall, mm
P_H_KPA = MAP_CCA * MMHG_TO_KPA

# ---- Bellini et al. (2014) mouse CCA: mass fractions + deposition stretch ---
PHI_ELASTIN, PHI_MEDIA_SMC, PHI_COLLAGEN = 0.249, 0.279, 0.458   # sum 0.986; renormalized below
G_CIRC_FIBER = 0.5 * (1.07 + 1.09)          # circumferential collagen+SMC family, shared value

# ---- Bersi/Ferruzzi/Eberth/Gleason/Humphrey (2014) Table 2, control fit -----
BERSI_CONTROL = dict(elastin_c=8.126, collagen_c1=4.782, collagen_c2=0.041)
BERSI_BANDING_35_56D = dict(elastin_c=4.025, collagen_c1=9.920, collagen_c2=11.579)   # reported only, NOT used as model input -- see docstring

_generic = {c.name: c for c in default_constituents()}   # only for the 2 fields no paper reports (elastin G; turnover k_d/gain)
_phi_sum = PHI_ELASTIN + PHI_MEDIA_SMC + PHI_COLLAGEN
model = Model(
    R=R_MID_MM,
    P_h=P_H_KPA,
    constituents=[
        Constituent("elastin", phi0=PHI_ELASTIN / _phi_sum, G=_generic["elastin"].G,
                    law=NeoHookean(c=BERSI_CONTROL["elastin_c"]),
                    k_d=0.0, gain=0.0, degradable=True),
        Constituent("collagen_smc", phi0=(PHI_MEDIA_SMC + PHI_COLLAGEN) / _phi_sum, G=G_CIRC_FIBER,
                    law=FungFiber(c1=BERSI_CONTROL["collagen_c1"], c2=BERSI_CONTROL["collagen_c2"]),
                    k_d=_generic["collagen"].k_d, gain=_generic["collagen"].gain),
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
    print("Mean-pressure-only stimulus (equilibrated CMM, closest-available real mouse-CCA parameters)")
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
a FIXED, real mouse-CCA material law (Bersi/Ferruzzi/Eberth/Gleason/Humphrey
2014 control fit) and real mouse-CCA composition (Bellini/Ferruzzi/Roccabianca/
Di Martino/Humphrey 2014), predicts (radius fold = lambda*, wall fold =
mass/lambda* from the equilibrated CMM); the Table1 columns are what Eberth
actually measured. Both banded carotids have LOWER MAP than baseline CCA, so a
mean-pressure-only stimulus predicts mild THINNING in both (wall fold < 1).
Table 1 shows the opposite in both cases -- pronounced THICKENING (RCCA-B wall
fold 3.55x, LCCA-B 1.68x) and diverging radius (RCCA-B dilates, LCCA-B
constricts).

A second, independent line of evidence points the same way. The equilibrated
CMM (like every theory in this package) assumes each constituent's intrinsic
stress-stretch behavior is FIXED over time -- only mass and natural
configuration evolve. But the SAME source paper's own "35-56 day post-banding"
fit is a genuinely different material:
    circumferential collagen+SMC (c1, c2):  control {BERSI_CONTROL['collagen_c1']:.3f}, {BERSI_CONTROL['collagen_c2']:.3f}  ->  banded {BERSI_BANDING_35_56D['collagen_c1']:.3f}, {BERSI_BANDING_35_56D['collagen_c2']:.3f}
    elastin c (kPa):                        control {BERSI_CONTROL['elastin_c']:.3f}  ->  banded {BERSI_BANDING_35_56D['elastin_c']:.3f}
A ~280x jump in the stiffening exponent c2 is not a fixed-law CMM predicting
"more of the same material" -- it is new collagen being laid down differently
(Eberth 2009 itself reports the new RCCA-B collagen fibers as thinner and less
ordered; fiber angle shifts 29.4 deg -> 37.1 deg). No fixed-material-law,
mean-stress-only theory in this package can produce that by construction.

Together: mean pressure (hence mean intramural stress) cannot be the sole
driver of this remodeling -- consistent with Eberth's own finding (wall
thickness correlates with pulse pressure, r*=0.632, but NOT mean pressure,
r*=0.020 NS). This motivates adding a flow/WSS/pulsatility channel (svGrowth's
existing wss gain, or a new PI-based one) as the next step.""")
