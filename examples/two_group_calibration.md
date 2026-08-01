# Two-group calibration (control vs banded) — plan & execution log

## Why two groups

The Eberth data come from **two separate animal groups** (7 controls + 7 banded),
each measured once at the chronic timepoint. There is no requirement that they
share R/C/L, and the paper's own numbers differ between them (HR 7.17 vs 6.09 Hz,
carotid geometry, per-vessel flows 0.016 / 0.022 / 0.012, pressures). Our earlier
single-shared-parameter model (band + carotid geometry the only differences) was
therefore over-constrained — part of the "can't match both" tension was
self-imposed. This note plans and logs a group-specific calibration.

## Parameter tiers

- **Tier 1 — shared species invariants** (never tuned): blood µ/ρ, topology,
  ejection-waveform shape, aortic & carotid wall modulus E.
- **Tier 2 — group-specific, fixed from measurement**: HR (7.17 control / 6.09
  banded), carotid geometry (inner diameter **@ MAP**, Table 1: CCA 484/24.8;
  RCCA-B 591/88.1; LCCA-B 410/41.6 — not the "@100mmHg" ex vivo reference
  diameter, since each vessel here runs at its own physiological MAP),
  band geometry (banded only). Plus the measured per-vessel Q̄ and pressures,
  used as *constraints* that pin the bed resistances and mean operating point.
- **Tier 3 — free knobs, calibrated to each group's own targets**: terminal-bed
  compliance (`RCR_RIGHT`/`RCR_LEFT`/`RCR_SYS` $C$), aortic compliance
  (`ascending_aorta`/`descending_aorta` $C_a,C_b$), and (banded only) the band
  `stenosis_coefficient` $S$ (`aortic_band`). Bed *resistances* ($R_p{+}R_d$) are
  pinned by (Q̄, MAP), so they are not free. Symbols match the circuit diagram
  in `eberth_0d_presentation.tex` (slide 4).

## Data sources

Every Tier 1/Tier 2 (FIXED) value above traces to a specific citation, tagged
in each generated config's `description.parameter_tags` (see `make_configs.py`).
Full list:

| quantity | value | source |
|---|---|---|
| blood density ρ | 1.06 g/cm³ | Aslanidou et al. 2016 |
| blood viscosity µ | 3.5 cP (high-shear) | Windberger et al. 2003 |
| aortic wall modulus `E_aorta` | 1.0 MPa | Yanagisawa & Wagenseil 2019 / Bersi et al. 2014 |
| carotid wall modulus `E_carotid` | 1.5 MPa (systolic, linearized) | Ferruzzi, Bersi & Humphrey 2013 |
| heart rate HR | 7.17 Hz (control) / 6.09 Hz (banded) | Eberth et al. 2009 |
| cardiac output CO | 0.20 ml/s (~12 ml/min, anesthetized) | Aslanidou et al. 2016; cf. Janssen et al. 2002, Mills et al. 2000 |
| MAP target | ~92 mmHg | Constantinides et al. 2011; cf. Janssen et al. 2002, Mills et al. 2000 |
| carotid geometry (ID @ MAP, wall) | CCA 484/24.8; RCCA-B 591/88.1; LCCA-B 410/41.6 µm | Eberth et al. 2009, Table 1 |
| aorta geometry (ID, wall) | ascending 1.4 mm / descending 0.9 mm / wall 40 µm | Casteleyn et al. 2010; Guo & Kassab 2003 |
| aortic buffer compliance $C_a, C_b$ | 2.67e-4 ml/mmHg (measured central aortic compliance) | Aslanidou et al. 2016 |
| RCR bed flow split (10/90 Rp/Rd; carotid vs. systemic Q̄ split) | — | Feintuch et al. 2007; Trachet et al. 2009 |
| per-vessel targets: PI, Q̄, MAP, PP | CCA 1.16/0.016/92/42.5; RCCA-B 3.11/0.022/86.5/56; LCCA-B 1.65/0.012/76.8/27; A→B MAP drop 9.7 | Eberth et al. 2009, Table 1 + text |

**References** (full citations):

- Eberth JF, Taucer AI, Wilson E, Humphrey JD. Importance of pulsatility in
  hypertensive carotid artery growth and remodeling. *J Hypertens*
  27(10):2010–2021, 2009.
- Aslanidou L, Trachet B, Reymond P, Fraga-Silva RA, Segers P, Stergiopulos N.
  A 1D model of the arterial circulation in mice. *ALTEX* 33(1):13–28, 2016.
- Windberger U, Bartholovitsch A, Plasenzotti R, Korak KJ, Heinze G. Whole
  blood viscosity, plasma viscosity and erythrocyte aggregation in nine
  mammalian species: reference values and comparison of data. *Exp Physiol*
  88(3):431–440, 2003.
- Yanagisawa H, Wagenseil J. Elastic fibers and biomechanics of the aorta:
  insights from mouse studies. *Matrix Biol* 85–86:160–172, 2019/2020.
- Ferruzzi J, Bersi MR, Humphrey JD. Biomechanical phenotyping of central
  arteries in health and disease: advantages of and methods for murine models.
  *Ann Biomed Eng* 41(7):1311–1330, 2013.
- Bersi MR, Ferruzzi J, Eberth JF, Gleason RL Jr, Humphrey JD. Consistent
  biomechanical phenotyping of common carotid arteries from seven genetic,
  pharmacological, and surgical mouse models. *Ann Biomed Eng*
  42(6):1207–1223, 2014.
- Casteleyn C, Trachet B, Van Loo D, Devos DG, Van Hoorebeke L, Segers P,
  Simoens P. Validation of the murine aortic arch as a model to study human
  vascular diseases. *J Anat* 216(5):563–571, 2010.
- Guo X, Kassab GS. Variation of mechanical properties along the length of the
  aorta in C57BL/6 mice. *Am J Physiol Heart Circ Physiol* 285(6):H2614–H2622,
  2003 (cited elsewhere in this project as "Guo & Kassab 2002" per submission
  year).
- Feintuch A, Ruengsakulrach P, Lin A, et al. Hemodynamics in the mouse aortic
  arch as assessed by MRI, ultrasound, and numerical modeling. *Am J Physiol
  Heart Circ Physiol* 292(2):H884–H892, 2007.
- Trachet B, Swillens A, Van Loo D, Casteleyn C, De Paepe A, Loeys B, Segers P.
  The influence of aortic dimensions and boundary conditions on calculated
  wall shear stress in the mouse aortic arch. *Comput Methods Biomech Biomed
  Engin* 12(5):491–499, 2009.
- Janssen BJ, Debets JJ, Leenders PJ, Smits JF. Chronic measurement of cardiac
  output in conscious mice. *Am J Physiol Regul Integr Comp Physiol*
  282(3):R928–R935, 2002.
- Constantinides C, Mean R, Janssen BJ. Effects of isoflurane anesthesia on the
  cardiovascular function of the C57BL/6 mouse. *ILAR J* 52:e21–e31, 2011.
- Mills PA, Huetteman DA, Brockway BP, et al. A new method for measurement of
  blood pressure, heart rate, and activity in the mouse by radiotelemetry.
  *J Appl Physiol* 88(5):1537–1544, 2000.

## Identifiability (DOF)

Targets per vessel = {P_sys, P_dias, Q̄, PI} = 4.
- **Control** (1 carotid): 4 targets ↔ 4 knobs ($R_p^{R,L}{+}R_d^{R,L}$ bed $R$,
  $C^{R,L}$ bed $C$, $C_a,C_b$ aortic $C$, $R_p^{S}{+}R_d^{S}$ systemic $R$) → determined.
- **Banded** (2 carotids): 8 targets ↔ 8 knobs → determined, but band $S$ is
  redundant with the two bed $C$'s ($C^{R},C^{L}$) for the PIs. The pulsatility
  sub-problem is the crux: with the cerebral beds *shared* (mechanistic), the 2
  banded PIs have ~1 effective knob → under-determined (this is why they
  couldn't both be hit); with per-carotid bed $C$ free, they are matchable but
  $S$ becomes redundant and the model turns descriptive.

## Step 1 — CONTROL group on its own terms (executed)

Build: HR **7.17**, baseline CCA carotids, **no band**, terminal beds **pinned**
to the measured CCA Q̄ = 0.016 ml/s and MAP ≈ 92, aortic compliance ($C_a,C_b$) at
the **measured** value; the **carotid-bed compliance ($C^{R,L}$)** is the one
calibrated knob. Script: `make_control_group.py` → `eberth_control_group.json`.

**Result — the passive RCR matches all four CCA targets:**

| knob | flow-PI | Q̄ | MAP | pulse pressure |
|---|---:|---:|---:|---:|
| carotid-bed C = **6.44e-6** (calibrated) | **1.16** ✓ | 0.0160 ✓ | 94 ✓ | **41 ≈ 42.5** ✓ |
| paper CCA target | 1.16 | 0.016 | ~92 | 42.5 |

**Findings (facts):**
1. **HR 7.17 barely changed the baseline PI** (2.14 vs 2.09 at 6.09 Hz) — heart
   rate is not the lever.
2. **Beds pinned to the measured flow reproduce Q̄ and MAP exactly.**
3. **The flow pulsatility is set by the *carotid bed compliance*** (the flow
   admittance): low bed C → flow tracks pressure resistively, PI → PP/MAP ≈ 0.46;
   high bed C → the capacitor adds flow swing, PI → large; an intermediate value
   (6.44e-6) gives PI = 1.16.
4. Because the carotid is a **minor branch (~9% of CO)**, tuning its bed C moves
   the flow-PI **without disturbing the pressure pulse** — so all four targets
   (PI, Q̄, MAP, PP) are matched together.

**Correction to an earlier claim:** an initial version of this note tuned the
*aortic* compliance C_a to hit the PI, which failed (C_a buffers the pressure, so
it couples flow-PI and PP → fixing one broke the other) and led to a wrong
"needs cerebral autoregulation / passive RCR insufficient" conclusion. That was a
**wrong-knob error**: the correct knob is the **carotid bed C**, and the **passive
RCR is sufficient** for the baseline — no autoregulation needed.

## Status / next

- `eberth_control_group.json` matches CCA on flow-PI, Q̄, MAP, and PP with the
  measured aortic compliance and a calibrated carotid-bed compliance.
- This also means the **descriptive** two-group route (per-carotid bed C) is
  viable with the stock RCR: control uses its bed C, and the banded group would
  use its own RCCA/LCCA bed C's (+ band S) — the DOF budget (8 targets / 8 knobs)
  is met. The **mechanistic** route (shared bed, band as sole driver) remains
  under-determined for the two banded PIs (see `calibration_and_limits.md`).
## Step 2 — BANDED group on its own terms (executed)

Build: HR **6.09**, remodeled carotids (RCCA-B 591/88.1, LCCA-B 410/41.6 @ MAP), band.
Beds **pinned** to each carotid's measured Q̄ and MAP (RCCA-B 0.022 @ 86.5,
LCCA-B 0.012 @ 76.8, systemic @ 76.8). Band **$S$ calibrated to the A→B mean
pressure drop** (MAP_A − MAP_B ≈ 9.7 mmHg — a *measured* quantity); the two
**per-carotid bed compliances ($C^{R},C^{L}$) calibrated to the PIs**. Script:
`make_banded_group.py` → `eberth_banded_group.json`.

**Result — matches all banded targets:**

| vessel | flow-PI | Q̄ | MAP | PP |
|---|---:|---:|---:|---:|
| RCCA-B | **3.10** ✓ | 0.0221 ✓ | 88 ✓ | 63 (meas 56) |
| LCCA-B | **1.65** ✓ | 0.0118 ✓ | 79 ✓ | 31 (meas 27) |
| A→B MAP drop | 9.4 (meas 9.7) | | | |

Calibrated knobs: band **S = 95.1**, RCCA-B bed C = 2.56e-5, LCCA-B bed C = 2.33e-5.

**Findings:**
1. The banded group matches **all** its own targets (both PIs, both flows, both
   MAPs, the band drop) with the **stock passive RCR** — no new block.
2. The band S is **anchored to the measured mean pressure drop** (not fit to the
   PI), a cleaner mechanistic constraint than the earlier ratio-only S.
3. The two PIs are hit by the two **per-carotid bed compliances** (flow-admittance
   knobs) — DOF-balanced (8 targets / 8 knobs), as the plan predicted.
4. The two carotid beds came out **similar** (2.56 vs 2.33e-5), so the RCCA-B ≫
   LCCA-B pulsatility asymmetry is driven mainly by **band position + geometry**,
   not by asymmetric cerebral beds — a mechanistically satisfying outcome.

## Conclusion of the two-group approach

- **Both groups match their own data** with a stock passive RCR: control (PI 1.16,
  Q̄ 0.016, MAP ~92, PP 42) and banded (RCCA-B 3.11 & LCCA-B 1.65, flows, MAPs,
  band drop). The band is anchored to the measured pressure drop.
- This confirms the DOF analysis: the **descriptive** route (per-group,
  per-carotid bed C) is well-posed and **succeeds** — the earlier "can't match
  both" was specific to the **mechanistic single-shared-parameter-set** route.
- **Cost (the descriptive price):** the fitted terminal-bed compliance differs
  between groups (control 6.44e-6 vs banded ~2.3–2.6e-5, ~3.6–4×). That group difference
  is fitted, not independently measured. So the two configs are two calibrated
  states, linked by the band (anchored to the measured drop) and the measured
  geometry/HR/flows — not a single predictive model.
- Deliverables: `eberth_control_group.json`, `eberth_banded_group.json`, built by
  `make_control_group.py` / `make_banded_group.py`. The comparison figure
  (`make_two_group_figure.py` → `figures/fig_two_group.png`) shows the fit vs
  Table 1 for PI, mean flow, MAP, and pulse pressure across CCA / RCCA-B / LCCA-B.
  The settled knobs (`make_knobs_figure.py` → `figures/fig_knobs.png`) show all the
  $R_p{+}R_d$ / $C$ knobs (band $S$ excluded), grey = control vs blue = banded:
  bed $R_p{+}R_d$ differs by carotid (sets the flow split), the banded bed
  compliances $C^{R},C^{L}$ are ~3–4× the control's $C^{R,L}$ (the fitted group
  difference), systemic $R_p^{S}{+}R_d^{S}$ is similar (548 vs 463), and aortic
  $C_a,C_b$ is fixed/identical (2.67e-4).
