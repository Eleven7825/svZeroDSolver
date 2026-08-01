# Growth & remodeling (G&R) extension — plan & progress log

## Biological question

Motivated directly by Eberth et al. (2009)'s own framing (and title): when a
vessel remodels its wall, what mechanobiological signal are the cells actually
responding to — the **pulsatility** of the loading (PI, pulse pressure, WSS
amplitude), the **time-averaged (mean)** value (MAP, mean flow, mean WSS), both,
or some other combined metric? Our already-calibrated 0D hemodynamics model
(`two_group_calibration.md`) gives us per-vessel Q(t)/P(t)/τ_w(t) time series
(`make_gr_timeseries_figure.py`) to drive a G&R model and test this against
Eberth's own measured remodeled geometry (Table 1) as the validation target.

### The evidence from the paper itself (Discussion, p.2016–2019)

- Wall thickness correlates with **pulse pressure** (r\* = 0.632) but **not**
  mean pressure (r\* = 0.020, NS).
- Carotid caliber correlates slightly better with **pulsatile flow**
  (r\* = 0.915) than mean flow (r\* = 0.834).
- They cite an in-vitro finding that ~3000 endothelial genes respond to *mean*
  flow, while only ~232 respond specifically to *cyclic* flow — mean and
  pulsatile stimuli are at least partially separable channels, biologically.

### The sharp test case in our own calibrated hemodynamics

Both banded carotids have **lower MAP than baseline CCA** (RCCA-B 86.5, LCCA-B
76.8, vs CCA 92 mmHg). A mean-pressure/mean-stress-only stimulus therefore
predicts wall **thinning** in both. Table 1 shows the opposite for RCCA-B
(upstream of the band, PI stays high): dramatic thickening and dilation
(wall 24.8→88.1 µm, diameter 484→591 µm). LCCA-B (downstream, PI damped)
thickens less and its lumen shrinks (→410 µm). This asymmetry is our in-silico
analog of the paper's own correlation numbers above, and a clean pass/fail test
for any stimulus hypothesis.

## Tools available

Two sibling repos, one level up from `svZeroDSolver`:

- **`Growth-Remodeling`** — teaching repo, 4 arterial G&R theories on one shared
  parameter set (`src/gr/`): kinematic growth, full constrained mixture (CMM),
  homogenized CMM (fast, ODE-based), equilibrated CMM (fastest, algebraic
  fixed-point only). **All four are driven by mean circumferential (intramural)
  Cauchy stress only** — no WSS term, no pulsatility term anywhere (confirmed by
  grepping `docs/`, `src/gr/`, `configs/`).
- **`svGrowth`** (`/home/shiyi/projects/svGrowth`) — research-grade biaxial
  thin-wall CMM (Stanford CBCL), cross-validated against `gr`'s full CMM in
  `Growth-Remodeling/comparison/`. Richer stimulus architecture:
  `k_alpha = k_alpha,h * (1 + sum_i K_i * delta_i^2)`, with named, independently
  gained stimuli (`intramural_stress`, `constituent_stress`, `wss`,
  `inflammation`) — already separates an endothelial/WSS channel (regulates
  lumen radius) from an intramural-stress/SMC-fibroblast channel (regulates wall
  thickness). Still no native PI/pulsatility channel, but the `sum_i K_i delta_i^2`
  architecture is built to accept a new named stimulus — the natural place to
  add one.

Neither repo has a pulsatility stimulus out of the box: that's the gap we're
filling, not something to discover already implemented.

## Where the constitutive parameters come from

Neither `gr` nor svGrowth ships mouse-specific parameters:

- `gr`'s `parameters.py` uses generic "large elastic artery" numbers — its own
  docstring says they're "in the range used by Humphrey, Cyron, Latorre and
  co-workers (order-of-magnitude, teaching-grade; **not a fit to any one
  animal**)." Reference scale `R = 1.0 mm`, `P_h = 13.3 kPa (~100 mmHg)` reads
  as human/large-artery scale, not a ~250 µm-radius mouse carotid.
- svGrowth ships one example, `latorre2018.yaml`, a generic "Cerebral artery"
  (`P_h = 14.18 kPa`, `Q_h = 1 mL/s`) reproducing Latorre & Humphrey (2018)'s
  own theoretical example — also not mouse, not carotid.
- **Eberth (2009) itself doesn't give us constituent parameters either.**
  Their Methods (p. 2012) compute mean wall stress from the simple
  thin/thick-wall formula `sigma_theta = P*r_i/(r_o - r_i)`, not a fitted
  strain-energy model. Their constituent-level data is histology only —
  elastin/collagen/GAG area fractions and collagen-to-elastin ratio (Fig. 7),
  reported **separately for RCCA-B and LCCA-B**, but not stiffness parameters.

**The two real sources, both genuine mouse common carotid artery (CCA) data:**

- Bersi, Ferruzzi, Eberth JF, Gleason, Humphrey (2014), *"Consistent
  Biomechanical Phenotyping of Common Carotid Arteries from Seven Genetic,
  Pharmacological, and Surgical Mouse Models,"* Ann Biomed Eng 42(6):1207–1223.
  Same J.F. Eberth, same aortic-banding mouse model, with an actual
  four-fiber-family constitutive fit at matched timepoints (35–56 days
  post-band ≈ Eberth 2009's 5–8 week chronic point). Table 2:

  | parameter | control | 35–56 day banding |
  |---|---|---|
  | elastin `c` (kPa) | 8.126 | 4.025 |
  | circumferential collagen+SMC `c1` (kPa) / `c2` | 4.782 / 0.041 | 9.920 / **11.579** |
  | fiber angle `α0` (°) | 29.4 | 37.1 |

- Bellini, Ferruzzi, Roccabianca, Di Martino, Humphrey (2014), *"A
  Microstructurally Motivated Model of Arterial Wall Mechanics with
  Mechanobiological Implications,"* Ann Biomed Eng 42(3):488–502. A **10-week
  C57BL/6 wild-type mouse CCA** (homeostatic MAP ≈ 93 mmHg — a different,
  healthy cohort, not the banding model, but same species/vessel), reporting
  what Bersi 2014 doesn't: mass fractions and deposition stretches.

  | quantity | value | note |
  |---|---|---|
  | elastin mass fraction `phi_e` | 0.249 | histology, N=10 images |
  | media (SMC + circumferential collagen) `phi_m` | 0.279 | — |
  | total collagen `phi_c` | 0.458 | sum with above ≈ 0.986 (renormalized to 1) |
  | circumferential collagen deposition stretch `G` | 1.07–1.09 | midpoint 1.08 used |
  | SMC deposition stretch | = circumferential collagen's | their model sets these literally equal |

**Why the model now has 2 constituents, not `gr`'s usual 3.** Both mouse-CCA
papers above treat circumferential collagen and smooth muscle as **one
mechanical family, not two** — Bellini's own notation sets their deposition
stretches equal; Bersi's four-fiber-family fit reports a single `(c1,c2)` pair
for "circumferential collagen+SMC" combined, never SMC alone. Forcing `gr`'s
usual separate elastin/collagen/smc split here would mean inventing an
uncited collagen-vs-SMC breakdown and double-counting that family's stiffness
across two constituents. So the model built here has exactly 2 constituents —
elastin, and one turnover fiber ("collagen_smc") representing the combined
family — which `gr`'s math (`parameters.py`, `equilibrated_cmm.py`) already
supports for any number of turnover constituents, no code change needed.

**Full parameter table used in `make_gr_meanstress_prototype.py`:**

| field | value | source |
|---|---|---|
| `R` (mid-wall radius) | 0.2544 mm | our own calibrated Eberth CCA baseline (`make_configs.CAROTID["CCA_baseline"]`: d=484, h=24.8 µm @ MAP) |
| `P_h` (homeostatic pressure) | 12.266 kPa (= 92 mmHg) | our own calibrated CCA MAP (`make_control_group.MAP`) |
| elastin `phi0` | 0.2525 | Bellini et al. 2014 (0.249, renormalized) |
| elastin `G` | 1.40 | **`gr` generic default — flagged.** Bellini's elastin deposition stretch is an anisotropic (circumferential/axial/radial) triad for a structurally different, anisotropic elastin law; not compatible with `gr`'s isotropic single-scalar elastin formulation. No compatible mouse-CCA value found. |
| elastin `c` (NeoHookean) | 8.126 kPa | Bersi et al. 2014, control fit |
| elastin `k_d`, `gain` | 0, 0 | standard non-turnover assumption (elastin isn't renewed in the adult) |
| collagen+SMC `phi0` | 0.7475 | Bellini et al. 2014 (0.279+0.458, renormalized) |
| collagen+SMC `G` | 1.08 | Bellini et al. 2014, circumferential family (range 1.07–1.09, midpoint used) |
| collagen+SMC `c1, c2` (Fung) | 4.782, 0.041 | Bersi et al. 2014, control fit |
| collagen+SMC `k_d`, `gain` | 1/20 day, 1.0 | **`gr` generic default — flagged.** Neither paper reports a turnover rate or mechanosensitivity gain. Doesn't affect the result below: the equilibrated CMM's algebraic equilibrium equation doesn't use `k_d`/`gain` at all (they only matter for *transient* dynamics — homogenized/full CMM, a later step). |

**We use ONLY the control fit as the model's fixed material law — for all
three vessels (CCA / RCCA-B / LCCA-B), never the post-band fit.** A
constrained-mixture theory's premise is that each constituent's intrinsic
stress-stretch behavior is fixed over time; only mass and natural
configuration evolve under a stimulus. Swapping in the post-band fit as a
different material law for the banded rows would hand the model the answer
instead of testing whether mean pressure predicts it. The post-band fit is
used only as an independent check (see Step 1 below), never as a solve input.

One more caveat carried forward: both real fits above are for combined/pooled
groups, not split RCCA-B vs. LCCA-B — Eberth's own Fig. 7 area fractions
remain the only per-side compositional data we have, not yet used here.

## Staged plan

1. **Fast prototype (mean-stress only)** — quick go/no-go: does the standard
   assumption (mean intramural stress alone) explain the RCCA-B/LCCA-B
   asymmetry? **Done, see below — no.**
2. **Extend svGrowth's stimulus** with a new named PI/pulsatility-amplitude (or
   WSS-amplitude) channel alongside its existing `wss` and `intramural_stress`
   gains, each independently fittable. Not started.
3. **Coupling**: quasi-static one-way — freeze geometry, run one cardiac cycle
   in svZeroDSolver, extract MAP/Q̄/WSS̄/PI per vessel as that G&R step's
   forcing, advance the G&R state, repeat (justified by the huge timescale
   separation: cycle ~0.15 s vs G&R over weeks). Two-way feedback (remodeled
   radius → bed R via Poiseuille, R∝1/a⁴) is a later refinement. Not started.

## Step 1 result: mean-pressure-only stimulus — REJECTED

`make_gr_meanstress_prototype.py` runs the equilibrated CMM
(`Growth-Remodeling/src/gr/equilibrated_cmm.py`) with a sustained mean-pressure
insult (load factor = calibrated group MAP / CCA MAP) on the fixed, real
mouse-CCA 2-constituent model above, against Table 1's "@ MAP" geometry
fold-changes:

| group | MAP factor | model: radius fold | model: wall fold | Table1: diam fold | Table1: wall fold |
|---|---|---|---|---|---|
| RCCA-B | 0.940 | 0.973 | 0.915 | **1.221** | **3.552** |
| LCCA-B | 0.835 | 0.932 | 0.778 | 0.847 | **1.677** |

Mean-pressure-only predicts **thinning in both vessels** (wall fold < 1, since
MAP is lower than baseline in both groups) and mild **shrinkage** in both
(radius fold < 1). Table 1 shows pronounced **thickening** in both (3.55× /
1.68×). Wrong sign on wall thickness for both vessels — mean intramural stress
alone is decisively ruled out as the sole G&R stimulus for this dataset,
consistent with the paper's own correlation finding. (Radius direction now
happens to match for LCCA-B — both predicted and actual shrink, though by very
different magnitudes, 0.932 vs 0.847 — see "note on parameter sensitivity"
below on why the radius-fold magnitude/sign is the less robust of the two
numbers and shouldn't be over-read either way.)

**A second, independent line of evidence points the same way.** Every theory
here assumes a constituent's material behavior is fixed over time — only mass
and configuration evolve. But Bersi et al. (2014)'s own "35–56 day banding" fit
(reported, not used as a model input — see above) is a genuinely different
material from their control fit: the circumferential stiffening exponent `c22`
jumps ~280× (0.041 → 11.579), and fiber angle shifts 29.4°→37.1°. That is not
what a fixed-material-law, mean-stress-driven CMM can produce by construction
— it's new collagen being deposited differently (consistent with Eberth 2009's
own observation that the new RCCA-B collagen fibers are thinner and less
ordered). This is further, independent evidence that something beyond mean
pressure is driving the remodeling.

### Note on parameter sensitivity (why we trust the *sign*, not the magnitude)

When this question came up mid-project — does the conclusion depend on the
specific constituent parameters? — we checked it two ways rather than assert
it:

1. **An exact structural identity, independent of any constituent parameter.**
   At equilibrium, "growth has stopped" means the mixture stress has returned
   exactly to `sigma_bar_h` (the theory's own definition of equilibrium); combined
   with the Laplace balance this forces `mass = gamma * lambda^n` (`n=2` for the
   artery), so `thickness_fold = gamma * lambda*`. This holds for *any*
   constituent parameters — it isn't something that happens to be true for
   the numbers in the file.
2. **What genuinely depends on the parameters is `lambda*` itself.** We
   swept 500 random draws of mass fractions, deposition stretches,
   stiffnesses, and gains at `gamma = 0.835` (our LCCA-B factor) and found
   **zero** cases where `thickness_fold > 1` despite `gamma < 1` — i.e., the
   *sign* of "mean-pressure-only predicts thinning" held in every draw, even
   though `lambda*` itself (and hence the exact magnitude) varied a lot,
   including flipping which direction the radius moved.

So: trust the qualitative conclusion (mean pressure alone predicts thinning,
Table 1 shows thickening) — that's robust. Don't over-read the exact model
percentages as calibrated mouse-carotid magnitudes; they weren't fit to
reproduce Table 1, only to use a real, cited material law instead of a
generic placeholder.

## Next up

Extend svGrowth's kinetics with a flow/WSS or PI-based channel, driven by our
RCCA-B/LCCA-B/CCA Q̄, WSS̄, and PI (already computed in
`make_gr_timeseries_figure.py`), and check whether that recovers Table 1's
thickening and the diameter divergence.

## References

- Eberth JF, Taucer AI, Wilson E, Humphrey JD. Importance of pulsatility in
  hypertensive carotid artery growth and remodeling. *J Hypertens*
  27(10):2010–2021, 2009.
- Bersi MR, Ferruzzi J, Eberth JF, Gleason RL Jr, Humphrey JD. Consistent
  biomechanical phenotyping of common carotid arteries from seven genetic,
  pharmacological, and surgical mouse models. *Ann Biomed Eng*
  42(6):1207–1223, 2014.
- Bellini C, Ferruzzi J, Roccabianca S, Di Martino ES, Humphrey JD. A
  microstructurally motivated model of arterial wall mechanics with
  mechanobiological implications. *Ann Biomed Eng* 42(3):488–502, 2014.
- Ferruzzi J, Bersi MR, Humphrey JD. Biomechanical phenotyping of central
  arteries in health and disease: advantages of and methods for murine models.
  *Ann Biomed Eng* 41(7):1311–1330, 2013.
- Cyron CJ, Aydin RC, Humphrey JD. A homogenized constrained mixture (and
  mechanical analog) model for growth and remodeling of soft tissue. *Biomech
  Model Mechanobiol* 15:1389–1403, 2016.
- Latorre M, Humphrey JD. Mechanobiological stability of biological soft
  tissues. *J Mech Phys Solids* 125:298–325, 2019 (equilibrated CMM: 2018
  preprint/earlier formulation used by `Growth-Remodeling/src/gr/equilibrated_cmm.py`).
