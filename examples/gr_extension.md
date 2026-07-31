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

`make_gr_meanstress_prototype.py` runs the equilibrated CMM (`Growth-Remodeling
/src/gr/equilibrated_cmm.py`) with a sustained mean-pressure insult (load factor
= calibrated group MAP / CCA MAP), generic teaching-grade constituent
parameters (not fit to mouse carotid tissue — so we compare **fold-change**,
which doesn't depend on those absolute constants, not absolute microns), against
Table 1's "@ MAP" geometry fold-changes:

| group | MAP factor | model: radius fold | model: wall fold | Table1: diam fold | Table1: wall fold |
|---|---|---|---|---|---|
| RCCA-B | 0.940 | 0.998 | 0.938 | **1.221** | **3.552** |
| LCCA-B | 0.835 | 0.994 | 0.830 | 0.847 | **1.677** |

Mean-pressure-only predicts mild **thinning** in both vessels (wall fold < 1,
since MAP is lower than baseline in both groups) and essentially unchanged
radius. Table 1 shows pronounced **thickening** in both (3.55× / 1.68×) with
radius diverging in opposite directions. Wrong sign on wall thickness for both
vessels — mean intramural stress alone is decisively ruled out as the sole G&R
stimulus for this dataset, consistent with the paper's own correlation finding.

## Next up

Extend svGrowth's kinetics with a flow/WSS or PI-based channel, driven by our
RCCA-B/LCCA-B/CCA Q̄, WSS̄, and PI (already computed in
`make_gr_timeseries_figure.py`), and check whether that recovers Table 1's
thickening and the diameter divergence.
