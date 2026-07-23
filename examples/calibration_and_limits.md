# Calibration attempts & the 0D limitation

This note records **what we tried** to make the Eberth aortic-banding 0D model
match *both* the healthy baseline and the banded carotids, and **why a 0D
lumped model cannot do both**. It backs the "Known limitation" in
[`eberth_aortic_band_model.md`](./eberth_aortic_band_model.md) §8 with the actual
experiments.

## The target

Eberth Table 1 (all velocity-based PI):

| | baseline CCA | RCCA-B (upstream) | LCCA-B (downstream) |
|---|---:|---:|---:|
| **PI** | 1.16 | 3.11 | 1.65 |
| pulse pressure (mmHg) | 42 | 56 | 27 |

The hard part is the **amplitude of the change**: banding must lift the upstream
carotid PI **2.7×** (1.16 → 3.11) and *also* lift the downstream one *above*
baseline (1.16 → 1.65), at nearly constant MAP.

What the 0D model reproduces well (all variants): the **pulse-pressure**
redistribution (≈42 → ≈60/≈36), **MAP similarity** across beds, mean carotid
flows in range, and the qualitative healthy → banded → remodeled story. The
gap is purely in the **flow/velocity PI amplitude**.

## What we tried (and what each did)

| # | Attempt | Effect on baseline PI | Effect on banded PI | Verdict |
|---|---------|----------------------|---------------------|---------|
| 1 | **Terminal-bed sweep** (Rp:Rd split, bed C, total R, aortic Ca) | any knob that lowers baseline… | …lowers banded by the same factor | Coupled — cannot decouple the two |
| 2 | **Velocity vs flow** (area pulsation, distensibility 3.76e-3/mmHg) | 2.23 → 1.99 (~5%) | negligible | Not the cause |
| 3 | **Ejection waveform** (systolic fraction) | broadening → 1.16 ✓ | but banded collapses to ~1.7 | Non-physiological fix; breaks banded |
| 4 | **Hartley waveform** (real mouse aortic Doppler, ET 38%) | 2.23 (physiological, peaky) | **3.0 / 1.6 ✓** (banded match) | Best 0D; baseline still ~2× high |
| 5 | **Contractile elastance heart** (mouse LV: Ees 9 mmHg/µL, EDV 45, ESV 18 µL) | **worse**, 3.2–4.9 | amplification 1.2–1.3× | Impulsive ejection worsens baseline |
| 6 | **Combined**: elastance heart × bed C × band inertance L × HR | plateaus ~1.5 (rigid bed) | with baseline low, RCCA maxes **~2.3** (even L=3.0, ~100×) | Confirms the cap |
| 7 | **State-specific HR** (7.17 baseline vs 6.09 banded) | small | small | Legitimate but second-order |

### Key numbers from the combined sweep (attempt 6)

Lowering bed compliance pulls the baseline down but **plateaus ~1.5** (can't reach
1.16); and once the baseline is low, the band can no longer amplify:

| bed C | baseline PI | banded RCCA (best, L=3.0) | banded LCCA |
|------:|------------:|--------------------------:|------------:|
| ×1.0 | 4.9 | (baseline already high) | — |
| ×0.25 | 3.2 | ~4.2 | ~1.2 |
| ×0.05 | 1.9 | **~2.3** (target 3.11) | ~0.7 (target 1.65) |
| ×0.02 | 1.6 | ~2.0 | ~0.6 |

The apparent "RCCA hits 3.11" at high bed C is an artifact of an **inflated
baseline** (3.2), not band amplification. When the baseline is correct, RCCA
tops out ~2.3.

## Two mechanistic facts

1. **Upstream amplification is capped at ~1.2–1.3× in 0D**, whereas the paper
   needs 2.7× (1.16 → 3.11).
2. **The downstream carotid always ends up *below* baseline in 0D** (LCCA
   0.6–0.9), whereas the paper shows it *above* (1.16 → 1.65).

## Why the contractile heart specifically did not help

- It made the **baseline worse**: a stiff mouse ventricle (Ees ≈ 9 mmHg/µL)
  ejects very impulsively → peakier carotid flow → higher PI (3.2 vs 2.2).
- Its hoped-for advantage — afterload-driven pulse amplification when banded —
  is **cancelled** because the higher afterload reduces stroke volume.
- Filling had to be tuned (preload 10 mmHg, low mitral resistance) so the stiff
  LV fills in the short mouse diastole (EDV 52, SV 33 µL, CO 0.20, MAP 93).

## Estimated size of a 1D wave-reflection effect (mouse)

A natural question is whether adding 1D wave propagation would supply the missing
amplification. An order-of-magnitude estimate for the mouse:

- Reflection coefficient at the band Γ = (Z_throat − Z_ao)/(Z_throat + Z_ao). With
  Z ∝ 1/A and Aₛ/A₀ ≈ 0.21, Z_throat/Z_ao ≈ 4.8 → **Γ ≈ 0.65** (a strong reflector).
- Mouse aortic pulse-wave velocity c ≈ 3.5 m/s. The RCCA/innominate takeoff sits
  essentially at the band's proximal face (L ≈ 0.5–1 mm), so the round-trip
  transit time τ = 2L/c ≈ **0.6 ms**, versus systole ≈ 60 ms → **τ/systole ≈ 1 %**.
- Fundamental wavelength λ = c/f ≈ 3.5/6.09 ≈ **0.57 m (57 cm)**, versus an arch of
  ~1 cm → **L/λ ≈ 0.02**. Only the highest pulse harmonics (~50 Hz, λ ≈ 6 cm) reach
  L/λ ≈ 0.1–0.2.

So although the band reflects strongly, the incident and reflected waves are
essentially **in phase** (the long-wavelength / lumped regime). The estimated
incremental 1D gain on the RCCA PI is therefore **~10–25 % (factor ≈ 1.1–1.25×)**,
concentrated in the sharpest systolic harmonics — **not** the ~2× needed to reach
3.11. Consistent with this, the 0D model already reproduces the banded *pressure*
pulse (PP ≈ 56–60 vs paper 56), i.e. the band's quasi-static pressure rise is
already captured.

## Conclusion (facts)

- Across four independent routes (bed sweep, waveform, contractile heart,
  combined sweep) **no physiological 0D parameterization matches both the healthy
  baseline (1.16) and the banded carotids (3.11/1.65)**.
- The band's 0D upstream amplification caps at **~1.3×**; the paper requires
  **2.7×**. The downstream carotid stays **below** baseline in 0D, but rises above
  it in the paper.
- For the mouse geometry, an added **1D wave-reflection effect is estimated at
  only ~10–25 %** — well short of the missing factor.
- 0D remains quantitatively valid for the **pulse-pressure redistribution, MAP,
  mean flows**, and the qualitative healthy → banded → remodeled story. The
  residual *flow-PI amplitude* gap is not accounted for by any factor tested here.

## Band-constancy assumption (sensitivity)

The models hold the band identical from acute to chronic (a rigid external
spacer), letting only the carotids remodel. The band region is living aortic
wall, though, so its effective throat could change over 5–8 weeks (neointima,
peri-band fibrosis, proximal dilation) — most plausibly tightening. Sweeping the
band throat (loose → tight, recomputing S, L, R with consistent stenosis physics)
on the chronic model:

| throat | area stenosis | S | RCCA-B PI | LCCA-B PI | MAP | band ΔP̄ |
|--------|--------------:|----:|----------:|----------:|----:|--------:|
| 528 µm (loose) | 64% | 12 | 2.39 | 1.84 | 94 | 2 |
| 467 µm (loose) | 72% | 26 | 2.59 | 1.73 | 96 | 4 |
| 406 µm (nominal) | 78% | 54 | 2.88 | 1.53 | 98 | 7 |
| 345 µm (tight) | 84% | 119 | 3.15 | 1.22 | 104 | 14 |
| 284 µm (tight) | 89% | 291 | 3.15 | 0.91 | 117 | 28 |

*(paper chronic: RCCA-B 3.11, LCCA-B 1.65, MAP similar)*

Facts from the sweep:

- **RCCA-B and LCCA-B move in opposite directions with band tightness.** A loose
  band matches LCCA-B (1.73–1.84 vs 1.65) but leaves RCCA-B too low; a tight band
  matches RCCA-B (3.15 vs 3.11) but collapses LCCA-B (~1.0). No single tightness
  matches both carotids.
- **RCCA-B saturates at ~3.15** even at 89% stenosis.
- **MAP excludes the tight end**: at 84–89% stenosis MAP rises to 104–117 mmHg,
  contradicting the paper's near-constant MAP.
- So relaxing band-constancy does not let 0D match both carotids and MAP; the
  three-way RCCA-B ↔ LCCA-B ↔ MAP tension persists at every band tightness. The
  nominal 78% band is the best compromise across the three targets.

## Why remodeling barely changes the pressure curve

The acute and chronic states give nearly identical carotid pressure (node A
77/137, MAP ~99, PP 60 in **both**) despite the ~3.5× RCCA-B wall thickening.
This is **not** self-compensating remodeling — the model is passive and linear,
with no autoregulation or feedback. The pressure is simply **insensitive to the
carotid wall geometry**, for two structural reasons:

- **The carotid segment is a negligible impedance in its own path.** The head RCR
  bed is Rp+Rd ≈ 5294 (held fixed); the carotid *segment* R is only **2.3%
  (acute) → 0.9% (chronic)** of that, and the wall C (~1e-6) is ~1% of the
  aortic/bed compliances. So remodeling the segment (R 124 → 47) shifts the
  carotid-path impedance by ~1% — the fixed distal bed sets the carotid's flow
  and pressure, not the wall geometry.
- **The carotid is a minor side branch** (~9% of cardiac output). Node-A pressure
  is set by the main path (heart → aorta → band → systemic, ~91% of CO), which is
  identical in both states; the carotid only "taps" node A and draws little
  current.

Caveat: only the carotid **artery-wall segment** remodels here (as Eberth
measured); the distal cerebral bed is held fixed. Since the bed dominates the
carotid impedance, segment remodeling has little leverage — if the cerebral bed
itself remodeled, the effect would be larger. This is the same reason acute ≈
chronic throughout the results.

## Reproducing the experiments

- `make_configs.py` — the prescribed-flow 0D models (prebanding/acute/chronic).
- `calibrate_stenosis.py` — band stenosis calibration.
- `make_elastance_config.py` — the **contractile-heart variant** (attempt 5/6):
  `preload → mitral valve → ChamberElastanceInductor (mouse LV) → aortic valve →
  aorta`, tuned to a physiological operating point. Writes
  `eberth_elastance_{prebanding,acute,chronic}.json`.
- LV literature values (attempt 5): Ees ≈ 8–10 mmHg/µL (Pacher 2008 *Nat Protoc*;
  Sci Rep 2019 sham 8.4), EDV ≈ 45 µL, ESV ≈ 18 µL, SV ≈ 27 µL, ESP ≈ 105,
  LVEDP ≈ 5 mmHg, filling ≈ 5 mmHg. TAC raises afterload ~2× (Sci Rep 2019).
- `make_pv_loops.py` — plots the elastance-heart LV pressure–volume loops
  (`figures/fig_pv_loops.png`). The loops confirm the afterload response: banding
  shifts the loop up-and-right (ESP 147→162 mmHg, SV 33→31 µL), the classic
  afterload signature.
