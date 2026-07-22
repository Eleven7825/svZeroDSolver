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

1. **Upstream amplification is capped at ~1.2–1.3× in 0D** (paper needs 2.7×).
   The missing factor is the **reflected pressure wave** returning from the
   stenosis to the RCCA takeoff — a wave-transit (1D) effect a lumped element
   has no representation for.
2. **The downstream carotid always ends up *below* baseline in 0D** (LCCA
   0.6–0.9), but the paper shows it *above* (1.16 → 1.65). For the *damped*
   downstream vessel to be more pulsatile than the healthy baseline, the whole
   systemic pulse must rise (afterload/reflection) — which prescribed flow
   doesn't produce and the elastance heart *offsets* (higher afterload → lower
   stroke volume).

## Why the contractile heart specifically did not help

- It made the **baseline worse**: a stiff mouse ventricle (Ees ≈ 9 mmHg/µL)
  ejects very impulsively → peakier carotid flow → higher PI (3.2 vs 2.2).
- Its hoped-for advantage — afterload-driven pulse amplification when banded —
  is **cancelled** because the higher afterload reduces stroke volume.
- Filling had to be tuned (preload 10 mmHg, low mitral resistance) so the stiff
  LV fills in the short mouse diastole (EDV 52, SV 33 µL, CO 0.20, MAP 93).

## Conclusion

Across **four independent routes** (bed sweep, waveform, contractile heart, and
the combined heart+bed+band-inertance+HR sweep) **no physiological 0D
parameterization matches both the healthy baseline (1.16) and the banded
carotids (3.11/1.65)**. The band's 0D amplification (~1.3×) is intrinsically
short of the required 2.7×, and the downstream vessel cannot rise above baseline
without a systemic pulse increase the lumped model does not generate.

**The missing physics is wave reflection off the stenosis — inherently 1D.** This
is exactly the boundary the source note flags: *"move to 1D only if wave-reflection
timing … becomes the question of interest."* For the full baseline→banded
*amplitude*, a 1D wave-propagation arch segment is required; 0D remains valid for
the pulse-pressure redistribution, MAP, mean flows, and the qualitative story.

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
