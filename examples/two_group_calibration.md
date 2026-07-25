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
  banded), carotid geometry (CCA 496/24.8; RCCA-B 633/88.1; LCCA-B 482/41.6),
  band geometry (banded only). Plus the measured per-vessel Q̄ and pressures,
  used as *constraints* that pin the bed resistances and mean operating point.
- **Tier 3 — free knobs, calibrated to each group's own targets**: terminal-bed
  compliance, aortic compliance, and (banded only) the band stenosis_coefficient.
  Bed *resistances* are pinned by (Q̄, MAP), so they are not free.

## Identifiability (DOF)

Targets per vessel = {P_sys, P_dias, Q̄, PI} = 4.
- **Control** (1 carotid): 4 targets ↔ 4 knobs (bed R, bed C, aortic C, systemic R) → determined.
- **Banded** (2 carotids): 8 targets ↔ 8 knobs → determined, but band S is
  redundant with the two bed C's for the PIs. The pulsatility sub-problem is the
  crux: with the cerebral beds *shared* (mechanistic), the 2 banded PIs have
  ~1 effective knob → under-determined (this is why they couldn't both be hit);
  with per-carotid bed C free, they are matchable but S becomes redundant and the
  model turns descriptive.

## Step 1 — CONTROL group on its own terms (executed)

Build: HR **7.17**, baseline CCA carotids, **no band**, terminal beds **pinned**
to the measured CCA Q̄ = 0.016 ml/s and MAP ≈ 92; aortic compliance is the free
knob. Script: `make_control_group.py` → `eberth_control_group.json`.

**Result:**

| aortic C | flow-PI | Q̄ | MAP | pulse pressure |
|---|---:|---:|---:|---:|
| **measured** (2.67e-4) | **2.14** | 0.0160 ✓ | 94 ✓ | **39 ≈ 42.5** ✓ |
| forced to PI 1.16 (6.3e-4, ×2.4) | 1.16 | 0.0160 | 94 | **25** ✗ |
| paper CCA target | 1.16 | 0.016 | ~92 | 42.5 |

**Findings (facts):**
1. **HR 7.17 barely changed the baseline PI** (2.14 vs 2.09 at 6.09 Hz) — heart
   rate is not the missing lever, contrary to the earlier hypothesis.
2. **Beds pinned to the measured flow reproduce Q̄ and MAP exactly** (0.016, 94).
3. **At the *measured* aortic compliance the pressure pulse matches** (PP 39 ≈
   42.5) — but the **flow-PI over-predicts** (2.14 vs 1.16).
4. **Forcing the flow-PI to 1.16 requires ~2.4× the measured aortic compliance,
   which then collapses the pressure pulse to 25** (vs measured 42.5).

**Refined diagnosis:** the discrepancy is *not* HR, bed resistance, the band, or
cross-group parameter sharing. Isolated cleanly in the control group, it is that
**the modeled carotid *flow* pulsatility (2.14) exceeds the measured *velocity* PI
(1.16) at the correct pressure pulse** — and flow-PI and pressure-PP cannot be
matched together (raising aortic compliance to fix one breaks the other). This
points to the carotid bed's *flow admittance* — how flow responds to the pressure
pulse — most plausibly **cerebral autoregulation** (active steadying of flow),
which a passive RCR windkessel cannot reproduce. This is a bed/BC modeling limit,
distinct from the band-amplification limit documented in `calibration_and_limits.md`.

## Status / next

- `eberth_control_group.json` is saved at the **measured** aortic compliance
  (matches Q̄, MAP, and pressure PP; flow-PI reported as an over-prediction) — not
  the compliance-inflated fit, which would break the pressure.
- Step 2 (banded group on its own terms: HR 6.09, remodeled carotids, beds pinned
  to 0.022 / 0.012, band S calibrated) is not yet built. Given the control-group
  finding, the open question there is the same flow-admittance issue plus the
  band-amplification cap.
