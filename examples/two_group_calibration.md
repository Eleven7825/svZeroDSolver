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
to the measured CCA Q̄ = 0.016 ml/s and MAP ≈ 92, aortic compliance at the
**measured** value; the **carotid-bed compliance** is the one calibrated knob.
Script: `make_control_group.py` → `eberth_control_group.json`.

**Result — the passive RCR matches all four CCA targets:**

| knob | flow-PI | Q̄ | MAP | pulse pressure |
|---|---:|---:|---:|---:|
| carotid-bed C = **6.25e-6** (calibrated) | **1.16** ✓ | 0.0161 ✓ | 94 ✓ | **41 ≈ 42.5** ✓ |
| paper CCA target | 1.16 | 0.016 | ~92 | 42.5 |

**Findings (facts):**
1. **HR 7.17 barely changed the baseline PI** (2.14 vs 2.09 at 6.09 Hz) — heart
   rate is not the lever.
2. **Beds pinned to the measured flow reproduce Q̄ and MAP exactly.**
3. **The flow pulsatility is set by the *carotid bed compliance*** (the flow
   admittance): low bed C → flow tracks pressure resistively, PI → PP/MAP ≈ 0.46;
   high bed C → the capacitor adds flow swing, PI → large; an intermediate value
   (6.25e-6) gives PI = 1.16.
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
- Step 2 (banded group on its own terms: HR 6.09, remodeled carotids, beds pinned
  to 0.022 / 0.012, per-carotid bed C + band S) is the next build.
