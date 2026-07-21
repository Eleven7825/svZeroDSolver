# Examples

## `eberth_aortic_band_pulsatile.json`

A pulsatile 0D lumped-parameter model of the **Eberth et al. (2009)** transverse
aortic-arch banding experiment in mice (J. Hypertens. 27:2010–2021). It
implements the electric-analogy circuit in the lab note *"0d pulsatile"*
(`main.pdf`, Fig. 1 and Eq. 6).

### Physiological question

A stiff band on the aortic arch sits **between** the two carotid takeoffs. The
right carotid (RCCA-B) is **upstream** of the band and stays pulsatile; the left
carotid (LCCA-B) is **downstream** and its pulse is damped by the
resistive–inertial band — at nearly constant mean pressure and flow. The band is
the *only* structural difference between the two carotids, so reproducing the
pulsatility asymmetry is the model's key qualitative test.

### Topology

```
heart Qin(t) → aortic valve → ascending aorta (Ra,La,Ca) → node A
   node A ├─→ RCCA-B  → RCR (right head)
          └─→ BAND ΔP_band(Q) → node B
   node B ├─→ LCCA-B  → RCR (left head)
          └─→ descending aorta (Ro,Lo) → RCR (systemic load)
```

### Block mapping (PDF → svZeroDSolver)

| Circuit element (Fig. 1)        | svZeroD block                                             |
|---------------------------------|-----------------------------------------------------------|
| Heart flow source `Qin(t)`      | `FLOW` boundary condition (ejection waveform, one period) |
| Diode aortic valve              | `ValveTanh`                                               |
| Ascending / descending aorta    | `BloodVessel` (`R_poiseuille`, `C`, `L`)                  |
| Band `ΔP_band(Q)` (Young–Tsai)  | `BloodVessel` with `stenosis_coefficient` (+`L`)          |
| Node A / node B                 | `NORMAL_JUNCTION`                                          |
| Terminal beds (Rp–C–Rd)         | `RCR` boundary conditions                                 |

The band uses svZeroD's `BloodVessel` stenosis term, whose governing law
`ΔP = (R + S|Q|)·Q + L·dQ/dt` is exactly the Young–Tsai form in Eq. (6): the
`stenosis_coefficient` **S = Kₜ·ρ/(2A₀²)·(A₀/Aₛ − 1)²** is the Bernoulli/turbulent
term, and `L` is the inertial (`dQ/dt`) term. This nonlinearity + inertance is
what damps the downstream pulse while leaving the upstream node pulsatile.

### Units

Consistent **mmHg / ml / s**: pressure [mmHg], flow [ml/s], time [s],
so R [mmHg·s/ml], C [ml/mmHg], L [mmHg·s²/ml].

### Run it

```bash
uv run svzerodsolver examples/eberth_aortic_band_pulsatile.json out.csv
# or the built binary:  ./svzerodsolver examples/eberth_aortic_band_pulsatile.json out.csv
```

### What it currently produces (illustrative parameters)

The band damps the pressure pulse ~2× across it, MAP stays similar, and mean
carotid flows fall in the paper's range:

| Node                 | MAP (mmHg) | Pulse pressure (mmHg) | Q̄ (ml/s) |
|----------------------|-----------:|----------------------:|----------:|
| A — RCCA-B (upstream)   | ~88 | ~26 | ~0.018 |
| B — LCCA-B (downstream) | ~84 | ~13 | ~0.017 |

The **direction** is correct (upstream more pulsatile than downstream). The
absolute pulsatility indices are **not** yet calibrated to Table 1
(PI_RCCA-B ≈ 3.11, PI_LCCA-B ≈ 1.65).

### Calibrating to the paper (starting point, not a fit)

Parameter values are illustrative mouse-scale numbers. To match the Table-1
targets, tune:

- **`stenosis_coefficient`** and **`L`** of `aortic_band` — set the pulse damping
  (calibrate against the paper's ΔP = 4·v²_jet regression → Kₜ).
- the three **RCR** beds (`Rp`, `C`, `Rd`) — set mean flow split and MAP.
- the **ejection waveform** (`Q`, `t` in the `INFLOW` BC) — heart rate / stroke
  volume (Eberth report ~6.09 Hz banded vs 7.17 Hz baseline).

See `tests/cases/steadyFlow_calibration.json` and the `svZeroDTuner` app
(`applications/svZeroDTuner`) for automated parameter estimation.
