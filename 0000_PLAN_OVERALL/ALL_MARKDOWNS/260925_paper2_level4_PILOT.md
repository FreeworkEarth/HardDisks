# Level 4 pilot — transmission through a free divider

2026-09-25. Geometry B. Two pilots: v1 in the master box (40 runs) and v2 in the two-compartment box
(40 runs), both **0 aborts, 0 health events**. Scripts `level4_pilot_20260925.sh`,
`level4_pilot_v2_20260925.sh`.

**No pass or fail is claimed. The pilot's job was to fix the protocol, and it did — by failing.**

---

## 1. Predictions, stated before the measurement

- **(i)** Quasi-statically the divider moves to pressure balance, so at first order each gas takes
  half the volume change and half the work: **W_far/W_total → 1/2**.
- **(ii)** T₁ − T₂ relaxes through the adiabatic-piston mechanism with τ_heat(M_d). A4 measured
  τ_heat = **48 σ-time at M = 50 and 154 at M = 200** in the spring geometry.
- **(iii)** The divider's settled displacement is **Δx/2 = 1.965 σ**, toward the far gas.

## 2. The pilot's main finding: v1 could not measure its own observable

v1 ran geometry B inside the master box and returned `ΔKE_far = 0.000 ± 0.000` for every seed. The
cause is in the trace schema, not the physics: `00ALLINONE.c:17137` computes

```c
if (sgi == 0) ke_left  += segment_ke[sgi];
else          ke_right += segment_ke[sgi];
```

`KE_gas_left` is **segment 0 only**; `KE_gas_right` is *every other segment*. In the master box
segment 0 is the empty spring compartment, so the two columns read "zero" and "both gases combined",
and the work split across the divider — the entire point of Level 4 — is invisible.

**Protocol fix, no code change required:** drop the empty spring compartment for Level 4 and use the
two-compartment box, where segment 0 *is* gas 1 and segment 1 *is* gas 2:

```
--particles=100 --particles-boxes=50,50 --l0=39.25 --num-walls=1 --wall-positions=39.25
```

box 0 … 78.5, divider centre 39.25, two compartments of 38.75 σ with 50 disks each, η = 0.1013 —
the same gas as every other level. Level 0b already established that an empty compartment behind a
held wall changes nothing the gas can see, so this costs no continuity with Levels 0–3. Grid-exact:
39.25 × 24 = 942, 78.5 × 24 = 1884.

## 3. v2, in the corrected geometry

| M_d | seeds | ⟨W_in⟩ | ΔKE far gas | ΔKE pushed gas | far/total | divider settled [σ] |
|---|---|---|---|---|---|---|
| 50 | 20 | 7.009 ± 0.107 | 3.202 ± 0.890 | 3.175 ± 0.865 | **0.502 ± 0.140** | **−2.001 ± 0.184** |
| 200 | 20 | 7.269 ± 0.120 | 4.735 ± 0.607 | 2.225 ± 0.588 | 0.680 ± 0.087 | −2.133 ± 0.110 |

**(i) — confirmed at M_d = 50, open at 200.** far/total = 0.502 ± 0.140 is **0.0σ from 1/2**. At
M_d = 200 it is 0.680 ± 0.087, **2.1σ high**. The heavy divider's record is 917 σ-time against
τ_heat ≈ 154, which should be ample, so this is not obviously a record-length artefact and is the
first thing the KOA ladder should resolve.

**(iii) — confirmed at both masses.** |−2.001| against 1.965 is **0.2σ**; |−2.133| is 1.5σ. The
divider takes half the compression, as pressure balance requires.

**(ii) — NOT measured. The estimator is broken.** A 1/e crossing of the ensemble-mean T₁ − T₂
returned τ = 0 and 6 σ-time against A4's 48 and 154 — the crossing fires immediately, which means
the statistic is not measuring a relaxation. This needs a proper exponential fit with a stated
baseline and a stated start time (the push ends at 79 σ-time), not a threshold crossing, and it
should be validated against A4's own curves before being trusted. **No τ_heat value from this pilot
should be quoted.**

## 4. What the KOA ladder needs, revised

- Two-compartment box, as §2.
- M_d ∈ {10, 50, 200, 1000, 5000}, 40 seeds, records ≥ 5 τ_heat(M_d).
- A **validated** τ estimator: fit A(t) = A_∞ + (A₀ − A_∞)e^{−t/τ} to the seed-averaged T₁ − T₂ from
  the end of the push, and reproduce A4's 48/154 on the A4 runs as the acceptance test before
  applying it here.
- Resolve the M_d = 200 work split: is 0.68 a real mass dependence, or an equilibration that the
  record has not completed?

## 5. Pictures

Not taken. The two `--render` shots were not run for this pilot because the geometry changed
mid-pilot; they belong with the corrected two-compartment geometry and are carried to the KOA
campaign setup.
