# Level 3 v5 — the peak was a statistic, not the gas

2026-09-23. **Post-processing only, no new runs.** Analysis on `level3_master_preload_20260921` and
`level3_v4_20260922`. Supersedes v4 §5's interpretation.

**v4's conclusion "the residual is the acoustic pulse" is WITHDRAWN.** Two tests kill it, and a
third shows what the excess actually was. On the settled state the model agrees to 2–4 %.

---

## 1. What "peak s" was, and what the control does with the same statistic

v4's "peak s" was statistic **(b): the maximum over time of the ensemble mean ⟨s(t)⟩**, taken over
the window from the end of the push to the end of the record, with **no bias correction** — the
same `argmax` bias I had measured and subtracted for E_coh, left uncorrected here.

Both statistics, on the push runs and on the A4 no-push control over a **time-matched window**:

| cell | (a) mean of per-seed max | (b) max of ⟨s⟩ | settled s̄ | (a)−s̄ | **(b)−s̄** | control (a)−mean | **control (b)−mean** |
|---|---|---|---|---|---|---|---|
| u=0.02, M=50 | 4.249 | 1.367 | 0.904 | +3.345 | **+0.463** | +3.026 | **+0.446** |
| u=0.02, M=200 | 3.079 | 1.241 | 0.901 | +2.178 | **+0.339** | +1.948 | **+0.354** |
| u=0.05, M=50 | 3.735 | 1.434 | 0.911 | +2.825 | **+0.524** | +2.323 | **+0.445** |
| u=0.05, M=200 | 2.319 | 1.214 | 0.922 | +1.397 | **+0.292** | +1.318 | **+0.126** |
| u=0.1, M=50 | 3.996 | 1.507 | 0.889 | +3.107 | **+0.618** | +2.496 | **+0.436** |
| u=0.1, M=200 | 2.741 | 1.697 | 1.008 | +1.734 | **+0.690** | +1.441 | **+0.182** |

**At u = 0.02 the control reproduces the entire excess** (+0.463 vs +0.446; +0.339 vs +0.354). The
control has no piston motion at all, so whatever produces that number is not the gas responding to
a push. The predicted sizes hold too: thermal rms √(kT/k_eff) = 1.462 σ with
k_eff = k_gas Δx/s_qs = 0.4677, giving 0.231 σ on a 40-seed mean; statistic (a) should sit 1–2 rms
above the mean and measures +1.3 to +3.3 σ; statistic (b) should sit a few tenths above and measures
+0.13 to +0.45 σ. Both as predicted.

A small genuine remainder survives at the faster end: (b) excess minus control excess is
+0.017 / −0.015 at u = 0.02, +0.079 / +0.166 at u = 0.05, +0.182 / +0.508 at u = 0.1. It grows with
u, but see §2 for how it compares with the acoustic estimate.

## 2. Acoustic scaling test — the excess does not scale like a pulse

A pulse launched by the piston start carries pressure ρ c_s u, i.e. a wall force
h ρ c_s u = (N m/L) c_s u, and can push the wall at most 2×(that)/k_eff beyond settled:

| u | pulse force | acoustic max overshoot | measured (b)−settled | control-subtracted |
|---|---|---|---|---|
| 0.02 | 0.0444 | **0.190 σ** | +0.401 | ≈ 0.00 |
| 0.05 | 0.1111 | **0.475 σ** | +0.408 | +0.12 |
| 0.1 | 0.2222 | **0.950 σ** | +0.654 | +0.35 |

Acoustic predicts a **factor 5 across the scan**; the raw measurement is flat (0.40, 0.41, 0.65),
and the control-subtracted remainder is consistent with zero at u = 0.02 where the acoustic estimate
is 0.19 σ. v4's own A3 said the same thing and I did not hear it: the ratio was flat from 1.8 to 8.8
sound traversals. **A transient that does not scale with u is not the acoustic pulse.**

## 3. ε on the settled state

ε_settled = [F s̄ + ½k s̄²]/W_qs, with s̄ the per-seed time average over t > τ_push + 3 T_w, and the
thermal contribution removed by the A4 control over a time-matched window with the same statistic.
W_qs = N(T_f/T_i − 1) = 100(1.1432 − 1) = **14.322 kT**, so the model's ε_qs = 1.500/14.322 =
**0.1047**.

| u | M_s | s̄ push | s̄ control | s̄ corrected | model s_qs | ΔE_settled | model | ε_settled | ε model | σ |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 50 | 0.911 ± 0.006 | 0.040 ± 0.005 | 0.871 ± 0.008 | 0.840 | 1.561 ± 0.016 | 1.500 | 0.1090 | 0.1047 | **3.9** |
| 0.02 | 200 | 0.901 ± 0.011 | 0.043 ± 0.007 | 0.858 ± 0.013 | 0.840 | 1.536 ± 0.026 | 1.500 | 0.1073 | 0.1047 | 1.4 |
| 0.05 | 50 | 0.904 ± 0.024 | 0.029 ± 0.013 | 0.875 ± 0.027 | 0.840 | 1.570 ± 0.055 | 1.500 | 0.1096 | 0.1047 | 1.3 |
| 0.05 | 200 | — record shorter than τ + 3 T_w = 519 σ-time; window empty | | | | | | | | |
| 0.1 | 50 | 0.888 ± 0.019 | 0.041 ± 0.013 | 0.847 ± 0.023 | 0.840 | 1.514 ± 0.046 | 1.500 | 0.1057 | 0.1047 | 0.3 |
| 0.1 | 200 | 1.026 ± 0.093 | 0.008 ± 0.037 | 1.018 ± 0.100 | 0.840 | 1.862 ± 0.209 | 1.500 | 0.1300 | 0.1047 | 1.7 |

**Level 3 does NOT pass the stated ≤1σ criterion:** 0.3σ, 1.3σ, 1.4σ, 1.7σ and one at 3.9σ. But the
disagreement has collapsed from 60 % in the peak to **2–4 % in s̄**
(s̄/s_qs = 1.037, 1.021, 1.042, 1.008, 1.212), and the 3.9σ cell is the one with the tightest error
bar (0.9 %), not the largest discrepancy.

**A candidate for the 2–4 %, stated not claimed:** the A1 held-wall run measured the standing force
as **1.5961 before any signal can arrive, against the EOS value F = N kT Z/L = 1.5749 — +1.3 %**.
A standing force 1.3 % above the value the model uses produces a settled displacement roughly 1.3 %
high, which is a third to a half of what is seen. The rest is not accounted for.

**Assumptions stated:** (i) the window t > τ + 3 T_w is where the wall has completed three free
periods; the M_s = 200, u = 0.05 record does not contain one and that cell is reported as missing
rather than computed on a shorter window; (ii) the control is subtracted as a per-seed time average
over a window of the same *length*, not the same absolute times, because the control has no push;
(iii) W_qs is the adiabat's N(T_f/T_i − 1) for the full Δx, not the measured ⟨W_in⟩, per the
criterion as written.

## 4. Radiation reaction — deliberately not built

It adds damping, which *reduces* overshoot, and cannot close an under-prediction. Recorded as a
Level 4 item instead: wall–gas coupling Z_g = N m c_s/L = 2.22, Q ≈ ω_w M/Z_g.

## 5. Recommendation for the plan

Drop the peak from the pass criteria and keep it as a diagnostic; state the criterion on
ε_settled. On that criterion Level 3 is one cell away, and the open question is a 2–4 % systematic
in the settled displacement rather than a factor 1.6 in a transient.

**One sentence for Susanne:** the parameter-free model reproduces where the spring-loaded wall
settles to 2–4 %; what remains is a small systematic in that settled position, and the large
"transient excess" reported earlier was a peak-picking statistic, reproduced by a control run with
no push at all.
