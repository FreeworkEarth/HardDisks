# Level 4 v3 — the Mac ladder

2026-09-26. `experiments_energy_transfer/level4_v3_20260926`, **88 runs, 0 aborts, 0 health events**.
Two-compartment box (the pilot's protocol fix). Record lengths set by the B0 prediction written
before the runs.

**Two of the three predictions hold at every mass. The third observable is empty — and it is empty
because the first two are true.**

---

## 1. Predictions, from B0, before the runs

Kinetic theory gives the Brownian-piston friction (derived in the plan, §"How long does the divider
take?")
$$\gamma = 4nh\sqrt{2m\kT/\pi} = 4.118 ,\qquad \tau_v = M/\gamma ,$$
so the *velocity* relaxes linearly in $M$: 2.4, 12.1, 48.6 σ-time at $M_d = 10, 50, 200$. The
*temperature* relaxation does not inherit it — the measured $\tau/\tau_v$ in the spring geometry is
3.7 at $M = 50$ and 13.8 at $M = 200$, itself $\propto M$ — giving $\tau_T \sim M^2/(c\,m\,\gamma)$
with $c \approx 14$, hence **2.4 (saturated) / 45.5 / 729**. Records were set at $\ge 5\tau_T$:
1000 σ-time at $M_d = 10, 50$ (40 seeds), 4000 at $M_d = 200$ (8 seeds).

## 2. What held

| $M_d$ | seeds | record [σ] | ⟨W_in⟩ | far/total | vs ½ | divider settled | vs Δx/2 |
|---|---|---|---|---|---|---|---|
| 10 | 40 | 916 | 7.089 ± 0.066 | 0.488 ± 0.126 | **0.1σ** | −1.867 ± 0.281 | **0.3σ** |
| 50 | 40 | 916 | 7.044 ± 0.075 | 0.438 ± 0.127 | **0.5σ** | −2.245 ± 0.214 | **1.3σ** |
| 200 | 8 | 3916 | 7.282 ± 0.194 | 0.276 ± 0.153 | **1.5σ** | −2.012 ± 0.377 | **0.1σ** |

**(i) The work splits in half** at all three masses, and **(iii) the divider takes half the
compression** at all three. The error bars are set by the per-gas kinetic-energy fluctuation,
$\sqrt{N_s} = 7.1\,\kT$ per seed; that is why the pilot's $2.1\sigma$ at $M_d = 200$ on 20 seeds was
not evidence, and on the longer record it is $1.5\sigma$.

**No mass dependence is detected** in either observable over a factor 20 in $M_d$.

## 3. What failed, and why it had to

The τ_heat fits returned $2.8\times10^9$ and $3.4\times10^8$ at $M_d = 50$ and 200 — i.e. no decay —
and $57.1 \pm 1.3$ at $M_d = 10$ against a predicted 2.4. Looking at the curve rather than the fit:

| $M_d$ | t=0 | 25 | 50 | 100 | 200 | 400 | 800 | late |
|---|---|---|---|---|---|---|---|---|
| 10 | −0.035 | −0.024 | +0.040 | −0.073 | −0.021 | −0.001 | −0.029 | −0.007 |
| 50 | +0.011 | −0.038 | +0.006 | +0.013 | +0.029 | +0.040 | +0.018 | +0.018 |
| 200 | +0.027 | −0.053 | −0.089 | −0.004 | −0.003 | −0.016 | −0.051 | +0.011 |

$T_1 - T_2$ **does not relax, because it was never created.** It wanders in the range
$\pm 0.09\,\kT$ with no systematic decay; the 57.1 at $M_d = 10$ is a fit to that noise.

The reason is geometric and follows from §2. At $u = 0.05$ the push lasts
$\tau_{\mathrm{push}} = 78.6\,\sigma$-time against $\tau_v = 2.4$, 12.1, 48.6 — ratios of 32, 6.5 and
1.6 — so the divider **follows the compression quasi-statically** and takes half of it. Both gases
therefore end at 36.78 and 36.79 σ, i.e. equal lengths, equal compression, equal temperature. Had
the divider been unable to move, gas 2 alone would be compressed 10.1 % and $T_1 - T_2$ would be
$-0.149\,\kT$; what is measured is 20–60 % of that at $M_d = 200$ and essentially nothing at 10,
tracking exactly how well the divider could follow.

**So B0's $\tau_T$ prediction is not refuted by this run — it is untested.** The protocol does not
populate the observable.

## 4. Protocol fix for the ladder

To create a temperature difference the push must outrun the divider, $\tau_{\mathrm{push}} \lesssim
\tau_v$, i.e.

| $M_d$ | 10 | 50 | 200 |
|---|---|---|---|
| $u$ needed | ≳ 1.62 | ≳ 0.32 | ≳ 0.08 |

Those speeds are in the regime Level 2 characterised, where the step protocol launches $Au^2$ of
acoustic energy, so the fast push will inject its own excess; that is a feature here, not a problem,
since the question is how the *difference* relaxes and not how large it starts. The alternative,
cleaner but a new initial condition, is to **start the two gases at different temperatures** and
watch them equalise with no piston motion at all — the exact analogue of the A4 no-push control,
and the measurement that isolates $\tau_T(M)$ from everything else.

**Recommended for KOA:** the no-push ΔT relaxation at $M_d \in \{10, 50, 200, 1000, 5000\}$, records
$\ge 5\tau_T$ from the $M^2$ law, 40 seeds. That is the measurement B0 predicts and this run did not
make.

## 5. Status of the three predictions

| | |
|---|---|
| (i) work splits in half | **holds**, 0.1–1.5σ, no mass dependence over ×20 in $M_d$ |
| (iii) divider takes Δx/2 | **holds**, 0.1–1.3σ |
| (ii) $\tau_T \propto M^2$ | **untested** — the protocol creates no ΔT; §4 says how to test it |

No pass or fail is claimed for Level 4 as a whole; (i) and (iii) are closed and (ii) needs the
protocol of §4.
