# 260909 — TASK C: fit-input manifest and refit of route-A c_s(η) on eligible trajectories (Claude Code)

Inputs: `campaign_r25_psi6_20260823/analysis/analysis/eta_*/speed_of_sound_runs.csv` (per-run ν and FFT quality), the leaf `run.log` (run → seed, health counters), `speed_of_sound_psi6.csv` (seed), and `analysis_paper1_20260908/trajectories.csv` (GPT structural flags, by seed). Script: `hspist3/validation/refit_sound_speed_manifest.py`. Outputs in `campaign_r25_psi6_20260823/analysis/refit_20260909/` and copied to `260909_plots/` (`routeA_fit_input_manifest_20260909.csv`, `routeA_refit_cs_vs_eta_20260909.csv`, `routeA_unmatched_joins_20260909.csv`). Originals untouched.

## Manifest

- rows: **7200** (32 η × 9 masses × 25 repeats); joins matched for all rows (unmatched: 0).
- eligible for the fit: **7186**; excluded **14**, all for `peak_on_search_boundary` (which is also the only `frequency_quality` failure): η=0.5500: 5, η=0.7600: 9.
- health: forced_advance = clamp_repair = overlap_repair = 0 in all 7200; `wall_overdue` ≠ 0 in 504, all in the documented t = 0 seeder class (see `260909_wall_overdue_and_temperature_resolution_CC.md`), retained.
- ψ₆ at dilute η is undefined for some trajectories (no neighbours inside 1.4 σ); that is not an exclusion for a sound-speed fit and is left blank.

## Fit

Román relation ν = c_s·K/(2π L_eff), K the fundamental root of cot K = αK by bisection (`k_root_bisect`, α = M/(2·N_side), N_side = 50), L_eff = L0 − 2r, per-mass mean ν over eligible repeats weighted by 1/sem². Two models are reported: **through the origin** (the relation has no intercept) and **free intercept** (the model the 2026-08-23 summary used; `force_zero_intercept` defaults to off in `plot_speed_of_sound_edmd.py`). EOS: Kolafa–Rottner 2006 (quoted for η ≤ 0.69 only) and Liu 2021 global, both through c_s² = (kT/m)[Z + ηZ′ + Z²] with kT/m = 1 (T = 1 exactly at t = 0 by construction, A.2).

| η | eligible/all | c_s (origin) | ± | c_s (free icpt) | intercept | 08-23 summary (free) | KR 2006 | dev vs KR | Liu 2021 | dev vs Liu |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0196 | 225/225 | 1.4641 | 0.0015 | 1.4646 | -0.00000 | 1.4582 | 1.47138 | -0.49% | 1.47138 | -0.49% |
| 0.0262 | 225/225 | 1.4871 | 0.0025 | 1.4956 | -0.00000 | 1.5039 | 1.49118 | -0.28% | 1.49118 | -0.28% |
| 0.0393 | 225/225 | 1.5288 | 0.0019 | 1.5329 | -0.00000 | 1.5282 | 1.53196 | -0.21% | 1.53196 | -0.21% |
| 0.0524 | 225/225 | 1.5741 | 0.0021 | 1.5684 | +0.00000 | 1.5621 | 1.57439 | -0.02% | 1.57439 | -0.02% |
| 0.0785 | 225/225 | 1.6622 | 0.0013 | 1.6680 | -0.00001 | 1.6658 | 1.66454 | -0.14% | 1.66455 | -0.14% |
| 0.1122 | 225/225 | 1.8012 | 0.0021 | 1.7934 | +0.00001 | 1.7923 | 1.79189 | +0.52% | 1.79194 | +0.52% |
| 0.1309 | 225/225 | 1.8781 | 0.0019 | 1.8797 | -0.00000 | 1.8831 | 1.86885 | +0.49% | 1.86893 | +0.49% |
| 0.1571 | 225/225 | 2.0024 | 0.0020 | 2.0029 | -0.00000 | 2.0007 | 1.98495 | +0.88% | 1.98509 | +0.87% |
| 0.1963 | 225/225 | 2.1992 | 0.0020 | 2.2036 | -0.00001 | 2.2027 | 2.17980 | +0.89% | 2.18006 | +0.88% |
| 0.2618 | 225/225 | 2.5985 | 0.0033 | 2.6178 | -0.00008 | 2.6169 | 2.57255 | +1.01% | 2.57300 | +0.99% |
| 0.3927 | 225/225 | 3.7877 | 0.0033 | 3.7797 | +0.00005 | 3.7753 | 3.74608 | +1.11% | 3.74615 | +1.11% |
| 0.5236 | 225/225 | 6.0790 | 0.0146 | 6.1111 | -0.00026 | 6.1282 | 5.93212 | +2.48% | 5.92901 | +2.53% |
| 0.5500 | 220/225 | 6.7313 | 0.0168 | 6.7139 | +0.00013 | 6.7824 | 6.60184 | +1.96% | 6.59786 | +2.02% |
| 0.5700 | 225/225 | 7.2379 | 0.0157 | 7.3273 | -0.00080 | 7.3752 | 7.18768 | +0.70% | 7.18345 | +0.76% |
| 0.5900 | 225/225 | 7.8528 | 0.0330 | 7.9897 | -0.00127 | 7.9349 | 7.85473 | -0.02% | 7.85141 | +0.02% |
| 0.6100 | 225/225 | 8.4916 | 0.0517 | 8.7815 | -0.00298 | 8.8065 | 8.61562 | -1.44% | 8.61587 | -1.44% |
| 0.6300 | 225/225 | 9.4125 | 0.0408 | 9.5814 | -0.00201 | 9.5809 | 9.48075 | -0.72% | 9.48837 | -0.80% |
| 0.6500 | 225/225 | 10.9525 | 0.0373 | 11.0717 | -0.00133 | 11.1478 | 10.43888 | +4.92% | 10.44997 | +4.81% |
| 0.6700 | 225/225 | 13.0195 | 0.0226 | 13.1473 | -0.00150 | 13.1449 | 11.37327 | +14.47% | 11.34499 | +14.76% |
| 0.6800 | 225/225 | 14.4885 | 0.0334 | 14.6077 | -0.00142 | 14.6666 | 11.67118 | +24.14% | 11.60169 | +24.88% |
| 0.6900 | 225/225 | 16.5303 | 0.0428 | 16.4917 | +0.00042 | 16.5371 | 11.55103 | +43.11% | 11.51687 | +43.53% |
| 0.6950 | 225/225 | 17.7803 | 0.0196 | 17.7963 | -0.00017 | 17.7167 | — | — | 11.26418 | +57.85% |
| 0.7000 | 225/225 | 19.3852 | 0.0262 | 19.3880 | -0.00003 | 19.3011 | — | — | 10.83136 | +78.97% |
| 0.7050 | 225/225 | 17.9750 | 0.2185 | 17.1757 | +0.01201 | 17.1757 | — | — | 10.23243 | +75.67% |
| 0.7100 | 225/225 | 16.2777 | 0.1648 | 15.9851 | +0.00428 | 15.9851 | — | — | 9.64650 | +68.74% |
| 0.7150 | 225/225 | 15.7272 | 0.0612 | 16.1732 | -0.00426 | 16.1142 | — | — | 9.71093 | +61.95% |
| 0.7200 | 225/225 | 16.5383 | 0.0256 | 16.6004 | -0.00074 | 16.5567 | — | — | 11.83387 | +39.75% |
| 0.7250 | 225/225 | 17.9122 | 0.0231 | 17.7752 | +0.00179 | 17.7635 | — | — | 12.16557 | +47.24% |
| 0.7300 | 225/225 | 18.5565 | 0.0096 | 18.5473 | +0.00012 | 18.5209 | — | — | 12.51138 | +48.32% |
| 0.7400 | 225/225 | 22.7531 | 0.0158 | 22.6997 | +0.00073 | 22.6958 | — | — | 13.26577 | +71.52% |
| 0.7500 | 225/225 | 27.8938 | 0.0365 | 27.8549 | +0.00055 | 27.8102 | — | — | 14.11561 | +97.61% |
| 0.7600 | 216/225 | 36.8617 | 0.5858 | 35.1740 | +0.02354 | 37.4530 | — | — | 15.08147 | +144.42% |

## Reading

- **Why the refit differs from the 2026-08-23 summary.** The old summary is a free-intercept fit; the refit's free-intercept column reproduces it to < 0.7 % everywhere except η = 0.55 (−1.0 %) and 0.76 (−6.1 %), where the refit drops boundary-flagged peaks the old fit kept. The remaining sub-percent scatter is the per-mass weighting. So nothing in the ν data changed; the model and the eligibility did.
- **Which value to quote.** The Román relation passes through the origin; the through-origin c_s is the prescribed estimator. Where origin and free-intercept values differ by more than their errors (η = 0.57, 0.59, 0.61, 0.65, 0.705–0.715, 0.76) the fitted intercept is absorbing a systematic — the ν(K) points do not extrapolate to zero, i.e. L_eff or the single-mode relation is off at that η. That is a diagnostic, reported here, not something to remove by choosing the model that agrees better with the EOS.
- **Dilute anchor (η ≤ 0.08).** Through-origin c_s is within −0.49 … −0.02 % of both EOS (they coincide there), with per-η errors 0.0013–0.0025. With the 504 flags resolved, the ideal-gas end of the c_s data stands.
- **η = 0.11–0.39:** +0.5 … +1.1 % above both EOS, consistently positive — larger than the fit errors (0.002–0.003). Not re-interpreted here; the same sign appeared in the 08-23 summary.
- **η ≥ 0.65:** the excess over the EOS is the known finite-size/structure effect (fixed-aspect campaign 2026-08-26); route-A N = 100 values above ≈ 0.65 are not bulk values and are plotted as such.
- Not done here (needs new runs or a different dataset): the fixed-aspect family-A finite-size fit with model sensitivity, per-mass residual tables, and the heavy-divider subset — the manifest carries the per-mass inputs for them.
