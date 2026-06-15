# Changelog

Notable changes to the **analysis code, data, and figures** in this repository.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning
follows [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html). The primary numerical
results (yearly prevalence, AAPC, prevalence ratios) and the canonical results JSON
(`data/canonical_aapc_table.json`) are **byte-identical across every version** —
later versions add analyses, independent checks, and presentation refinements only.

## v1.0.5
- Added `code/methodology_upgrade.py` — deterministic (seed 20260520) script that reproduces external age-standardisation to the WHO World Standard Population and a harmonic seasonal Poisson / negative-binomial GLM trend per genotype; writes Supplementary Figure S6 and Supplementary Tables S9–S10, with a self-check against the canonical yearly prevalence.
- Added `figures/FigS6_methodology_upgrade` (2-panel: WHO standardisation; genotype two-method forest).
- Figure number formatting harmonised (thin-space thousands separators, period decimals) on the STROBE-flow, forest, and age figures — **formatting only; no data values changed.**
- Manuscript and supplement snapshots refreshed. Primary results unchanged; canonical JSON byte-identical.

## v1.0.4
- Documentation and metadata refresh; manuscript and supplement snapshots updated. No code, data, figures, or numerical results changed.

## v1.0.3
- Self-consistency documentation fix. No code, data, figures, or numerical results changed.

## v1.0.2
- Metadata refresh. No code, data, figures, or numerical results changed.

## v1.0.1
- Added `code/qc_verification_18185.py` — a stand-alone, Jupyter-free verifier that recomputes the full result set from the raw export across 10 independent layers (412 assertions; exits 0 on success; ~5 s runtime).
- Added Supplementary Methods Notes documenting: data-quality handling (three weakly-positive free-text labels — a sensitivity exclusion changes no AAPC by > 0.01 percentage points); bootstrap implementation (woman-level resampling; NumPy `default_rng` MT19937; seed 20260520; two-sided percentile p values); and the HR-HPV multiplicity (co-infection) structure (0 / 1 / 2 / ≥3 HR types = 15 867 / 1 879 / 345 / 94 women).
- Figure-presentation refinements (axis padding, log-axis tick labels, STROBE-flow layout). No numerical result changed.

## v1.0.0
- Initial release: canonical analysis pipeline (`code/hpv_timeseries_complete.ipynb`), canonical results JSON (`data/canonical_aapc_table.json`; 1500-replicate bootstrap, seed 20260520, N = 18 185), main and supplementary figures (PDF + PNG), and the manuscript + supplementary snapshots.
- Cohort: 18 185 women, Norinco General Hospital, Xi'an, 5 Jan 2021 – 31 Dec 2025.
