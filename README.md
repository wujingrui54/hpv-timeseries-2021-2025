# HPV high-risk genotype time-series — Xi'an, China, 2021–2025 (analysis code & reproducibility archive)

Analysis code, de-identified aggregate data, and output figures and tables underlying the manuscript *"Shifting oncogenic HPV ecology after COVID-19: emergence of HPV-51 and decoupling of HPV-16 in a cervical-cancer screening cohort in Xi'an, China (2021–2025)"* (Wu et al, 2026; manuscript in preparation).

> **Persistent identifier.** This GitHub repository is the working home of the analysis code and aggregate outputs. A permanent, versioned **DOI will be minted via Zenodo from a tagged GitHub release upon journal acceptance**, and added here. (The `.zenodo.json` in this repository pre-configures that deposit.)

- **Cohort**: 18 185 women, Norinco General Hospital (Xi'an, Shaanxi, China), 5 Jan 2021 – 31 Dec 2025
- **IRB**: Norinco General Hospital, Approval No. 202601271338000343089 (issued 27 Jan 2026) for a low-risk retrospective analysis of de-identified clinical records under project 26z03; the de-identified HPV record stream analysed spans 5 Jan 2021 – 31 Dec 2025.
- **Funding**: Xi'an Medical College Scientific Research Project, Grant No. **26z03** (PI 吴静蕊 / Jingrui Wu).

## Authors (8 ICMJE-qualified)

Jingrui Wu¹,⁷; Haihui Yang¹; Jinhua Zhi¹; Jinjin Yan²; Wen Wang³; Ki-Hyun Kim⁴; Wenping Yang⁵; Le Fang⁶,\*

¹ College of Medical Technology, Xi'an Medical College, Xi'an, Shaanxi, China
² Department of Basic Sciences, Xi'an Medical College, Xi'an, Shaanxi, China
³ Clinical Medical College, Xi'an Medical College, Xi'an, Shaanxi, China
⁴ Department of Physical Therapy, Kyungwoon University, Gumi, Gyeongsangbuk-do, Republic of Korea
⁵ Department of Nursing, Kyungwoon University, Gumi, Gyeongsangbuk-do, Republic of Korea
⁶ Department of Clinical Laboratory, Norinco General Hospital, Xi'an, Shaanxi, China
⁷ Department of Public Health, Kyungwoon University, Gumi, Gyeongsangbuk-do, Republic of Korea

\* Co-corresponding authors: Le Fang (289278162@qq.com) and Jingrui Wu (253207006@ikw.ac.kr)

## Get the code

```bash
git clone https://github.com/wujingrui54/hpv-timeseries-2021-2025.git
```

## Contents

```
hpv-timeseries-2021-2025/
├── README.md                                  this file
├── CHANGELOG.md                               version history
├── MANIFEST.txt                               file listing with SHA-256 checksums
├── LICENSE                                    licensing terms (see License section)
├── .zenodo.json                               metadata for the Zenodo DOI to be minted at acceptance
├── code/
│   ├── hpv_timeseries_complete.ipynb          canonical analysis pipeline (Jupyter)
│   ├── methodology_upgrade.py                 reproduces Fig S6 + Tables S9–S11 (WHO std / harmonic GLM / monthly joinpoint)
│   └── qc_verification_18185.py               independent publication-grade QC (10 layers, 412+ assertions)
├── data/
│   └── canonical_aapc_table.json              machine-readable canonical results (seed=20260520, B=1500, N=18185)
├── manuscript/
│   ├── manuscript.docx                        main manuscript snapshot
│   └── supplementary_appendix.docx            supplement (3 Methods Notes; Supplementary Tables S0–S10; Supplementary Figures S1–S6)
└── figures/
    └── 9 figures (PDF + PNG); FigS6_methodology_upgrade.* is the methodology-upgrade figure, the rest are the main + extended-data figures
```

> **Figure-naming note.** The deposited figure files retain their analysis-pipeline output names; the manuscript renumbers them. Mapping: Fig1_main_4panel = **Figure 1** (temporal dynamics); Fig2_genotype_landscape = **Figure 2** (genotype landscape); Fig5_age_3panel = **Figure 3** (age structure); EDF1_STROBE_flow_18185 = **Supplementary Figure S1**; EDF3_joinpoint_4panel = **Supplementary Figure S2**; EDF2_forest_18metric = **Supplementary Figure S3**; Fig3_monthly_ITS = **Supplementary Figure S4**; Fig4_SARIMA_forecast = **Supplementary Figure S5**; FigS6_methodology_upgrade = **Supplementary Figure S6**. The manuscript presents 3 main figures + 3 main tables, with the remaining display items as Supplementary Figures S1–S6 and Supplementary Tables S0–S10.

## Reproducing every number in the manuscript

Independent QC verification (no Jupyter required, ~5 s runtime):

```bash
# Requires: pandas, numpy, scipy, openpyxl
python code/qc_verification_18185.py \
    --raw "/path/to/2021_2025HPV _12_31.xlsx" \
    --canonical data/canonical_aapc_table.json
```

Expected output: `OVERALL: ALL QC LAYERS PASSED` (exit code 0). The script checks **10 QC layers** end-to-end (412 independent assertions; Wilson CI worst diff < 1.78 × 10⁻¹⁵, the floating-point limit).

## Full pipeline reproduction

```bash
jupyter nbconvert --to notebook --execute code/hpv_timeseries_complete.ipynb \
    --output executed_pipeline.ipynb
```

## Data-quality handling (disclosed)

The raw 18 670-record export contains three records with the laboratory free-text label `弱Positive` (one each for HPV-52, HPV-66, HPV-81). Per the manufacturer's reporting convention these are treated as positive throughout; a sensitivity exclusion changes no AAPC point estimate by more than 0.01 percentage points. See `manuscript/supplementary_appendix.docx` → Supplementary Methods Note 1.

## Bootstrap implementation (disclosed)

All 95% confidence intervals were obtained from 1 500 **binomial bootstrap** resamples — for each genotype-year the positive count was redrawn from a binomial distribution defined by that year's woman-level denominator and observed prevalence — using NumPy's `default_rng` PRNG (MT19937 Mersenne Twister) seeded at **20260520**. See `manuscript/supplementary_appendix.docx` → Supplementary Methods Note 2.

## Co-infection (HR-HPV multiplicity) structure

Among the 18 185 women in the analytic cohort:

| HR types carried | Women | Cohort % |
|---|---|---|
| 0 | 15 867 | 87.25% |
| 1 | 1 879 | 10.33% |
| 2 | 345 | 1.90% |
| ≥3 | 94 | 0.52% |

## Not included in this repository

To comply with patient privacy and IRB conditions, the individual-level raw genotyping export (`2021_2025HPV _12_31.xlsx`, 18 670 rows) is **not** included here. It is available from the corresponding authors on reasonable request, conditional on Norinco General Hospital data-sharing approval.

## Relationship to the prior Wu et al BMC Infectious Diseases 2026 paper

A subset of the present cohort (2023–2025; n = 11 413, 62.8% of the present 18 185) overlaps with our published cross-sectional study:

> Wu J, Yi J, Fang L, *et al*. Urban–rural disparities in HPV prevalence and genotype distribution: a large-scale cervical cancer screening study in Xi'an, China. *BMC Infect Dis* 2026; **26**: 13647. doi:10.1186/s12879-026-13647-2

The two studies differ in research question (5-year temporal dynamics vs single-cross-section urban-versus-rural disparity), design (single-centre longitudinal vs cross-sectional), statistical core (AAPC + interrupted time-series + SARIMA vs age-standardised prevalence ratios), and eligibility (age ≥21 vs 15–100 years). No result, figure, or table of the prior paper is reproduced here.

## Citation

Wu J, Yang H, Zhi J, *et al*. *Shifting oncogenic HPV ecology after COVID-19: emergence of HPV-51 and decoupling of HPV-16 in a cervical-cancer screening cohort in Xi'an, China (2021–2025)*. Manuscript in preparation, 2026.

Reproducibility code & data: this repository — https://github.com/wujingrui54/hpv-timeseries-2021-2025. A citable, versioned DOI will be minted via Zenodo from a tagged release upon journal acceptance.

## License

- **Code** (`code/`): MIT License
- **Manuscript files** (`manuscript/`, `figures/`): CC BY 4.0 (author copyright retained until publication)
- **Canonical JSON** (`data/`): CC0 (public-domain dedication)

See `LICENSE` for details.

## Contact

The manuscript has two designated **co-corresponding authors**:

- **Jingrui Wu** (co-corresponding author; also first and submitting author; pre-acceptance correspondence): College of Medical Technology, Xi'an Medical College, No. 1 Weiwu Road, Huyi District, Xi'an, Shaanxi 710309, China; and Department of Public Health, Kyungwoon University, Gumi 39160, Republic of Korea. E-mail: 253207006@ikw.ac.kr. ORCID: 0009-0006-0123-808X.
- **Le Fang, MSc** (co-corresponding author; post-acceptance correspondence): Department of Clinical Laboratory, Norinco General Hospital, Xi'an, Shaanxi 710043, China. E-mail: 289278162@qq.com. ORCID: 0009-0009-3486-6076.
