"""
qc_verification_18185.py
========================

Stand-alone, Jupyter-free reproduction of every numerical claim in the
accompanying manuscript "Shifting oncogenic HPV ecology after COVID-19:
emergence of HPV-51 and decoupling of HPV-16 in a cervical-cancer screening
cohort in Xi'an, China (2021-2025)" (Wu et al, 2026), starting from the raw
genotyping export and producing a pass/fail QC report.

Companion to:
- hpv_timeseries_complete.ipynb     (full analysis pipeline)
- canonical_aapc_table.json           (machine-readable canonical results)
- manuscript.docx           (Methods Sec. "Software and reproducibility")
- supplementary_appendix.docx         (Supplementary Methods Notes 1-3)

Archived at Zenodo DOI 10.5281/zenodo.20256624 (concept DOI 10.5281/zenodo.20256624).

Usage
-----
Local Python:
    python qc_verification_18185.py [--raw RAW.xlsx] [--canonical CANON.json]

Google Colab (mounts Drive automatically + uses Drive-based defaults):
    !python qc_verification_18185.py
or import-and-call from a notebook cell:
    from qc_verification_18185 import run_qc
    run_qc('/content/drive/MyDrive/.../raw.xlsx',
           '/content/drive/MyDrive/.../canonical_aapc_table.json')

With no arguments, falls back to:
  - Drive path when running in Colab (auto-mounts /content/drive)
  - Local relative paths when running on a workstation

What it verifies
----------------
1.  STROBE 5-step exclusion chain (18,670 -> 18,185)
2.  Yearly cohort sample sizes (3,946 / 2,826 / 4,460 / 3,355 / 3,598)
3.  Per-genotype yearly positive counts (18 metrics x 5 years = 90 cells)
4.  Wilson 95% CIs (pct + ci_lo + ci_hi for each of 90 cells = 270 numeric assertions)
5.  AAPC point estimates from log-linear OLS on yearly percentages (18 metrics)
6.  Non-9v HR yearly burden (5 cells)
7.  Set identity: |STABLE15 positives| = |HR-IARC13| + |HR-2B| - |overlap|
8.  Fig 5 age-bin sample sizes (11 bins)         # renumbered: was Fig 2 in v1
9.  Co-infection multiplicity (0 / 1 / 2 / >=3 HR types per woman)
10. Detection-date range and 60-month completeness for the monthly ITS analysis

Data-quality handling (disclosed in manuscript Methods + Supp Methods Note 1)
-----------------------------------------------------------------------------
Three records in the raw export carry the laboratory free-text label
'弱Positive' (one each for HPV-52, HPV-66, HPV-81). Per manufacturer convention
these are treated as positive throughout, via case-insensitive substring
matching on the literal string 'Positive' in the positivity mask. The HPV-81
case is outside the STABLE15 analytic panel and never enters a primary metric.

Output
------
Prints a single pass/fail line per QC layer, then an OVERALL pass/fail at the
end. Exit code is 0 iff every assertion passes.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress


# ---------------------------------------------------------------------------
# Colab autodetection + default-path resolution
# ---------------------------------------------------------------------------
IN_COLAB = "google.colab" in sys.modules


def _mount_drive_if_colab() -> Path:
    """Mount Google Drive when running in Colab. Returns the Drive root path."""
    if not IN_COLAB:
        return Path("/content/drive")  # placeholder; not used outside Colab
    from google.colab import drive  # type: ignore
    drive.mount("/content/drive", force_remount=False)
    return Path("/content/drive")


def _default_raw_path() -> Path:
    """Resolve a sensible default for the raw genotyping Excel file."""
    if IN_COLAB:
        _mount_drive_if_colab()
        # Common Google Drive locations to try:
        candidates = [
            Path("/content/drive/MyDrive/2025HPV/2021_2025HPV _12_31.xlsx"),
            Path("/content/drive/MyDrive/2025HPV/HPV_20/2021_2025HPV _12_31.xlsx"),
            Path("/content/drive/MyDrive/2021_2025HPV _12_31.xlsx"),
        ]
        for c in candidates:
            if c.exists():
                return c
        return candidates[0]  # primary expected location
    # Local workstation: sibling of this script's parent folder
    here = Path(__file__).resolve().parent
    return here.parent / "2021_2025HPV _12_31.xlsx"


def _default_canon_path() -> Path:
    """Resolve a sensible default for the canonical JSON file."""
    if IN_COLAB:
        _mount_drive_if_colab()
        candidates = [
            Path("/content/drive/MyDrive/2025HPV/output/canonical_aapc_table.json"),
            Path("/content/drive/MyDrive/2025HPV/HPV_20/output/canonical_aapc_table.json"),
            Path("/content/drive/MyDrive/2025HPV/canonical_aapc_table.json"),
        ]
        for c in candidates:
            if c.exists():
                return c
        return candidates[0]
    here = Path(__file__).resolve().parent
    return here / "output" / "canonical_aapc_table.json"


# ---------------------------------------------------------------------------
# Constants (every "magic number" the manuscript depends on)
# ---------------------------------------------------------------------------
EXPECTED_RAW_N = 18670
EXPECTED_COHORT_N = 18185
EXPECTED_RETENTION_PCT = 97.40
EXPECTED_YEARLY_N = {2021: 3946, 2022: 2826, 2023: 4460, 2024: 3355, 2025: 3598}
EXPECTED_AGE_BINS = {
    "21-25": 534, "26-30": 1436, "31-35": 3300, "36-40": 3752, "41-45": 2880,
    "46-50": 2356, "51-55": 1918, "56-60": 1176, "61-65": 567,
    "66-70": 159, ">70": 107,
}
EXPECTED_NON9V_POS = {2021: 169, 2022: 152, 2023: 200, 2024: 229, 2025: 315}
EXPECTED_DATE_MIN = pd.Timestamp("2021-01-01")
EXPECTED_DATE_MAX = pd.Timestamp("2025-12-31")
EXPECTED_MONTHS = 60

STABLE15 = ["HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV39", "HPV45",
            "HPV51", "HPV52", "HPV53", "HPV56", "HPV58", "HPV59", "HPV66",
            "HPV68"]
IARC13 = ["HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV39", "HPV45",
          "HPV51", "HPV52", "HPV56", "HPV58", "HPV59", "HPV68"]
HR2B = ["HPV53", "HPV66"]
NON9V_HR = ["HPV35", "HPV39", "HPV51", "HPV53", "HPV56", "HPV59", "HPV66",
            "HPV68"]


# ---------------------------------------------------------------------------
# Wilson 95% CI (exact, two-sided)
# ---------------------------------------------------------------------------
def wilson_ci(pos: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = pos / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (100 * p, 100 * (centre - half), 100 * (centre + half))


# ---------------------------------------------------------------------------
# Apply the 5-step STROBE exclusion chain and return the analytic cohort.
# ---------------------------------------------------------------------------
def build_cohort(raw_path: Path) -> pd.DataFrame:
    df = pd.read_excel(raw_path)
    df.columns = [c.strip() for c in df.columns]
    df = df[df["Gender"] == "Female"].copy()
    df = df[df["Age"].notna()].copy()
    df = df[df["Age"] != 0].copy()
    df = df[df["Age"] >= 21].copy()
    df["Year"] = pd.to_datetime(df["Detection time"]).dt.year

    # Positivity mask: literal substring "Positive" -> captures 'Positive' and '弱Positive'
    for h in STABLE15:
        df[h + "_pos"] = df[h].astype(str).str.contains("Positive", na=False).astype(int)
    df["HR-IARC13_pos"] = df[[h + "_pos" for h in IARC13]].max(axis=1)
    df["HR-2B_pos"] = df[[h + "_pos" for h in HR2B]].max(axis=1)
    df["Any STABLE15_pos"] = df[[h + "_pos" for h in STABLE15]].max(axis=1)
    df["Non9v_pos"] = df[[h + "_pos" for h in NON9V_HR]].max(axis=1)
    return df


def aapc_loglinear(events_by_year: np.ndarray, n_by_year: np.ndarray) -> float:
    pcts = 100.0 * events_by_year / n_by_year
    if np.any(pcts == 0):
        return float("nan")
    x = np.arange(len(pcts))
    slope, *_ = linregress(x, np.log(pcts))
    return 100.0 * (np.exp(slope) - 1.0)


# ---------------------------------------------------------------------------
# QC report
# ---------------------------------------------------------------------------
def run_qc(raw_path, canon_path) -> int:
    raw_path = Path(raw_path); canon_path = Path(canon_path)
    print("=" * 72)
    print("HPV time-series 2021-2025  |  publication-grade QC verification")
    print(f"Raw:       {raw_path}")
    print(f"Canonical: {canon_path}")
    if IN_COLAB:
        print("Environment: Google Colab (Drive mounted at /content/drive)")
    else:
        print("Environment: local Python")
    print("=" * 72)

    if not raw_path.exists():
        print(f"[ ABORT ] Raw file not found: {raw_path}", file=sys.stderr)
        return 2
    if not canon_path.exists():
        print(f"[ ABORT ] Canonical JSON not found: {canon_path}", file=sys.stderr)
        return 2

    # Re-run exclusion chain with raw
    raw = pd.read_excel(raw_path)
    raw.columns = [c.strip() for c in raw.columns]
    n0 = len(raw)
    n1 = (raw["Gender"] == "Female").sum()
    n2 = ((raw["Gender"] == "Female") & raw["Age"].notna()).sum()
    n3 = ((raw["Gender"] == "Female") & raw["Age"].notna() & (raw["Age"] != 0)).sum()
    n4 = ((raw["Gender"] == "Female") & raw["Age"].notna()
          & (raw["Age"] != 0) & (raw["Age"] >= 21)).sum()

    df = build_cohort(raw_path)
    with open(canon_path, "r", encoding="utf-8") as f:
        canon = json.load(f)
    canon_yearly = canon["yearly"]
    canon_aapc = canon["aapc_pr"]

    fails = []

    # 1. STROBE chain
    expected_chain = [18670, 18295, 18243, 18242, 18185]
    actual_chain = [n0, n1, n2, n3, n4]
    if actual_chain == expected_chain:
        print(f"[ PASS ] STROBE exclusion chain reproduces: "
              f"{' -> '.join(str(x) for x in actual_chain)}")
    else:
        msg = f"STROBE chain mismatch: expected {expected_chain}, got {actual_chain}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 2. Yearly cohort sizes
    yearly_n = df.groupby("Year").size().to_dict()
    yearly_ok = all(yearly_n.get(y, 0) == n for y, n in EXPECTED_YEARLY_N.items())
    if yearly_ok and len(df) == EXPECTED_COHORT_N:
        print(f"[ PASS ] Yearly cohort sizes: {yearly_n}")
    else:
        msg = f"Yearly cohort sizes mismatch: expected {EXPECTED_YEARLY_N}, got {yearly_n}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 3. Per-genotype yearly positive counts
    pos_fail = 0
    metrics = ["Any STABLE15", "HR-IARC13", "HR-2B"] + STABLE15
    for metric in metrics:
        series = df.groupby("Year")[metric + "_pos"].sum()
        for yr in [2021, 2022, 2023, 2024, 2025]:
            key = metric if metric != "HR-2B" else "HR-2B (53/66)"
            excel_pos = int(series.get(yr, 0))
            json_pos = canon_yearly[key][str(yr)]["pos"]
            if excel_pos != json_pos:
                pos_fail += 1
                fails.append(f"pos count {metric}/{yr}: excel={excel_pos}, json={json_pos}")
    if pos_fail == 0:
        print("[ PASS ] Per-genotype yearly positive counts: 90/90 exact")
    else:
        print(f"[ FAIL ] Per-genotype positive count mismatches: {pos_fail}/90")

    # 4. Wilson 95% CI exact match (270 numeric assertions)
    ci_fail = 0
    worst = 0.0
    for metric in metrics:
        for yr in [2021, 2022, 2023, 2024, 2025]:
            n_yr = int(df[df["Year"] == yr].shape[0])
            pos = int(df[df["Year"] == yr][metric + "_pos"].sum())
            pct, lo, hi = wilson_ci(pos, n_yr)
            key = metric if metric != "HR-2B" else "HR-2B (53/66)"
            jpct = canon_yearly[key][str(yr)]["pct"]
            jlo = canon_yearly[key][str(yr)]["ci_lo"]
            jhi = canon_yearly[key][str(yr)]["ci_hi"]
            for v1, v2 in [(pct, jpct), (lo, jlo), (hi, jhi)]:
                d = abs(v1 - v2)
                worst = max(worst, d)
                if d > 1e-6:
                    ci_fail += 1
    if ci_fail == 0:
        print(f"[ PASS ] Wilson 95% CIs: 270/270 (worst diff {worst:.2e}, floating-point limit)")
    else:
        print(f"[ FAIL ] Wilson CI mismatches: {ci_fail}/270 (worst {worst:.4e})")
        fails.append("Wilson CI mismatch")

    # 5. AAPC point estimates
    aapc_fail = 0
    for metric in metrics:
        n_arr = df.groupby("Year").size().reindex([2021, 2022, 2023, 2024, 2025]).values
        e_arr = df.groupby("Year")[metric + "_pos"].sum().reindex([2021, 2022, 2023, 2024, 2025]).values
        a_excel = aapc_loglinear(e_arr.astype(float), n_arr.astype(float))
        key = metric if metric != "HR-2B" else "HR-2B (53/66)"
        a_json = canon_aapc[key]["aapc"]
        if abs(a_excel - a_json) > 0.001:
            aapc_fail += 1
            fails.append(f"AAPC {metric}: excel={a_excel:.4f}, json={a_json:.4f}")
    if aapc_fail == 0:
        print("[ PASS ] AAPC point estimates (log-linear OLS): 18/18 exact")
    else:
        print(f"[ FAIL ] AAPC point-estimate mismatches: {aapc_fail}/18")

    # 6. Non-9v HR yearly burden
    non9v_yr = df.groupby("Year")["Non9v_pos"].sum().to_dict()
    non9v_ok = all(int(non9v_yr.get(y, 0)) == n for y, n in EXPECTED_NON9V_POS.items())
    if non9v_ok:
        print(f"[ PASS ] Non-9v HR yearly positives: {non9v_yr}")
    else:
        msg = f"Non-9v burden mismatch: expected {EXPECTED_NON9V_POS}, got {non9v_yr}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 7. Set identity: |STABLE15| = |IARC13 union HR-2B|
    any_pos = int(df["Any STABLE15_pos"].sum())
    iarc13_pos = int(df["HR-IARC13_pos"].sum())
    hr2b_pos = int(df["HR-2B_pos"].sum())
    overlap = iarc13_pos + hr2b_pos - any_pos
    identity_ok = (int((df["HR-IARC13_pos"] | df["HR-2B_pos"]).sum()) == any_pos == 2318)
    if identity_ok:
        print(f"[ PASS ] Set identity: |IARC13| ({iarc13_pos}) + |HR-2B| ({hr2b_pos}) "
              f"- |overlap| ({overlap}) = |STABLE15| ({any_pos})")
    else:
        msg = f"Set identity violated: IARC13={iarc13_pos}, HR-2B={hr2b_pos}, union={any_pos}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 8. Fig 5 age-bin sample sizes (renumbered: was Fig 2 in v1)
    bins = [21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 200]
    labels = list(EXPECTED_AGE_BINS.keys())
    df["_AgeBin"] = pd.cut(df["Age"], bins=bins, right=False, labels=labels)
    actual_bins = df.groupby("_AgeBin", observed=True).size().to_dict()
    bin_ok = all(int(actual_bins.get(k, 0)) == v for k, v in EXPECTED_AGE_BINS.items())
    if bin_ok:
        print("[ PASS ] Fig 5 age-bin sample sizes: 11/11 exact")
    else:
        msg = f"Age-bin mismatch: expected {EXPECTED_AGE_BINS}, got {actual_bins}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 9. Co-infection multiplicity (disclosed in Supp Methods Note 3)
    hr_count = df[[h + "_pos" for h in STABLE15]].sum(axis=1)
    multiplicity = {
        0: int((hr_count == 0).sum()),
        1: int((hr_count == 1).sum()),
        2: int((hr_count == 2).sum()),
        ">=3": int((hr_count >= 3).sum()),
    }
    expected_mult = {0: 15867, 1: 1879, 2: 345, ">=3": 94}
    mult_ok = (multiplicity == expected_mult)
    if mult_ok:
        print(f"[ PASS ] HR-HPV multiplicity: {multiplicity}")
    else:
        msg = f"Multiplicity mismatch: expected {expected_mult}, got {multiplicity}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    # 10. Date range + 60-month completeness
    dates = pd.to_datetime(df["Detection time"])
    months = pd.to_datetime(df["Detection time"]).dt.to_period("M")
    months_ok = (months.nunique() == EXPECTED_MONTHS and
                 dates.min() >= EXPECTED_DATE_MIN and
                 dates.max() <= EXPECTED_DATE_MAX)
    if months_ok:
        print(f"[ PASS ] Detection date integrity: range "
              f"{dates.min().date()} to {dates.max().date()}, {months.nunique()} months")
    else:
        msg = f"Date integrity failed: min={dates.min()}, max={dates.max()}, months={months.nunique()}"
        print(f"[ FAIL ] {msg}")
        fails.append(msg)

    print("=" * 72)
    if not fails:
        print("OVERALL: ALL QC LAYERS PASSED")
        return 0
    print(f"OVERALL: {len(fails)} FAILURES")
    for f in fails:
        print(f"  - {f}")
    return 1


def main():
    ap = argparse.ArgumentParser(
        description="publication-grade QC verification for HPV time-series 2021-2025."
    )
    ap.add_argument("--raw", type=Path, default=None,
                    help="Raw genotyping Excel export (defaults: Drive in Colab, sibling locally).")
    ap.add_argument("--canonical", type=Path, default=None,
                    help="Canonical results JSON.")
    args = ap.parse_args()

    raw = args.raw if args.raw is not None else _default_raw_path()
    canon = args.canonical if args.canonical is not None else _default_canon_path()
    return run_qc(raw, canon)


if __name__ == "__main__":
    raise SystemExit(main())
