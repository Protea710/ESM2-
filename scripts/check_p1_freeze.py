"""Validate P1 preregistration without reading residual/result data."""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("results/p1")
formation = pd.read_csv(ROOT / "P1_FORMATION.csv")
test = pd.read_csv(ROOT / "P1_TEST.csv")
prop = (ROOT / "P1_PROPOSITION.md").read_text(encoding="utf-8")
failures = []

if len(formation) != 126 or len(test) != 65:
    failures.append(f"wrong split sizes: {len(formation)} + {len(test)}")
if set(formation.DMS_id) & set(test.DMS_id):
    failures.append("DMS_id overlap")
if set(formation.UniProt_ID) & set(test.UniProt_ID):
    failures.append("UniProt_ID overlap")
if set(formation.DMS_id) | set(test.DMS_id) != set(
    pd.read_csv("data/assay_list_proposed.csv").DMS_id
):
    failures.append("split union differs from frozen 191-assay list")

required = [
    "至少高 0.10",
    "ddof=0",
    "至少有 8 个有效单点突变",
    "θ = median_a(Δ_a)",
    "10,000 次",
    "区间下限 ≤ 0.10",
    "P1_TEST.csv",
    "检验集确认",
]
for text in required:
    if text not in prop:
        failures.append(f"proposition missing: {text}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    "P1 freeze OK: 126 formation assays / 65 test assays; "
    "no DMS or UniProt overlap; proposition threshold 0.10 present"
)
