"""Validate P3 preregistration without reading residual/result data."""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("results/p3")
formation = pd.read_csv(ROOT / "P3_FORMATION.csv")
test = pd.read_csv(ROOT / "P3_TEST.csv")
prop = (ROOT / "P3_PROPOSITION.md").read_text(encoding="utf-8")
p2f = set(pd.read_csv("results/p2/P2_FORMATION.csv").DMS_id)
failures = []

if len(test) != 26 or len(formation) != 165:
    failures.append(f"wrong split sizes: {len(test)} + {len(formation)}")
if set(test.DMS_id) & set(formation.DMS_id):
    failures.append("DMS_id overlap")
if not set(test.DMS_id) <= p2f:
    failures.append("P3 test set not confined to the never-touched pool")
if (test["n_doubles_both_singles"] < 32).any():
    failures.append("test assay below the 32-double eligibility floor")

required = [
    "−0.1390",
    "[−0.1767, −0.0839]",
    "+0.0699",
    "[+0.0012, +0.1289]",
    "无法判断",
    "ε_m = DMS_double,m − DMS_s1,m − DMS_s2,m",
    "LLR_norm_dbl = LLR_{i1}/R_{i1} + LLR_{i2}/R_{i2}",
    "ddof=0",
    "n_dbl ≥ 32",
    "26 个",
    "130,340",
    "76%",
    "不放宽阈值",
    "P3_TEST.csv",
    "10,000 次",
    "CI 下限 ≤ 0.10",
]
for text in required:
    if text not in prop:
        failures.append(f"proposition missing: {text}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    "P3 freeze OK: 26 test / 165 formation assays; test confined to clean pool; "
    "cites P1/P2 numbers, states 76% concentration and power limits, threshold 0.10"
)
