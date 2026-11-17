"""Validate P2 preregistration without reading residual/result data."""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("results/p2")
formation = pd.read_csv(ROOT / "P2_FORMATION.csv")
test = pd.read_csv(ROOT / "P2_TEST.csv")
prop = (ROOT / "P2_PROPOSITION.md").read_text(encoding="utf-8")
p1_test = set(pd.read_csv("results/p1/P1_TEST.csv").DMS_id)
p1_formation = set(pd.read_csv("results/p1/P1_FORMATION.csv").DMS_id)
failures = []

if len(formation) != 78 or len(test) != 48:
    failures.append(f"wrong split sizes: {len(formation)} + {len(test)}")
if set(formation.DMS_id) & set(test.DMS_id):
    failures.append("DMS_id overlap")
if set(formation.UniProt_ID) & set(test.UniProt_ID):
    failures.append("UniProt_ID overlap")
if not set(test.DMS_id) <= p1_formation:
    failures.append("P2 test set not contained in P1 formation pool")
if not (set(formation.DMS_id) | set(test.DMS_id)) <= (p1_formation | p1_test):
    failures.append("P2 lists contain assays outside the frozen 191")

required = [
    "−0.1390",
    "[−0.1767, −0.0839]",
    "50/61",
    "R_i = max",
    "LLR_norm = LLR / R_i",
    "ddof=0",
    "|θ_norm| ≤ 0.05",
    "ci_high < −0.05 或 ci_low > 0.05",
    "P2_TEST.csv",
    "10,000 次",
]
for text in required:
    if text not in prop:
        failures.append(f"proposition missing: {text}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    "P2 freeze OK: 78 formation / 48 test assays; test pool confined to "
    "P1 formation; proposition cites P1 numbers and sets 0.05 null band"
)
