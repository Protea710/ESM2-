"""Validate the final report against the archived result files (no recomputation)."""

import json
import sys
from pathlib import Path

import pandas as pd

failures = []
report = Path("docs/FINAL_REPORT.md").read_text(encoding="utf-8")

p1 = json.loads(Path("results/p1/P1_summary.json").read_text())
p2 = json.loads(Path("results/p2/P2_summary.json").read_text())
p3 = json.loads(Path("results/p3/P3_summary.json").read_text())
compare = pd.read_csv("results/compare.csv")

required = [
    "−0.1390",
    "[−0.1767, −0.0839]",
    "REFUTED",
    "+0.0699",
    "[+0.0012, +0.1289]",
    "INDETERMINATE",
    "+0.1366",
    "[+0.0022, +0.2904]",
    "3b8ae8f3536c39a0ddbc87a6132e6444113b9bc16a6399a1cf90c977b17bbf89",
    "ef0e814e43008da4befddd02b7b88f628042476f0842ff9c585b94feee3a3865",
    "618840c25f90209a8e275f05c75391d31eb8d18b559727f5bdf1bab32c46b607",
    "CSV 误读",
    "门禁误删",
    "run_all.sh",
    "requirements-lock.txt",
    "20261117",
    "P4",
    "0.45042",
    "0.000506",
]
for text in required:
    if text not in report:
        failures.append(f"report missing: {text}")

if abs(compare["rho_all"].mean() - 0.45042) > 5e-5:
    failures.append("report mean rho drift vs compare.csv")

verdicts = {"P1": p1["verdict"], "P2": p2["verdict"], "P3": p3["verdict"]}
expected = {"P1": "REFUTED", "P2": "INDETERMINATE", "P3": "REFUTED"}
if verdicts != expected:
    failures.append(f"verdict drift: {verdicts}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    "FINAL_REPORT OK: verdicts, hashes, corrections, rerun command, seeds, P4 present"
)
