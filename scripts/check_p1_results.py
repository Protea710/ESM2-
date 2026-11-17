"""Validate P1 results against the frozen proposition (no recomputation of data)."""

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

failures = []

freeze = json.loads(Path("results/p1/P1_FREEZE.json").read_text(encoding="utf-8"))
for art in freeze["artifacts"]:
    got = hashlib.sha256(Path(art["path"]).read_bytes()).hexdigest()
    if got != art["sha256"]:
        failures.append(f"frozen file modified: {art['path']}")

summary = json.loads(Path("results/p1/P1_summary.json").read_text(encoding="utf-8"))
per = pd.read_csv("results/p1/P1_per_assay.csv")
excluded = pd.read_csv("results/p1/P1_excluded.csv")

if summary["n_eligible"] != len(per):
    failures.append(f"summary n_eligible {summary['n_eligible']} != rows {len(per)}")
if summary["n_eligible"] + summary["n_excluded"] != 65:
    failures.append("eligible + excluded != 65 test assays")
theta = float(per["delta_a"].median())
if abs(theta - summary["theta"]) > 1e-9:
    failures.append(
        f"theta drift: per-assay median {theta} != summary {summary['theta']}"
    )
expected = "REFUTED" if summary["ci_low"] <= 0.10 else "SUPPORTED"
if summary["verdict"] != expected:
    failures.append(
        f"verdict {summary['verdict']} inconsistent with rule at ci_low {summary['ci_low']}"
    )

res_mtime = Path("results/p1/P1_RESULTS.md").stat().st_mtime
freeze_mtime = Path("results/p1/P1_FREEZE.json").stat().st_mtime
if res_mtime <= freeze_mtime:
    failures.append("P1_RESULTS.md not written after the freeze manifest")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    f"P1 results OK: frozen artifacts intact; theta={summary['theta']:.4f} "
    f"CI=[{summary['ci_low']:.4f},{summary['ci_high']:.4f}]; verdict {summary['verdict']}; "
    f"excluded={len(excluded)}"
)
