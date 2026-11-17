"""Validate P2 results against the frozen proposition (no data recomputation)."""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

failures = []

for manifest, artifacts in [
    ("results/p2/P2_FREEZE.json", None),
]:
    freeze = json.loads(Path(manifest).read_text(encoding="utf-8"))
    for art in freeze["artifacts"]:
        got = hashlib.sha256(Path(art["path"]).read_bytes()).hexdigest()
        if got != art["sha256"]:
            failures.append(f"frozen file modified: {art['path']}")

s = json.loads(Path("results/p2/P2_summary.json").read_text(encoding="utf-8"))
per = pd.read_csv("results/p2/P2_per_assay.csv")
excluded = pd.read_csv("results/p2/P2_excluded.csv")

if s["n_eligible"] != len(per):
    failures.append(f"n_eligible {s['n_eligible']} != rows {len(per)}")
if s["n_eligible"] + s["n_excluded"] != 48:
    failures.append("eligible + excluded != 48")
if abs(float(per["delta_a"].median()) - s["theta_norm"]) > 1e-9:
    failures.append("theta_norm drift vs per-assay median")
rule_verdict = (
    "SUPPORTED"
    if s["ci_low"] >= -0.05 and s["ci_high"] <= 0.05
    else "REFUTED"
    if (s["ci_high"] < -0.05 or s["ci_low"] > 0.05)
    else "INDETERMINATE"
)
if s["verdict"] != rule_verdict:
    failures.append(f"verdict {s['verdict']} != rule outcome {rule_verdict}")

res_mtime = Path("results/p2/P2_RESULTS.md").stat().st_mtime
freeze_mtime = Path("results/p2/P2_FREEZE.json").stat().st_mtime
if res_mtime <= freeze_mtime:
    failures.append("P2_RESULTS.md not written after the freeze manifest")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    f"P2 results OK: frozen artifacts intact; theta_norm={s['theta_norm']:.4f} "
    f"CI=[{s['ci_low']:.4f},{s['ci_high']:.4f}]; verdict {s['verdict']}; excluded={len(excluded)}"
)
