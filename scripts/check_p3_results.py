"""Validate P3 results against the frozen proposition (no data recomputation)."""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

failures = []

freeze = json.loads(Path("results/p3/P3_FREEZE.json").read_text(encoding="utf-8"))
for art in freeze["artifacts"]:
    got = hashlib.sha256(Path(art["path"]).read_bytes()).hexdigest()
    if got != art["sha256"]:
        failures.append(f"frozen file modified: {art['path']}")

s = json.loads(Path("results/p3/P3_summary.json").read_text(encoding="utf-8"))
per = pd.read_csv("results/p3/P3_per_assay.csv")
excluded = pd.read_csv("results/p3/P3_excluded.csv")

if s["n_eligible"] != len(per):
    failures.append(f"n_eligible {s['n_eligible']} != rows {len(per)}")
if s["n_eligible"] + s["n_excluded"] != 26:
    failures.append("eligible + excluded != 26 test assays")
if abs(float(per["delta_a"].median()) - s["theta_3"]) > 1e-9:
    failures.append("theta_3 drift vs per-assay median")
rule_verdict = "REFUTED" if s["ci_low"] <= 0.10 else "SUPPORTED"
if s["verdict"] != rule_verdict:
    failures.append(f"verdict {s['verdict']} != rule outcome {rule_verdict}")
total_dbl = int(per["n_dbl"].sum())
if abs(total_dbl - 130340) > 0:
    failures.append(f"total valid doubles {total_dbl} != frozen census 130340")

res_mtime = Path("results/p3/P3_RESULTS.md").stat().st_mtime
freeze_mtime = Path("results/p3/P3_FREEZE.json").stat().st_mtime
if res_mtime <= freeze_mtime:
    failures.append("P3_RESULTS.md not written after the freeze manifest")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    f"P3 results OK: frozen artifacts intact; theta_3={s['theta_3']:.4f} "
    f"CI=[{s['ci_low']:.4f},{s['ci_high']:.4f}]; verdict {s['verdict']}; "
    f"excluded={len(excluded)}; doubles={total_dbl}"
)
