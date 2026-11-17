"""Validate the synced A-segment results against the official benchmark.
Runs locally on results/ copies synced from the lab machine."""

import sys

import pandas as pd

failures = []

summary = pd.read_csv("results/assay_summary.csv")
compare = pd.read_csv("results/compare.csv")
official = pd.read_csv(
    "data/reference/DMS_substitutions_Spearman_DMS_level_esm2_650M.csv"
)
col = [c for c in official.columns if "650M" in c][0]

if len(summary) != 191:
    failures.append(f"summary rows {len(summary)} != 191")
if summary["DMS_id"].duplicated().any():
    failures.append("duplicated DMS_id in summary")
if compare["rho_official"].isna().any():
    failures.append("unmatched rows in compare")

mean_ours = compare["rho_all"].mean()
mean_off = compare["rho_official"].mean()
if abs(mean_ours - mean_off) > 0.01:
    failures.append(
        f"mean rho drift {mean_ours:.4f} vs official same-assay mean {mean_off:.4f}"
    )

n_bad = int((compare["delta"].abs() > 0.05).sum())
if n_bad > 0:
    failures.append(f"{n_bad} assays with |delta| > 0.05 (unexplained)")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    f"A-segment OK: 191 assays, mean rho {mean_ours:.4f} vs official {mean_off:.4f}, "
    f"max |delta| {compare['delta'].abs().max():.4f}"
)
