"""Execute frozen proposition P1 on the confirmed test set.

Steps, in order:
1. verify frozen artifact hashes against results/p1/P1_FREEZE.json (abort on mismatch)
2. per test assay: recompute singles residuals per the frozen definition
   (z over that assay's valid singles, ddof=0), eligibility >=8 singles per
   position and >=8 eligible positions per assay
3. quartile split by (margin, position) ascending; Delta_a = mean(e_low)-mean(e_high)
4. theta = median Delta_a; 10,000 assay-level bootstrap (seed 20261117); verdict:
   REFUTED iff CI lower bound <= 0.10
"""

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline.scoring_math import zscore  # noqa: E402


def verify_freeze(freeze_path: pathlib.Path, artifact_dir: pathlib.Path):
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    for art in freeze["artifacts"]:
        name = pathlib.Path(art["path"]).name
        p = artifact_dir / name
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != art["sha256"]:
            raise SystemExit(f"FREEZE MISMATCH {name}: {got}")
    print(f"freeze verified: {len(freeze['artifacts'])} artifacts intact")
    return freeze


def assay_delta(dms_id: str, mutants_csv: pathlib.Path, positions_csv: pathlib.Path):
    mut = pd.read_csv(mutants_csv)
    singles = mut[
        ~mut["mutant"].astype(str).str.contains(":")
        & mut["DMS_score"].notna()
        & mut["model_score"].notna()
    ].copy()
    if len(singles) < 8:
        return None, f"too_few_singles({len(singles)})"
    singles["r"] = zscore(singles["DMS_score"].to_numpy(float)) - zscore(
        singles["model_score"].to_numpy(float)
    )
    singles["absr"] = singles["r"].abs()
    singles["pos"] = singles["mutant"].map(lambda m: int(str(m)[1:-1]))
    agg = singles.groupby("pos")["absr"].agg(["mean", "count"])
    agg = agg[agg["count"] >= 8].rename(columns={"mean": "e_i", "count": "n_single"})
    if len(agg) < 8:
        return None, f"too_few_eligible_positions({len(agg)})"

    posf = pd.read_csv(positions_csv)[["position", "margin"]]
    mg = agg.join(posf.set_index("position"), how="left").reset_index()
    mg = mg.rename(columns={mg.columns[0]: "position"})
    if mg["margin"].isna().any():
        return None, "missing_margin"
    mg = mg.sort_values(["margin", "position"], kind="stable")
    k = int(np.ceil(len(mg) / 4))
    low, high = mg.head(k), mg.tail(k)
    return {
        "DMS_id": dms_id,
        "n_eligible": len(mg),
        "k": k,
        "mean_e_low": float(low["e_i"].mean()),
        "mean_e_high": float(high["e_i"].mean()),
        "delta_a": float(low["e_i"].mean() - high["e_i"].mean()),
    }, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutants-dir", default="/work/results_a/mutants")
    ap.add_argument("--positions-dir", default="/work/results_a/positions")
    ap.add_argument("--p1-dir", default="/work/p1")
    ap.add_argument("--out", default="/work/results_p1")
    args = ap.parse_args()
    p1 = pathlib.Path(args.p1_dir)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    verify_freeze(p1 / "P1_FREEZE.json", p1)
    test = pd.read_csv(p1 / "P1_TEST.csv")

    rows, excluded = [], []
    for did in test["DMS_id"]:
        d, reason = assay_delta(
            did,
            pathlib.Path(args.mutants_dir) / f"{did}.csv",
            pathlib.Path(args.positions_dir) / f"{did}.csv",
        )
        if d is None:
            excluded.append({"DMS_id": did, "reason": reason})
        else:
            rows.append(d)
    per = pd.DataFrame(rows).sort_values("DMS_id").reset_index(drop=True)
    per.to_csv(out / "P1_per_assay.csv", index=False)
    pd.DataFrame(excluded).to_csv(out / "P1_excluded.csv", index=False)

    deltas = per["delta_a"].to_numpy()
    theta = float(np.median(deltas))
    rng = np.random.default_rng(20261117)
    n = len(deltas)
    boots = np.array([np.median(deltas[rng.integers(0, n, n)]) for _ in range(10000)])
    ci_low, ci_high = (
        float(np.percentile(boots, 2.5)),
        float(np.percentile(boots, 97.5)),
    )
    verdict = "REFUTED" if ci_low <= 0.10 else "SUPPORTED"

    summary = {
        "n_test_assays": int(len(test)),
        "n_eligible": int(n),
        "n_excluded": len(excluded),
        "theta": theta,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "rule": "REFUTED iff ci_low <= 0.10",
        "verdict": verdict,
        "bootstrap_seed": 20261117,
        "bootstrap_reps": 10000,
    }
    (out / "P1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("top deltas:")
    print(per.sort_values("delta_a").head(5).to_string(index=False))
    print(per.sort_values("delta_a").tail(5).to_string(index=False))


if __name__ == "__main__":
    main()
