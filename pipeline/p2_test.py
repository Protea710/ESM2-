"""Execute frozen proposition P2 on the confirmed test set.

Order of operations, exactly as frozen in results/p2/P2_PROPOSITION.md:
1. verify frozen artifact hashes against results/p2/P2_FREEZE.json (abort on mismatch)
2. recompute WT-context masked logits for the test assays' mutated positions
   (model-internal only; produces no experimental residual), derive
   R_i = max over the 19 non-WT aa of |logp(aa) - logp(wt)|, and cross-check
   the recomputed logp(wt) against the stored positions table (determinism guard)
3. per assay: normalized model score = LLR / R_i; residuals r = z_exp - z_model
   over valid singles (ddof=0); eligibility >=8 singles/position, >=8 positions/assay
4. margin quartiles (same grouping as P1); Delta_a = mean(e_low) - mean(e_high);
   theta_norm = median; 10,000-assay bootstrap (seed 20261117)
5. verdict: SUPPORTED iff ci_low >= -0.05 and ci_high <= 0.05;
   REFUTED iff ci_high < -0.05 or ci_low > 0.05; else INDETERMINATE
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline.p1_test import verify_freeze  # noqa: E402
from pipeline.scoring_math import zscore  # noqa: E402


def assay_delta_norm(dms_id, mut_csv, pos_csv, ri_csv):
    mut = pd.read_csv(mut_csv)
    singles = mut[
        ~mut["mutant"].astype(str).str.contains(":")
        & mut["DMS_score"].notna()
        & mut["model_score"].notna()
    ].copy()
    if len(singles) < 8:
        return None, f"too_few_singles({len(singles)})"
    singles["pos"] = singles["mutant"].map(lambda m: int(str(m)[1:-1]))
    ri = pd.read_csv(ri_csv)  # position, R_i
    bad = ri[(ri["R_i"] <= 0) | ~np.isfinite(ri["R_i"])]
    if len(bad):
        singles = singles[~singles["pos"].isin(set(bad["position"]))]
    ri = ri[(ri["R_i"] > 0) & np.isfinite(ri["R_i"])].rename(
        columns={"position": "pos"}
    )
    singles = singles.merge(ri, on="pos", how="inner")
    if len(singles) < 8:
        return None, "too_few_singles_after_ri"
    singles["model_norm"] = singles["model_score"] / singles["R_i"]
    singles["r"] = zscore(singles["DMS_score"].to_numpy(float)) - zscore(
        singles["model_norm"].to_numpy(float)
    )
    singles["absr"] = singles["r"].abs()
    agg = singles.groupby("pos")["absr"].agg(["mean", "count"])
    agg = agg[agg["count"] >= 8].rename(columns={"mean": "e_i", "count": "n_single"})
    if len(agg) < 8:
        return None, f"too_few_eligible_positions({len(agg)})"

    posf = pd.read_csv(pos_csv)[["position", "margin"]]
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


def compute_ri(engine, dms_id, ref, data_dir, out_csv):
    """Recompute WT-context logits for mutated positions; return (R_i table, max
    |recomputed logp(wt) - stored logp_wt|) as a determinism guard."""
    row = ref.loc[ref.DMS_id == dms_id].iloc[0]
    seq = str(row.target_seq).upper()
    df = pd.read_csv(data_dir / f"{dms_id}.csv", usecols=["mutant"])
    positions = sorted(
        {int(p[1:-1]) - 1 for m in df["mutant"] for p in str(m).split(":")}
    )
    logits = engine.position_logits(seq, positions)
    posf = pd.read_csv(
        pathlib.Path("/work/results_a/positions") / f"{dms_id}.csv",
        usecols=["position", "logp_wt"],
    )
    stored = posf.set_index("position")["logp_wt"]
    rows, max_dev = [], 0.0
    for r, p in enumerate(positions):
        lp = logits[r]
        wt_lp = lp[engine.aa_to_col[seq[p]]]
        others = np.delete(lp, engine.aa_to_col[seq[p]])
        R = float(np.max(np.abs(others - wt_lp)))
        stored_lp = stored.get(p + 1, np.nan)
        if np.isfinite(stored_lp):
            max_dev = max(max_dev, abs(wt_lp - stored_lp))
        rows.append({"position": p + 1, "R_i": R})
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    return out, max_dev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="/work/esm2_model")
    ap.add_argument(
        "--data-dir",
        default="/work/DMS_ProteinGym_substitutions/DMS_ProteinGym_substitutions",
    )
    ap.add_argument("--ref-csv", default="/work/reference/DMS_substitutions.csv")
    ap.add_argument("--positions-dir", default="/work/results_a/positions")
    ap.add_argument("--p2-dir", default="/work/p2")
    ap.add_argument("--out", default="/work/results_p2")
    args = ap.parse_args()
    p2 = pathlib.Path(args.p2_dir)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ri").mkdir(exist_ok=True)

    verify_freeze(p2 / "P2_FREEZE.json", p2)
    from pipeline.esm_masked import MaskedMarginalEngine

    test = pd.read_csv(p2 / "P2_TEST.csv")
    ref = pd.read_csv(args.ref_csv)
    engine = MaskedMarginalEngine(args.model_dir)

    rows, excluded = [], []
    max_dev_overall = 0.0
    for did in test["DMS_id"]:
        ri_csv = out / "ri" / f"{did}.csv"
        if not ri_csv.exists():
            _, dev = compute_ri(engine, did, ref, pathlib.Path(args.data_dir), ri_csv)
            max_dev_overall = max(max_dev_overall, dev)
        d, reason = assay_delta_norm(
            did,
            pathlib.Path("/work/results_a/mutants") / f"{did}.csv",
            pathlib.Path(args.positions_dir) / f"{did}.csv",
            ri_csv,
        )
        if d is None:
            excluded.append({"DMS_id": did, "reason": reason})
        else:
            rows.append(d)
    per = pd.DataFrame(rows).sort_values("DMS_id").reset_index(drop=True)
    per.to_csv(out / "P2_per_assay.csv", index=False)
    pd.DataFrame(excluded).to_csv(out / "P2_excluded.csv", index=False)

    deltas = per["delta_a"].to_numpy()
    theta = float(np.median(deltas))
    rng = np.random.default_rng(20261117)
    n = len(deltas)
    boots = np.array([np.median(deltas[rng.integers(0, n, n)]) for _ in range(10000)])
    ci_low, ci_high = (
        float(np.percentile(boots, 2.5)),
        float(np.percentile(boots, 97.5)),
    )
    if ci_low >= -0.05 and ci_high <= 0.05:
        verdict = "SUPPORTED"
    elif ci_high < -0.05 or ci_low > 0.05:
        verdict = "REFUTED"
    else:
        verdict = "INDETERMINATE"

    summary = {
        "n_test_assays": int(len(test)),
        "n_eligible": int(n),
        "n_excluded": len(excluded),
        "ri_logpwt_max_dev": max_dev_overall,
        "theta_norm": theta,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "null_band": [-0.05, 0.05],
        "rule": "SUPPORTED iff ci_low>=-0.05 and ci_high<=0.05; REFUTED iff ci_high<-0.05 or ci_low>0.05; else INDETERMINATE",
        "verdict": verdict,
        "bootstrap_seed": 20261117,
        "bootstrap_reps": 10000,
    }
    (out / "P2_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(per.sort_values("delta_a").head(5).to_string(index=False))
    print(per.sort_values("delta_a").tail(5).to_string(index=False))


if __name__ == "__main__":
    main()
