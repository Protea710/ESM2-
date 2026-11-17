"""Execute frozen proposition P3 (multi-site non-additivity) on the confirmed test set.

Frozen protocol (results/p3/P3_PROPOSITION.md), executed in order:
1. verify frozen artifact hashes against results/p3/P3_FREEZE.json
2. recompute R_i for the test assays' mutated positions (model-internal only),
   with logp(wt) determinism cross-check against the stored positions table
3. per eligible assay (n_dbl >= 32 valid doubles whose both component singles
   exist and all three mutations have valid DMS_score and model_score):
     epsilon_m = DMS_double - DMS_s1 - DMS_s2 (raw scale)
     eps_z = (eps - mean)/sd over the assay's double set (ddof=0; sd=0 -> exclude)
     LLR_norm_dbl = LLR_1/R_i1 + LLR_2/R_i2
     r_m = z(DMS_double) - z(LLR_norm_dbl) over the double set (ddof=0)
4. sort doubles by (|eps_z|, mutant) ascending; low/high quartiles = head/tail
   ceil(n/4); Delta_a = mean|r|_high - mean|r|_low; theta = median over assays;
   10,000-assay bootstrap (seed 20261117)
5. verdict: REFUTED iff ci_low <= 0.10, else SUPPORTED
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline.p1_test import verify_freeze  # noqa: E402
from pipeline.p2_test import compute_ri  # noqa: E402
from pipeline.scoring_math import zscore  # noqa: E402

MIN_DOUBLES = 32


def assay_delta_nonadd(dms_id, mut_csv, ri_csv):
    mut = pd.read_csv(mut_csv)
    singles = set(
        mut.loc[~mut["mutant"].astype(str).str.contains(":"), "mutant"].astype(str)
    )
    dbl = mut[mut["mutant"].astype(str).str.count(":") == 1].copy()
    ok = (
        dbl["DMS_score"].notna()
        & dbl["model_score"].notna()
        & dbl["mutant"].astype(str).map(lambda m: m.split(":")[0] in singles)
        & dbl["mutant"].astype(str).map(lambda m: m.split(":")[1] in singles)
    )
    dbl = dbl[ok].copy()
    if len(dbl) < MIN_DOUBLES:
        return None, f"too_few_valid_doubles({len(dbl)})"

    dms_of = mut.dropna(subset=["DMS_score"]).set_index("mutant")["DMS_score"]
    score_of = mut.dropna(subset=["model_score"]).set_index("mutant")["model_score"]
    dbl["s1"] = dbl["mutant"].map(lambda m: m.split(":")[0])
    dbl["s2"] = dbl["mutant"].map(lambda m: m.split(":")[1])
    dbl["eps"] = dbl["DMS_score"] - dbl["s1"].map(dms_of) - dbl["s2"].map(dms_of)
    sd = dbl["eps"].std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return None, f"degenerate_eps_sd({sd})"
    dbl["eps_z"] = (dbl["eps"] - dbl["eps"].mean()) / sd

    ri = pd.read_csv(ri_csv).rename(columns={"position": "pos"})
    ri = ri[(ri["R_i"] > 0) & np.isfinite(ri["R_i"])]
    ri_of = ri.set_index("pos")["R_i"]

    def norm_score(m: str) -> float:
        total = 0.0
        for part in m.split(":"):
            pos = int(part[1:-1])
            if pos not in ri_of.index or part not in score_of.index:
                return np.nan
            total += float(score_of[part]) / float(ri_of[pos])
        return total

    dbl["model_norm"] = dbl["mutant"].map(norm_score)
    dbl = dbl.dropna(subset=["model_norm"])
    if len(dbl) < MIN_DOUBLES:
        return None, f"too_few_doubles_after_ri({len(dbl)})"

    dbl["r"] = zscore(dbl["DMS_score"].to_numpy(float)) - zscore(
        dbl["model_norm"].to_numpy(float)
    )
    dbl["absr"] = dbl["r"].abs()
    dbl = dbl.assign(eps_z_abs=dbl["eps_z"].abs()).sort_values(
        ["eps_z_abs", "mutant"], kind="stable"
    )
    k = int(np.ceil(len(dbl) / 4))
    low, high = dbl.head(k), dbl.tail(k)
    return {
        "DMS_id": dms_id,
        "n_dbl": len(dbl),
        "k": k,
        "mean_absr_low": float(low["absr"].mean()),
        "mean_absr_high": float(high["absr"].mean()),
        "delta_a": float(high["absr"].mean() - low["absr"].mean()),
    }, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="/work/esm2_model")
    ap.add_argument("--mutants-dir", default="/work/results_a/mutants")
    ap.add_argument(
        "--data-dir",
        default="/work/DMS_ProteinGym_substitutions/DMS_ProteinGym_substitutions",
    )
    ap.add_argument("--ref-csv", default="/work/reference/DMS_substitutions.csv")
    ap.add_argument("--p3-dir", default="/work/p3")
    ap.add_argument("--out", default="/work/results_p3")
    args = ap.parse_args()
    p3 = pathlib.Path(args.p3_dir)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ri").mkdir(exist_ok=True)

    verify_freeze(p3 / "P3_FREEZE.json", p3)
    from pipeline.esm_masked import MaskedMarginalEngine

    test = pd.read_csv(p3 / "P3_TEST.csv")
    ref = pd.read_csv(args.ref_csv)
    engine = MaskedMarginalEngine(args.model_dir)

    rows, excluded = [], []
    max_dev_overall = 0.0
    for did in test["DMS_id"]:
        ri_csv = out / "ri" / f"{did}.csv"
        if not ri_csv.exists():
            _, dev = compute_ri(engine, did, ref, pathlib.Path(args.data_dir), ri_csv)
            max_dev_overall = max(max_dev_overall, dev)
        d, reason = assay_delta_nonadd(
            did, pathlib.Path(args.mutants_dir) / f"{did}.csv", ri_csv
        )
        if d is None:
            excluded.append({"DMS_id": did, "reason": reason})
        else:
            rows.append(d)
    per = pd.DataFrame(rows).sort_values("DMS_id").reset_index(drop=True)
    per.to_csv(out / "P3_per_assay.csv", index=False)
    pd.DataFrame(excluded).to_csv(out / "P3_excluded.csv", index=False)

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
        "ri_logpwt_max_dev": max_dev_overall,
        "theta_3": theta,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "rule": "REFUTED iff ci_low <= 0.10, else SUPPORTED",
        "verdict": verdict,
        "bootstrap_seed": 20261117,
        "bootstrap_reps": 10000,
    }
    (out / "P3_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(per.sort_values("delta_a").head(5).to_string(index=False))
    print(per.sort_values("delta_a").tail(5).to_string(index=False))


if __name__ == "__main__":
    main()
