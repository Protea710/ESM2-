"""Score one assay under the official ProteinGym ESM masked-marginals protocol,
emitting per-mutant scores+residuals and per-position model-internal features."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .scoring_math import llr_at, parse_mutant, features_from_logprobs, zscore


def score_assay(
    engine, dms_id: str, seq: str, df: pd.DataFrame, mutant_col: str = "mutant"
):
    """df must contain mutant + DMS_score columns. Returns (mutant_df, position_df, info)."""
    L = len(seq)
    positions = sorted(
        {int(part[1:-1]) - 1 for m in df[mutant_col] for part in str(m).split(":")}
    )
    pos_arr = np.array(positions, dtype=int)
    logits = engine.position_logits(seq, pos_arr)  # (P, 20)

    pos_feature_rows = []
    for r, p in enumerate(pos_arr):
        f = features_from_logprobs(logits[r], seq[p], engine.aa_to_col)
        f["position"] = p + 1  # 1-based, ProteinGym convention
        pos_feature_rows.append(f)
    pos_feat = pd.DataFrame(pos_feature_rows)

    rows = []
    for mutant, dms_score in zip(df[mutant_col], df["DMS_score"]):
        try:
            subs = parse_mutant(str(mutant), seq)
        except (ValueError, IndexError):
            rows.append(
                {
                    "mutant": mutant,
                    "DMS_score": dms_score,
                    "model_score": np.nan,
                    "error": "parse",
                }
            )
            continue
        score = 0.0
        for wt, idx0, mt in subs:
            r = np.searchsorted(pos_arr, idx0)
            score += llr_at(logits[r], engine.aa_to_col, wt, mt)
        rows.append(
            {
                "mutant": mutant,
                "DMS_score": dms_score,
                "model_score": score,
                "error": "",
            }
        )
    mut = pd.DataFrame(rows)

    valid = mut["DMS_score"].notna() & mut["model_score"].notna()
    v = mut[valid]
    rho_all = spearmanr(v["model_score"], v["DMS_score"]).statistic
    singles = v[~v["mutant"].str.contains(":")]
    rho_single = spearmanr(singles["model_score"], singles["DMS_score"]).statistic

    mut["r"] = np.nan
    if len(v):
        mut.loc[valid, "r"] = zscore(v["DMS_score"].to_numpy()) - zscore(
            v["model_score"].to_numpy()
        )
    mut["n_pos"] = v["mutant"].map(
        lambda m: len(str(m).split(":")) if isinstance(m, str) else np.nan
    )

    # residue-level error over singles (recompute from mut after r exists);
    # pos is 1-based to align with pos_feat["position"]
    sv = mut[~mut["mutant"].astype(str).str.contains(":")].dropna(subset=["r"]).copy()
    sv["pos"] = sv["mutant"].map(lambda m: int(str(m)[1:-1]))
    sv["absr"] = sv["r"].abs()
    agg = sv.groupby("pos")["absr"].agg(["mean", "count"])
    pos_feat = pos_feat.merge(
        agg.rename(columns={"mean": "e_i", "count": "n_single"}),
        left_on="position",
        right_index=True,
        how="left",
    )
    pos_feat["pos_norm"] = (pos_feat["position"] - 1) / max(L - 1, 1)
    multi_touch = {}
    for m in v["mutant"]:
        if ":" in str(m):
            for part in str(m).split(":"):
                multi_touch[int(part[1:-1])] = multi_touch.get(int(part[1:-1]), 0) + 1
    pos_feat["n_multi"] = pos_feat["position"].map(lambda p: multi_touch.get(p, 0))

    info = {
        "DMS_id": dms_id,
        "seq_len": L,
        "n_rows": len(df),
        "n_valid": int(valid.sum()),
        "rho_all": float(rho_all),
        "rho_single": float(rho_single),
    }
    return mut, pos_feat, info
