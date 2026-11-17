"""Local unit tests for the pure scoring layer (no torch needed). Run: python3 scripts/test_scoring_math.py"""

import sys
import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from pipeline.scoring_math import (
    parse_mutant,
    features_from_logprobs,
    zscore,
    STANDARD_AA,
)
from pipeline.assay_scoring import score_assay


def test_parse():
    seq = "MAGKV"
    assert parse_mutant("A2G", seq) == [("A", 1, "G")]
    assert parse_mutant("M1A:V5G", seq) == [("M", 0, "A"), ("V", 4, "G")]
    try:
        parse_mutant("A2G", seq), parse_mutant("A1G", seq)
        raise SystemExit("expected wt mismatch error")
    except ValueError:
        pass
    print("parse OK")


def test_features_and_zscore():
    row = np.log(np.array([0.5, 0.2, 0.2] + [0.6 / 17] * 17))
    aa_col = {aa: j for j, aa in enumerate(STANDARD_AA)}
    f = features_from_logprobs(row, "A", aa_col)
    assert abs(f["logp_wt"] - np.log(0.5)) < 1e-9
    assert abs(f["top1_prob"] - 0.5) < 1e-9
    z = zscore(np.array([1.0, 2.0, 3.0]))
    assert abs(z.mean()) < 1e-12 and abs(z.std() - 1.0) < 1e-12
    print("features/zscore OK")


def test_score_assay_synthetic():
    """Engine stub: identity mapping of masked position -> hand-set distribution."""

    class StubEngine:
        aa_to_col = {aa: j for j, aa in enumerate(STANDARD_AA)}

        def position_logits(self, seq, positions):
            out = np.full((len(positions), 20), np.log(0.8 / 19))
            for r, p in enumerate(positions):
                out[r, self.aa_to_col[seq[p]]] = np.log(
                    0.2 + 0.01 * p
                )  # wt prob grows with position
            return out

    seq = "MAGKV"
    df = pd.DataFrame(
        {
            "mutant": ["M1A", "A2G", "K4R:V5L"],
            "DMS_score": [0.1, 0.9, 0.5],
        }
    )
    mut, posf, info = score_assay(StubEngine(), "TEST", seq, df)
    assert info["n_rows"] == 3 and info["n_valid"] == 3
    # M1A: wt prob log(0.2), mut log(0.8/19); LLR negative -> ranked against A2G (positive)
    m1 = mut.loc[mut.mutant == "M1A", "model_score"].iloc[0]
    assert m1 < 0
    assert (
        abs(posf.loc[posf.position == 2, "logp_wt"].iloc[0] - np.log(0.2 + 0.01)) < 1e-9
    )
    assert posf["pos_norm"].iloc[-1] == 1.0
    # alignment: e_i must equal the mean |r| of singles at that same 1-based position
    singles = (
        mut[~mut["mutant"].astype(str).str.contains(":")].dropna(subset=["r"]).copy()
    )
    singles["pos"] = singles["mutant"].map(lambda m: int(str(m)[1:-1]))
    ref = singles.groupby("pos")["r"].apply(lambda s: s.abs().mean())
    for p, val in ref.items():
        assert abs(posf.loc[posf.position == p, "e_i"].iloc[0] - val) < 1e-12, (
            f"e_i misaligned at {p}"
        )
    print("score_assay synthetic OK (incl. e_i alignment)")


test_parse()
test_features_and_zscore()
test_score_assay_synthetic()
print("ALL UNIT TESTS PASS")
