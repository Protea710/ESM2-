"""Local synthetic test for the P2 statistics layer (no GPU, no real data)."""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from pipeline.p2_test import assay_delta_norm


def main():
    tmp = Path(tempfile.mkdtemp())
    # 12 positions, margins ascending; LLR amplitude large on high-margin half.
    # Normalization by R_i removes the amplitude difference: after dividing by
    # R_i, low- and high-margin groups have identical |r| -> delta_a = 0.
    rows_m, rows_p, rows_ri = [], [], []
    for j in range(12):
        pos = j + 1
        llr = 4.0 if j < 6 else 1.0  # high-margin positions: big |LLR|
        ri = llr / 1.0  # choose R_i = |LLR| so normalized score = -1 for all
        for m in range(10):
            dms = 1.0 if m % 2 == 0 else -1.0
            rows_m.append({"mutant": f"A{pos}G", "DMS_score": dms, "model_score": -llr})
        rows_p.append({"position": pos, "margin": 0.01 * j})
        rows_ri.append({"position": pos, "R_i": ri})
    pd.DataFrame(rows_m).to_csv(tmp / "mut.csv", index=False)
    pd.DataFrame(rows_p).to_csv(tmp / "pos.csv", index=False)
    pd.DataFrame(rows_ri).to_csv(tmp / "ri.csv", index=False)
    d, reason = assay_delta_norm(
        "SYN", tmp / "mut.csv", tmp / "pos.csv", tmp / "ri.csv"
    )
    assert reason == "", reason
    # all normalized scores are -1 -> z_model identical for all -> |r| constant
    assert abs(d["mean_e_low"] - d["mean_e_high"]) < 1e-9, d
    assert abs(d["delta_a"]) < 1e-9, d
    print(
        f"P2 synthetic: normalization flattens amplitude effect, delta_a={d['delta_a']:.2e} — OK"
    )

    # degenerate R_i position excluded
    rows_ri2 = rows_ri + [{"position": 13, "R_i": 0.0}]
    for m in range(10):
        rows_m.append({"mutant": "A13G", "DMS_score": 1.0, "model_score": -2.0})
    pd.DataFrame(rows_m).to_csv(tmp / "mut2.csv", index=False)
    pd.DataFrame(rows_ri2).to_csv(tmp / "ri2.csv", index=False)
    d2, reason2 = assay_delta_norm(
        "SYN2", tmp / "mut2.csv", tmp / "pos.csv", tmp / "ri2.csv"
    )
    assert reason2 == "" and d2["n_eligible"] == 12, (reason2, d2)
    print("P2 degenerate R_i exclusion — OK")


main()
print("ALL P2 STATS TESTS PASS")
