"""Local synthetic test for the P1 statistics layer (no GPU, no real data)."""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from pipeline.p1_test import assay_delta


def main():
    tmp = Path(tempfile.mkdtemp())
    # 12 eligible positions: margin ascending with position; error constructed
    # so low-margin positions have e_i = 1.0 and high-margin have e_i = 0.5.
    n_pos, per = 12, 10
    rows_m, rows_p = [], []
    for j in range(n_pos):
        pos = j + 1
        e = 1.0 if j < 6 else 0.5
        for m in range(per):
            # residuals: low-margin group |r|=1.0, high-margin group |r|=0.5
            val = e if m % 2 == 0 else -e
            rows_m.append(
                {"mutant": f"A{pos}G", "DMS_score": float(val), "model_score": 0.0}
            )
        rows_p.append({"position": pos, "margin": 0.01 * j})
    pd.DataFrame(rows_m).to_csv(tmp / "mut.csv", index=False)
    pd.DataFrame(rows_p).to_csv(tmp / "pos.csv", index=False)
    d, reason = assay_delta("SYN", tmp / "mut.csv", tmp / "pos.csv")
    assert reason == "", reason
    assert d["n_eligible"] == 12 and d["k"] == 3
    # z is over ALL singles in the assay: pooled sd = sqrt(0.625), so the
    # constructed e values (1.0 low / 0.5 high) scale to these exact values
    sd = np.sqrt(0.625)
    assert abs(d["mean_e_low"] - 1.0 / sd) < 1e-12
    assert abs(d["mean_e_high"] - 0.5 / sd) < 1e-12
    assert abs(d["delta_a"] - 0.5 / sd) < 1e-12
    print(f"P1 synthetic: delta_a = {d['delta_a']:.4f} matches constructed effect — OK")

    # eligibility: a position with 7 singles must be excluded
    rows_m2 = rows_m + [
        {"mutant": "A13G", "DMS_score": 1.0, "model_score": 0.0} for _ in range(7)
    ]
    rows_p2 = rows_p + [{"position": 13, "margin": 0.5}]
    pd.DataFrame(rows_m2).to_csv(tmp / "mut2.csv", index=False)
    pd.DataFrame(rows_p2).to_csv(tmp / "pos2.csv", index=False)
    d2, reason2 = assay_delta("SYN2", tmp / "mut2.csv", tmp / "pos2.csv")
    assert reason2 == "" and d2["n_eligible"] == 12, (reason2, d2)
    print("P1 eligibility (>=8 singles) enforced — OK")


main()
print("ALL P1 STATS TESTS PASS")
