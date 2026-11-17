"""Local synthetic test for the P3 non-additivity statistics layer (no GPU)."""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from pipeline.p3_test import assay_delta_nonadd


def main():
    tmp = Path(tempfile.mkdtemp())
    n_dbl, per_group = 40, 10
    rows = []
    # eps: small values for the first 30 doubles, large for the last 10
    eps = np.concatenate(
        [
            np.linspace(-0.2, 0.2, 30),
            np.array([3.0, -3.0, 2.5, -2.5, 2.0, -2.0, 1.5, -1.5, 1.2, -1.2]),
        ]
    )
    for j in range(n_dbl):
        pos = j + 1
        s1, s2 = f"A{pos}G", f"C{pos + 100}D"
        rows.append({"mutant": s1, "DMS_score": 0.0, "model_score": -1.0})
        rows.append({"mutant": s2, "DMS_score": 0.0, "model_score": -1.0})
        rows.append(
            {"mutant": f"{s1}:{s2}", "DMS_score": float(eps[j]), "model_score": -2.0}
        )
    pd.DataFrame(rows).to_csv(tmp / "mut.csv", index=False)
    ri_rows = [{"position": p, "R_i": 2.0} for p in range(1, n_dbl + 101)]
    pd.DataFrame(ri_rows).to_csv(tmp / "ri.csv", index=False)

    d, reason = assay_delta_nonadd("SYN", tmp / "mut.csv", tmp / "ri.csv")
    assert reason == "", reason
    assert d["n_dbl"] == 40 and d["k"] == 10
    # s1 = s2 = 0, model_norm constant -> r = z(DMS_double) = eps_z exactly,
    # so |r| ordering equals |eps_z| ordering and delta is exactly computable
    eps_z = (eps - eps.mean()) / eps.std(ddof=0)
    order = np.argsort(np.abs(eps_z), kind="stable")
    expected_high = np.abs(eps_z[order[-10:]]).mean()
    expected_low = np.abs(eps_z[order[:10]]).mean()
    assert abs(d["mean_absr_high"] - expected_high) < 1e-9, d
    assert abs(d["mean_absr_low"] - expected_low) < 1e-9, d
    assert d["delta_a"] > 0
    print(f"P3 synthetic: delta_a = {d['delta_a']:.4f} matches construction — OK")

    # eligibility: drop one single -> its double is excluded
    df = pd.DataFrame(rows)
    df = df[df["mutant"] != "A1G"]  # remove s1 of the first double
    df.to_csv(tmp / "mut2.csv", index=False)
    d2, reason2 = assay_delta_nonadd("SYN2", tmp / "mut2.csv", tmp / "ri.csv")
    assert reason2 == "" and d2["n_dbl"] == 39, (reason2, d2)
    print("P3 single-presence filter — OK")

    # degenerate eps: all doubles share one eps value -> excluded
    rows3 = []
    for j in range(40):
        pos = j + 1
        s1, s2 = f"A{pos}G", f"C{pos + 100}D"
        rows3.append({"mutant": s1, "DMS_score": 1.0, "model_score": -1.0})
        rows3.append({"mutant": s2, "DMS_score": 1.0, "model_score": -1.0})
        rows3.append({"mutant": f"{s1}:{s2}", "DMS_score": 2.0, "model_score": -2.0})
    pd.DataFrame(rows3).to_csv(tmp / "mut3.csv", index=False)
    d3, reason3 = assay_delta_nonadd("SYN3", tmp / "mut3.csv", tmp / "ri.csv")
    assert reason3.startswith("degenerate_eps_sd"), reason3
    print("P3 degenerate eps exclusion — OK")


main()
print("ALL P3 STATS TESTS PASS")
