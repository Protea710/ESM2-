"""Build P3 formation/test lists.

P3 test set = all census-eligible assays (>=32 doubles with both component
singles present, label-only census over the full frozen list before any split
existed) that lie inside the never-touched pool (P2 formation, 78 assays —
those have had no per-position or per-mutation residual statistics).
Formation = every other frozen assay. Only frozen metadata is read.
"""

from pathlib import Path

import pandas as pd

CENSUS = Path("data/reference/P3_census.csv")
CLEAN_POOL = Path("results/p2/P2_FORMATION.csv")
LIST = Path("data/assay_list_proposed.csv")
OUT = Path("results/p3")
MIN_DOUBLES = 32
KEEP = [
    "DMS_index",
    "DMS_id",
    "UniProt_ID",
    "coarse_selection_type",
    "seq_len",
    "DMS_total_number_mutants",
    "DMS_number_single_mutants",
    "DMS_number_multiple_mutants",
    "n_doubles_both_singles",
]


def main():
    cen = pd.read_csv(CENSUS)
    lst = pd.read_csv(LIST)
    pool = pd.read_csv(CLEAN_POOL)
    assert len(cen) == 191 and len(lst) == 191 and len(pool) == 78
    m = lst.merge(cen, on="DMS_id", validate="1:1")
    clean = set(pool["DMS_id"])

    test = m[(m["n_doubles_both_singles"] >= MIN_DOUBLES) & m["DMS_id"].isin(clean)]
    formation = m[~m["DMS_id"].isin(set(test["DMS_id"]))]

    assert set(test["DMS_id"]).isdisjoint(set(formation["DMS_id"]))
    assert len(test) + len(formation) == 191
    assert set(test["DMS_id"]) <= clean, "test must stay inside the never-touched pool"

    OUT.mkdir(parents=True, exist_ok=True)
    test[KEEP].sort_values("DMS_index").to_csv(OUT / "P3_TEST.csv", index=False)
    formation[KEEP].sort_values("DMS_index").to_csv(
        OUT / "P3_FORMATION.csv", index=False
    )
    print(
        f"test_assays: {len(test)}; doubles_both_singles: {int(test['n_doubles_both_singles'].sum())}"
    )
    print(f"formation_assays: {len(formation)}")


if __name__ == "__main__":
    main()
