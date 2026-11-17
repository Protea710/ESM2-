"""Build P2 formation/test lists.

Test-set pool is restricted to P1's formation assays (126) — those assays have
never had any per-position or per-mutation residual statistics computed. P1's
test assays (65) go to P2 formation and are declared as seen data in the
proposition. Only frozen list metadata is read; no result files are touched.
"""

import hashlib
from pathlib import Path

import pandas as pd

SEED = "P2-split-v1-2026-11-17"
POOL = Path("results/p1/P1_FORMATION.csv")
OUT = Path("results/p2")
KEEP = [
    "DMS_index",
    "DMS_id",
    "UniProt_ID",
    "coarse_selection_type",
    "seq_len",
    "DMS_total_number_mutants",
    "DMS_number_single_mutants",
    "DMS_number_multiple_mutants",
]


def stable_key(text: str) -> str:
    return hashlib.sha256(f"{SEED}|{text}".encode()).hexdigest()


def main():
    pool = pd.read_csv(POOL)
    assert len(pool) == 126, f"pool rows {len(pool)} != 126"
    assert pool["DMS_id"].is_unique

    groups = []
    for uniprot, g in pool.groupby("UniProt_ID", sort=True):
        stratum = "+".join(sorted(g["coarse_selection_type"].fillna("NA").unique()))
        groups.append(
            {
                "UniProt_ID": uniprot,
                "stratum": stratum,
                "n_assays": len(g),
                "key": stable_key(uniprot),
            }
        )
    groups = pd.DataFrame(groups)

    test_proteins = set()
    for _, sg in groups.groupby("stratum", sort=True):
        sg = sg.sort_values("key").reset_index(drop=True)
        target = max(1, round(sg["n_assays"].sum() / 3))
        chosen = 0
        for _, row in sg.iterrows():
            if chosen >= target:
                break
            test_proteins.add(row["UniProt_ID"])
            chosen += int(row["n_assays"])

    is_test = pool["UniProt_ID"].isin(list(test_proteins))
    test = pool.loc[is_test, KEEP].sort_values("DMS_index")
    formation = pool.loc[~is_test, KEEP].sort_values("DMS_index")

    assert set(formation["DMS_id"]).isdisjoint(set(test["DMS_id"]))
    assert set(formation["UniProt_ID"]).isdisjoint(set(test["UniProt_ID"]))
    assert len(formation) + len(test) == 126

    OUT.mkdir(parents=True, exist_ok=True)
    formation.to_csv(OUT / "P2_FORMATION.csv", index=False)
    test.to_csv(OUT / "P2_TEST.csv", index=False)
    print(f"formation_assays: {len(formation)}")
    print(f"test_assays: {len(test)}")
    print(f"formation_proteins: {formation['UniProt_ID'].nunique()}")
    print(f"test_proteins: {test['UniProt_ID'].nunique()}")
    print("test_selection_counts:")
    print(
        test["coarse_selection_type"]
        .value_counts(dropna=False)
        .sort_index()
        .to_string()
    )


if __name__ == "__main__":
    main()
