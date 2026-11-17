"""Build P1 formation/test lists from frozen assay metadata only.

No result, residual, score, or per-position file is read. Assays sharing a
UniProt_ID are kept together. coarse_selection_type is used only to stratify.
"""

import hashlib
from pathlib import Path

import pandas as pd

SEED = "P1-split-v1-2026-11-17"
SOURCE = Path("data/assay_list_proposed.csv")
OUT = Path("results/p1")
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
    assays = pd.read_csv(SOURCE)
    assert len(assays) == 191
    assert assays["DMS_id"].is_unique

    # A protein may have assays under multiple selection types. The sorted tuple
    # is its stratum; the protein remains indivisible across the split.
    groups = []
    for uniprot, g in assays.groupby("UniProt_ID", sort=True):
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

    is_test = assays["UniProt_ID"].isin(test_proteins)
    formation = assays.loc[~is_test, KEEP].sort_values("DMS_index")
    test = assays.loc[is_test, KEEP].sort_values("DMS_index")

    assert set(formation["DMS_id"]).isdisjoint(test["DMS_id"])
    assert set(formation["UniProt_ID"]).isdisjoint(test["UniProt_ID"])
    assert len(formation) + len(test) == 191

    OUT.mkdir(parents=True, exist_ok=True)
    formation.to_csv(OUT / "P1_FORMATION.csv", index=False)
    test.to_csv(OUT / "P1_TEST.csv", index=False)
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
