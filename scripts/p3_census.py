"""Label-only census of exactly-two-substitution mutants whose both component
singles are present, over the frozen 191-assay list. Reads mutant strings only
(never DMS_score); runs on the full list BEFORE any P3 split exists."""

import pathlib
import sys

import pandas as pd

sys.path.insert(0, "/work")


def main():
    lst = pd.read_csv("/work/reference/assay_list_proposed.csv")
    rows = []
    for did in lst["DMS_id"]:
        mut = pd.read_csv(
            pathlib.Path(
                "/work/DMS_ProteinGym_substitutions/DMS_ProteinGym_substitutions"
            )
            / f"{did}.csv",
            usecols=["mutant"],
        )["mutant"].astype(str)
        singles = set(mut[~mut.str.contains(":")])
        doubles = mut[mut.str.count(":") == 1]
        n_both = 0
        for m in doubles:
            a, b = m.split(":")
            if a in singles and b in singles:
                n_both += 1
        rows.append(
            {"DMS_id": did, "n_doubles": len(doubles), "n_doubles_both_singles": n_both}
        )
    out = pd.DataFrame(rows)
    out.to_csv("/work/reference/P3_census.csv", index=False)
    print(out["n_doubles_both_singles"].describe().to_string())
    print(
        "assays with >=32 doubles-both-singles:",
        int((out["n_doubles_both_singles"] >= 32).sum()),
    )
    print("total doubles-both-singles:", int(out["n_doubles_both_singles"].sum()))


if __name__ == "__main__":
    main()
