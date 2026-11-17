import pandas as pd

df = pd.read_csv("data/reference/DMS_substitutions.csv")
assert len(df) == 217, f"manifest rows {len(df)} != 217"

# R2: ESM2 max input 1022 residues; longer sequences would need windowing,
# which introduces a non-model choice into per-residue signals.
R_SEQLEN = df.seq_len <= 1022
# R3: residue-level tests need enough mutations per assay for stable
# per-residue error statistics and bootstrap resampling.
R_COUNT = df.DMS_total_number_mutants >= 500

sel = R_SEQLEN & R_COUNT
cols = [
    "DMS_index",
    "DMS_id",
    "UniProt_ID",
    "seq_len",
    "DMS_total_number_mutants",
    "DMS_number_single_mutants",
    "DMS_number_multiple_mutants",
    "includes_multiple_mutants",
    "coarse_selection_type",
]
out = df.loc[sel, cols].copy()
out["n_unique_positions_upper"] = (
    pd.concat(
        [df.loc[sel, "DMS_number_single_mutants"], df.loc[sel, "seq_len"]],
        axis=1,
    ).min(axis=1)
    + df.loc[sel, "DMS_number_multiple_mutants"]
)
out.sort_values("DMS_index").to_csv("data/assay_list_proposed.csv", index=False)

total_mut = int(out["DMS_total_number_mutants"].sum())
total_pos = int(out["n_unique_positions_upper"].sum())
summary = {
    "assays_total": len(df),
    "pass_seqlen_1022": int(R_SEQLEN.sum()),
    "pass_count_500": int(R_COUNT.sum()),
    "assays_pass": int(sel.sum()),
    "total_mutants_pass": total_mut,
    "total_forward_upper": total_pos,
    "multi_assays": int(df.loc[sel, "includes_multiple_mutants"].sum()),
    "seq_len_max": int(out["seq_len"].max()),
    "seq_len_min": int(out["seq_len"].min()),
}
pd.Series(summary).to_csv("data/assay_list_summary.csv", header=False)
for k, v in summary.items():
    print(f"{k}: {v}")
