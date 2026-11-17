"""Compare per-assay Spearman against the official ProteinGym ESM2-650M column."""

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--official", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ours = pd.read_csv(args.summary)
    off = pd.read_csv(args.official)
    col = [c for c in off.columns if "650M" in c][0]
    o = off[["DMS ID", col]].rename(columns={"DMS ID": "DMS_id", col: "rho_official"})
    m = ours.merge(o, on="DMS_id", how="left")
    assert m["rho_official"].notna().all(), "unmatched DMS_id against official CSV"
    m["delta"] = m["rho_all"] - m["rho_official"]
    m[["DMS_id", "rho_all", "rho_official", "delta", "n_valid"]].to_csv(
        args.out, index=False
    )
    stats = {
        "n_assays": len(m),
        "mean_rho_ours": m["rho_all"].mean(),
        "mean_rho_official": m["rho_official"].mean(),
        "mean_abs_delta": m["delta"].abs().mean(),
        "max_abs_delta": m["delta"].abs().max(),
        "n_gt_005": int((m["delta"].abs() > 0.05).sum()),
        "n_gt_002": int((m["delta"].abs() > 0.02).sum()),
    }
    for k, v in stats.items():
        print(f"{k}: {v if not isinstance(v, float) else round(v, 5)}")


if __name__ == "__main__":
    main()
