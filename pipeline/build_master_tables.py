"""Build per-mutation master tables for all 191 frozen assays, plus a manifest.

Outputs (in --out, default /work/master_tables):
  singles.csv.gz   one row per single mutation
  doubles.csv.gz   one row per double mutation (exactly 2 substitutions)
  higher.csv.gz    one row per mutation with >=3 substitutions
  MANIFEST.json    row counts, column lists, SHA256 of each artifact

R_i is recomputed for every mutated position of every assay (model-internal
forward passes only); the recomputation cross-checks logp(wt) against the
stored positions table. Split membership columns come from the frozen
P1/P2/P3 formation/test lists.
"""

import argparse
import gzip
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/work")
from pipeline.p2_test import compute_ri  # noqa: E402
from pipeline.scoring_math import zscore  # noqa: E402

SPLIT_COLS = ["p1_set", "p2_set", "p3_set"]


def sha256_file(p: pathlib.Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def split_flags(dms_id: str, splits: dict) -> dict:
    out = {}
    for r, (fset, tset) in splits.items():
        out[f"{r}_set"] = (
            "test" if dms_id in tset else ("formation" if dms_id in fset else "-")
        )
    return out


def rounds_test(df: pd.DataFrame, kind: str) -> pd.Series:
    p1 = df["p1_set"] == "test"
    p2 = df["p2_set"] == "test"
    p3 = (df["p3_set"] == "test") & (kind != "single")
    tags = np.where(p1, "P1;", "") + np.where(p2, "P2;", "") + np.where(p3, "P3;", "")
    return pd.Series(tags, index=df.index).str.rstrip(";")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="/work/esm2_model")
    ap.add_argument(
        "--data-dir",
        default="/work/DMS_ProteinGym_substitutions/DMS_ProteinGym_substitutions",
    )
    ap.add_argument("--ref-csv", default="/work/reference/DMS_substitutions.csv")
    ap.add_argument("--results-a", default="/work/results_a")
    ap.add_argument("--p1-dir", default="/work/p1")
    ap.add_argument("--p2-dir", default="/work/p2")
    ap.add_argument("--p3-dir", default="/work/p3")
    ap.add_argument("--ri-all-dir", default="/work/ri_all")
    ap.add_argument("--out", default="/work/master_tables")
    args = ap.parse_args()
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ri_dir = pathlib.Path(args.ri_all_dir)
    ri_dir.mkdir(parents=True, exist_ok=True)

    from pipeline.esm_masked import MaskedMarginalEngine

    ref = pd.read_csv(args.ref_csv)
    seq_of = dict(zip(ref["DMS_id"], ref["target_seq"].astype(str).str.upper()))
    engine = MaskedMarginalEngine(args.model_dir)

    splits = {}
    for r, pdir in [("p1", args.p1_dir), ("p2", args.p2_dir), ("p3", args.p3_dir)]:
        fset = set(
            pd.read_csv(pathlib.Path(pdir) / f"{r.upper()}_FORMATION.csv")["DMS_id"]
        )
        tset = set(pd.read_csv(pathlib.Path(pdir) / f"{r.upper()}_TEST.csv")["DMS_id"])
        splits[r] = (fset, tset)

    lst = pd.read_csv("/work/reference/assay_list_proposed.csv")
    rows_s, rows_d, rows_h = [], [], []
    max_dev = 0.0
    for _, a in lst.sort_values("DMS_index").iterrows():
        did = str(a["DMS_id"])
        seq = seq_of[did]
        mut = pd.read_csv(pathlib.Path(args.results_a) / "mutants" / f"{did}.csv")
        pos = pd.read_csv(pathlib.Path(args.results_a) / "positions" / f"{did}.csv")
        feat = pos.set_index("position")
        ri_csv = ri_dir / f"{did}.csv"
        if not ri_csv.exists():
            _, dev = compute_ri(engine, did, ref, pathlib.Path(args.data_dir), ri_csv)
            max_dev = max(max_dev, dev)
        ri = pd.read_csv(ri_csv).rename(columns={"position": "pos"})
        ri_of = ri.set_index("pos")["R_i"]
        flags = split_flags(did, splits)

        mut = mut[mut["DMS_score"].notna() & mut["model_score"].notna()].copy()
        kind = (
            mut["mutant"]
            .astype(str)
            .str.count(":")
            .map(lambda c: "single" if c == 0 else ("double" if c == 1 else "higher"))
        )
        mut["kind"] = kind
        mut["DMS_id"] = did
        score_of = mut.set_index("mutant")["model_score"]

        # singles
        s = mut[mut["kind"] == "single"].copy()
        s["position"] = s["mutant"].map(lambda m: int(str(m)[1:-1]))
        s["wt_aa"] = s["position"].map(lambda p: seq[p - 1])
        s["mut_aa"] = s["mutant"].map(lambda m: str(m)[-1])
        s["R_i"] = s["position"].map(ri_of)
        s["llr_norm"] = s["model_score"] / s["R_i"]
        s["margin"] = s["position"].map(feat["margin"])
        s["logp_wt"] = s["position"].map(feat["logp_wt"])
        s["entropy"] = s["position"].map(feat["entropy"])
        s["top1_prob"] = s["position"].map(feat["top1_prob"])
        s["pos_norm"] = s["position"].map(feat["pos_norm"])
        s["e_i"] = s["position"].map(feat["e_i"])
        s["n_single_at_pos"] = s["position"].map(feat["n_single"])
        for c in SPLIT_COLS:
            s[c] = flags[c]
        s["kind"] = "single"
        rows_s.append(
            s[
                [
                    "DMS_id",
                    "mutant",
                    "position",
                    "wt_aa",
                    "mut_aa",
                    "DMS_score",
                    "model_score",
                    "r",
                    "R_i",
                    "llr_norm",
                    "margin",
                    "logp_wt",
                    "entropy",
                    "top1_prob",
                    "pos_norm",
                    "e_i",
                    "n_single_at_pos",
                ]
                + SPLIT_COLS
            ]
        )

        # doubles
        d = mut[mut["kind"] == "double"].copy()
        parts = d["mutant"].str.split(":", expand=True).reindex(columns=[0, 1])
        for j, col in [(0, "1"), (1, "2")]:
            d[f"pos{col}"] = parts[j].map(lambda m: int(m[1:-1]))
            d[f"wt{col}"] = d[f"pos{col}"].map(lambda p: seq[p - 1])
            d[f"mut{col}"] = parts[j].map(lambda m: m[-1])
            d[f"llr{col}"] = parts[j].map(score_of)
            d[f"R_i{col}"] = d[f"pos{col}"].map(ri_of)
            d[f"llr_norm{col}"] = d[f"llr{col}"] / d[f"R_i{col}"]
            d[f"margin{col}"] = d[f"pos{col}"].map(feat["margin"])
            d[f"entropy{col}"] = d[f"pos{col}"].map(feat["entropy"])
        # experimental non-additivity (raw scale, both singles must exist)
        dms_of = mut.set_index("mutant")["DMS_score"]
        singles_set = set(mut.loc[mut["kind"] == "single", "mutant"])
        has1 = parts[0].isin(singles_set) & parts[0].map(dms_of).notna()
        has2 = parts[1].isin(singles_set) & parts[1].map(dms_of).notna()
        d["eps"] = np.where(
            has1 & has2,
            d["DMS_score"] - parts[0].map(dms_of) - parts[1].map(dms_of),
            np.nan,
        )
        for c in SPLIT_COLS:
            d[c] = flags[c]
        rows_d.append(
            d[
                [
                    "DMS_id",
                    "mutant",
                    "pos1",
                    "wt1",
                    "mut1",
                    "pos2",
                    "wt2",
                    "mut2",
                    "DMS_score",
                    "model_score",
                    "r",
                    "llr1",
                    "llr2",
                    "R_i1",
                    "R_i2",
                    "llr_norm1",
                    "llr_norm2",
                    "eps",
                    "margin1",
                    "margin2",
                    "entropy1",
                    "entropy2",
                ]
                + SPLIT_COLS
            ]
        )

        # higher (>=3)
        h = mut[mut["kind"] == "higher"].copy()
        h["k"] = h["mutant"].astype(str).str.count(":") + 1
        h["positions"] = h["mutant"].map(
            lambda m: ";".join(p[1:-1] for p in str(m).split(":"))
        )
        h["wt_list"] = h["mutant"].map(
            lambda m: ";".join(seq[int(p[1:-1]) - 1] for p in str(m).split(":"))
        )
        h["mut_list"] = h["mutant"].map(
            lambda m: ";".join(p[-1] for p in str(m).split(":"))
        )
        h["llr_norm_total"] = h["mutant"].map(
            lambda m: sum(
                (float(score_of.get(p, np.nan)) / ri_of[int(p[1:-1])])
                if int(p[1:-1]) in ri_of.index
                else np.nan
                for p in str(m).split(":")
            )
        )
        for c in SPLIT_COLS:
            h[c] = flags[c]
        rows_h.append(
            h[
                [
                    "DMS_id",
                    "mutant",
                    "k",
                    "positions",
                    "wt_list",
                    "mut_list",
                    "DMS_score",
                    "model_score",
                    "r",
                    "llr_norm_total",
                ]
                + SPLIT_COLS
            ]
        )

    tables = {
        "singles": pd.concat(rows_s, ignore_index=True),
        "doubles": pd.concat(rows_d, ignore_index=True),
        "higher": pd.concat(rows_h, ignore_index=True),
    }
    manifest = {"ri_logpwt_max_dev": max_dev, "files": {}}
    for name, df in tables.items():
        df["rounds_test"] = rounds_test(df, name)
        path = out / f"{name}.csv.gz"
        df.to_csv(path, index=False, compression={"method": "gzip", "compresslevel": 9})
        manifest["files"][f"{name}.csv.gz"] = {
            "rows": int(len(df)),
            "columns": list(df.columns),
            "sha256": sha256_file(path),
        }
        print(name, len(df), "rows")
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {k: v["sha256"][:16] for k, v in manifest["files"].items()}, indent=2
        )
    )


if __name__ == "__main__":
    main()
