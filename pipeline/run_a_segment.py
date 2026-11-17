"""A-segment driver: run masked-marginals over the frozen 191-assay list with
per-assay checkpointing. Usage: python -m pipeline.run_a_segment --data-dir ... --out-dir ..."""

import argparse
import pathlib
import sys
import time

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline.esm_masked import MaskedMarginalEngine
from pipeline.assay_scoring import score_assay


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="/work/esm2_model")
    ap.add_argument("--data-dir", default="/work/DMS_ProteinGym_substitutions")
    ap.add_argument("--ref-csv", required=True)
    ap.add_argument("--list-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = pathlib.Path(args.out_dir)
    (out / "mutants").mkdir(parents=True, exist_ok=True)
    (out / "positions").mkdir(parents=True, exist_ok=True)

    ref = pd.read_csv(args.ref_csv)
    lst = pd.read_csv(args.list_csv)
    engine = MaskedMarginalEngine(args.model_dir)

    summary_path = out / "assay_summary.csv"
    done = set()
    if summary_path.exists():
        prev = pd.read_csv(summary_path)
        done = set(prev["DMS_id"])
        print(f"resume: {len(done)} assays already done", flush=True)

    t0 = time.time()
    infos = []
    if summary_path.exists():
        infos = pd.read_csv(summary_path).to_dict("records")
    for _, a in lst.sort_values("DMS_index").iterrows():
        did = str(a["DMS_id"])
        if did in done:
            continue
        f = pathlib.Path(args.data_dir) / (a["DMS_id"] + ".csv")
        if not f.exists():
            print(f"MISSING FILE {did}", flush=True)
            continue
        df = pd.read_csv(f)
        wt_seq = str(ref.loc[ref.DMS_id == did, "target_seq"].iloc[0])
        mut, posf, info = score_assay(
            engine,
            did,
            wt_seq,
            df,
        )
        info["seconds"] = round(time.time() - t0, 1)
        mut.to_csv(out / "mutants" / f"{did}.csv", index=False)
        posf.to_csv(out / "positions" / f"{did}.csv", index=False)
        infos.append(info)
        pd.DataFrame(infos).to_csv(summary_path, index=False)
        print(
            f"{did} rho={info['rho_all']:.4f} n={info['n_valid']} t={info['seconds']}",
            flush=True,
        )
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
