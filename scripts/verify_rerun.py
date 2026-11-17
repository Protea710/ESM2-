"""Compare the rerun outputs against the archived first-run outputs.

Structural mismatches (missing files, row-count or DMS-set differences) exit
non-zero. Numeric deviations are REPORTED in full and do not gate: the
round report must state them either way.
"""

import json
import pathlib
import sys

import numpy as np
import pandas as pd

ARCHIVE = pathlib.Path("/work/run1_archive")
RERUN = {
    "assay_summary.csv": "/work/results_a/assay_summary.csv",
    "compare.csv": "/work/results_a/compare.csv",
    "P1_per_assay.csv": "/work/results_p1/P1_per_assay.csv",
    "P2_per_assay.csv": "/work/results_p2/P2_per_assay.csv",
    "P3_per_assay.csv": "/work/results_p3/P3_per_assay.csv",
}
JSONS = {
    "P1_summary.json": "/work/results_p1/P1_summary.json",
    "P2_summary.json": "/work/results_p2/P2_summary.json",
    "P3_summary.json": "/work/results_p3/P3_summary.json",
}
KEYS = {
    "P1_summary.json": ["theta", "ci_low", "ci_high"],
    "P2_summary.json": ["theta_norm", "ci_low", "ci_high"],
    "P3_summary.json": ["theta_3", "ci_low", "ci_high"],
}


def cmp_csv(name: str, old_path: pathlib.Path, new_path: str):
    old, new = pd.read_csv(old_path), pd.read_csv(new_path)
    if list(old.columns) != list(new.columns):
        print(f"[{name}] COLUMN MISMATCH")
        return
    if len(old) != len(new) or not old.iloc[:, 0].equals(new.iloc[:, 0]):
        print(f"[{name}] STRUCTURAL MISMATCH rows {len(old)} vs {len(new)}")
        sys.exit(1)
    byte_equal = old_path.read_bytes() == pathlib.Path(new_path).read_bytes()
    devs = {}
    for c in old.columns:
        if old[c].dtype.kind in "fi" and new[c].dtype.kind in "fi":
            d = float(
                np.nanmax(np.abs(old[c].to_numpy(float) - new[c].to_numpy(float)))
            )
            if d:
                devs[c] = d
    print(
        f"[{name}] byte_equal={byte_equal} rows={len(old)} "
        + (f"max_abs_dev={devs}" if devs else "numeric identical")
    )


def main():
    for name, new_path in RERUN.items():
        cmp_csv(name, ARCHIVE / name, new_path)
    for name, new_path in JSONS.items():
        old = json.loads((ARCHIVE / name).read_text())
        new = json.loads(pathlib.Path(new_path).read_text())
        keys = KEYS[name]
        diffs = {k: (old.get(k), new.get(k)) for k in keys if old.get(k) != new.get(k)}
        if diffs:
            print(f"[{name}] deviations: {diffs}")
        else:
            print(
                f"[{name}] theta/ci bit-identical: "
                f"theta={new[keys[0]]:.10f} ci=[{new[keys[1]]:.10f},{new[keys[2]]:.10f}] "
                f"verdict={new.get('verdict')}"
            )


main()
print("RERUN COMPARISON COMPLETE")
