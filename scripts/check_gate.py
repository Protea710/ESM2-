import sys

import pandas as pd

failures = []

df = pd.read_csv("data/reference/DMS_substitutions.csv")
if len(df) != 217:
    failures.append(f"manifest rows {len(df)} != 217")

lst = pd.read_csv("data/assay_list_proposed.csv")
if len(lst) != 191:
    failures.append(f"proposed list rows {len(lst)} != 191")
if lst.seq_len.max() > 1022:
    failures.append("R2 violated: seq_len > 1022 in proposed list")
if (lst.DMS_total_number_mutants < 500).any():
    failures.append("R3 violated: mutant count < 500 in proposed list")
if lst.DMS_id.duplicated().any():
    failures.append("duplicated DMS_id in proposed list")

summary = pd.read_csv("data/assay_list_summary.csv", header=None, index_col=0)[1]
expected = {
    "assays_total": "217",
    "assays_pass": "191",
    "total_mutants_pass": "2411093",
    "total_forward_upper": "1822479",
    "multi_assays": "69",
}
for k, v in expected.items():
    if str(summary[k]) != v:
        failures.append(f"summary {k}={summary[k]} != {v}")

gate = open("docs/GATE1.md", encoding="utf-8").read()
for section in [
    "测定清单及筛选规则",
    "一轮命题的格式",
    "「误差」的定义",
    "通道验证",
    "三分类判定规则",
]:
    if section not in gate:
        failures.append(f"GATE1 missing section: {section}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print(
    f"gate integrity OK: {len(lst)} assays, filters satisfied, GATE1 sections present"
)
