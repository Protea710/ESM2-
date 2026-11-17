"""Independent check: figure A's numbers against the frozen split files."""

import sys

import pandas as pd

sys.path.insert(0, "figures")
import build_alloc as B


def main():
    df, splits = B.load()
    df = B.order(df)
    rows = {}
    for r in ("P1", "P2", "P3"):
        rows[r] = {k: int(v) for k, v in df[r].value_counts().items()}

    p1t = set(pd.read_csv("results/p1/P1_TEST.csv").DMS_id)
    p2t = set(pd.read_csv("results/p2/P2_TEST.csv").DMS_id)
    p3t = set(pd.read_csv("results/p3/P3_TEST.csv").DMS_id)
    lm = pd.read_csv("data/assay_list_proposed.csv")
    d = pd.read_csv("results/master_tables/doubles.csv.gz", usecols=["DMS_id", "eps"])
    elig = set(d.groupby("DMS_id")["eps"].count()[lambda s: s >= 32].index)

    checks = [
        ("P1 检验 65", rows["P1"].get("test") == len(p1t) == 65),
        ("P1 成题 126", rows["P1"].get("form") == 126),
        ("P2 检验 48", rows["P2"].get("test") == len(p2t) == 48),
        ("P2 成题 78", rows["P2"].get("form") == 78),
        ("P2 已见 65 = P1 检验", rows["P2"].get("seen") == len(p1t)),
        ("P3 检验 26", rows["P3"].get("test") == len(p3t) == 26),
        ("P3 无数据 124 = 191-67", rows["P3"].get("none") == 191 - len(elig) == 124),
        ("P3 已见 41 = |合格 ∩ 往轮检验|", rows["P3"].get("seen") == len(elig & (p1t | p2t)) == 41),
        ("检验集两两交集 0", not (p1t & p2t) and not (p1t & p3t) and not (p2t & p3t)),
        ("191 测定 / 164 蛋白", len(lm) == 191 and lm.UniProt_ID.nunique() == 164),
    ]
    g = lm.groupby("UniProt_ID").size()
    checks.append(("20 个多测定蛋白 / 47 测定", int((g > 1).sum()) == 20 and int(g[g > 1].sum()) == 47))
    multi = set(g[g > 1].index)
    checks.append(("P3 检验集不含多测定蛋白",
                   len(p3t & set(lm[lm.UniProt_ID.isin(multi)].DMS_id)) == 0))

    bad = 0
    for name, ok in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
        bad += 0 if ok else 1
    print(f"\n{len(checks) - bad}/{len(checks)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
