"""Audit: every number drawn on the plate must be reproduced from results/ at build time."""

import json
import os
import sys

import pandas as pd

sys.path.insert(0, "figures")
import build_plate as B
import vizlib as V


def counts():
    s = pd.read_csv(
        "results/master_tables/singles.csv.gz",
        usecols=["DMS_id", "n_single_at_pos", "p1_set", "p2_set"],
    )
    out = {}
    for tag, r in (("P1", "p1"), ("P2", "p2")):
        per = pd.read_csv(f"results/{r}/{tag}_per_assay.csv")
        m = (
            (s[f"{r}_set"] == "test")
            & s.DMS_id.isin(set(per.DMS_id))
            & (s.n_single_at_pos >= 8)
        )
        out[tag] = int(m.sum())
    p3 = pd.read_csv("results/p3/P3_per_assay.csv")
    out["P3"] = int(p3.n_dbl.sum())
    return out


def expected():
    man = json.loads(open("results/master_tables/MANIFEST.json").read())
    cmp = pd.read_csv("results/compare.csv")
    rho = pd.read_csv("results/assay_summary.csv").set_index("DMS_id")["rho_single"]
    p1 = pd.read_csv("results/p1/P1_per_assay.csv")
    p3 = pd.read_csv("results/p3/P3_per_assay.csv")
    n = counts()

    exp = {
        "总突变数": f"{man['total_rows_all_tables']:,}",
        "P1 行数": f"{n['P1']:,}",
        "P2 行数": f"{n['P2']:,}",
        "P3 行数": f"{n['P3']:,}",
        "A 段均值": f"{cmp.rho_all.mean():.5f}",
        "官方均值": f"{cmp.rho_official.mean():.5f}",
        "最大绝对差": f"{abs(cmp.delta).max():.6f}",
        "P1 方向": f"{int((p1.delta_a < 0).sum())}/{len(p1)}",
        "P3 方向": f"{int((p3.delta_a > 0).sum())} 正 / "
        f"{int((p3.delta_a < 0).sum())} 负",
        "SUMO1 ρ": f"{rho['SUMO1_HUMAN_Weile_2017']:.5f}",
        "TNKS2 ρ": f"{rho['TNKS2_HUMAN_Tsuboyama_2023_5JRT']:.5f}",
        "MLAC ρ": f"{rho['MLAC_ECOLI_MacRae_2023']:.5f}",
    }
    for tag, lo, hi in (
        ("P1", "mean_e_low", "mean_e_high"),
        ("P2", "mean_e_low", "mean_e_high"),
        ("P3", "mean_absr_low", "mean_absr_high"),
    ):
        p = pd.read_csv(f"results/{tag.lower()}/{tag}_per_assay.csv")
        summ = json.loads(open(f"results/{tag.lower()}/{tag}_summary.json").read())
        theta = summ.get("theta", summ.get("theta_norm", summ.get("theta_3")))
        lvl = float(((p[lo] + p[hi]) / 2).median())
        exp[f"{tag} 残差量级"] = f"{lvl:.2f}"
        exp[f"{tag} 比值"] = f"{lvl / abs(theta):.1f}"
    for tag, e in B.event_times().items():
        exp[f"{tag} 命题时刻"] = f"{e['prop']:%H:%M:%S}"
        exp[f"{tag} 哈希时刻"] = f"{e['freeze']:%H:%M:%S}"
        exp[f"{tag} 统计时刻"] = f"{e['stats']:%H:%M:%S}"
    return exp


def main():
    V.apply_style()
    fig = V.mm_figure(B.PAGE_W, B.PAGE_H)
    fig.text(
        B.MARGIN / B.PAGE_W,
        1 - 14 / B.PAGE_H,
        "模型知道自己错在哪吗：三轮预登记检验的执行链",
        fontsize=18,
        va="top",
    )
    fig.text(
        B.MARGIN / B.PAGE_W,
        1 - 25 / B.PAGE_H,
        "全部数字只读自 results/ 与 docs/，逐项出处见复核说明",
        fontsize=10,
        va="top",
        color=V.PAL["dai"],
    )
    B.panel_a(fig, V.mm_axes(fig, B.X_L, B.Y_TOP, B.COL_L, B.ROW_H))
    B.panel_b(fig)
    B.panel_c(fig)
    B.panel_d(fig)

    texts = "\n".join(s for s, _ in V.text_boxes(fig))

    bad = []
    for label, val in expected().items():
        ok = val in texts
        print(f"{'OK ' if ok else 'MISS'} {label:12s} {val}")
        if not ok:
            bad.append((label, val))
    print(f"\n{len(bad)} missing of {len(expected())}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
