"""Figure A: where each of the 191 assays went, round by round, plus the sealed window."""

import datetime as dt
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "figures")
import vizlib as V

PAGE_W, PAGE_H = 594, 420
M = 18
X0, X1 = 62, 470
ROW_Y = {"P1": 322, "P2": 274, "P3": 226}
BAR_H = 44
LOG_LO, LOG_HI = 1.2, 5.8

ST = {
    "test": ("检验", V.PAL["carmine"]),
    "form": ("成题", V.PAL["indigo"]),
    "seen": ("已见", V.PAL["seen"]),
    "none": ("无该轮数据", None),
}


def load():
    lst = pd.read_csv("data/assay_list_proposed.csv")[["DMS_id", "UniProt_ID"]]
    s = pd.read_csv("results/master_tables/singles.csv.gz", usecols=["DMS_id"])
    d = pd.read_csv("results/master_tables/doubles.csv.gz", usecols=["DMS_id", "eps"])
    df = lst.merge(
        s.groupby("DMS_id").size().rename("n_single"), on="DMS_id", how="left"
    )
    df = df.merge(
        d.groupby("DMS_id")["eps"].count().rename("n_dbl"), on="DMS_id", how="left"
    )
    df = df.fillna({"n_single": 0, "n_dbl": 0})

    splits = {}
    for r in ("P1", "P2", "P3"):
        splits[r] = (
            set(pd.read_csv(f"results/{r.lower()}/{r}_FORMATION.csv").DMS_id),
            set(pd.read_csv(f"results/{r.lower()}/{r}_TEST.csv").DMS_id),
        )

    nd = dict(zip(df.DMS_id, df.n_dbl))

    def state(dms, r):
        form, test = splits[r]
        if r == "P3" and nd[dms] < 32:
            return "none"
        if dms in test:
            return "test"
        earlier = set()
        for q in ("P1", "P2"):
            if q < r:
                earlier |= splits[q][1]
        if dms in earlier:
            return "seen"
        return "form" if dms in form else "none"

    for r in ("P1", "P2", "P3"):
        df[r] = [state(x, r) for x in df.DMS_id]
    return df, splits


def order(df):
    return df.sort_values(
        ["UniProt_ID", "n_single"], ascending=[True, False]
    ).reset_index(drop=True)


def ypos(v, base):
    if v <= 0:
        return base
    t = (np.log10(v) - LOG_LO) / (LOG_HI - LOG_LO)
    return base + BAR_H * max(t, 0.02)


def strip(fig, df):
    ax = V.mm_axes(fig, 0, 190, PAGE_W, 220)
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0, 220)
    n = len(df)
    pitch = (X1 - X0) / n
    bw = pitch * 0.82

    for tag in ("P1", "P2", "P3"):
        base = ROW_Y[tag] - 190
        col = "n_dbl" if tag == "P3" else "n_single"
        for i, row in df.iterrows():
            x = X0 + i * pitch
            state = row[tag]
            if state == "none":
                ax.plot([x, x + bw], [base, base], color=V.PAL["grey"], lw=0.35,
                        zorder=1)
                continue
            h = ypos(row[col], base) - base
            if state == "form":
                ax.add_patch(plt_rect(x, base, bw, h, facecolor=V.PAL["paper"],
                                      edgecolor=V.PAL["indigo"], lw=0.45))
            else:
                ax.add_patch(plt_rect(x, base, bw, h, facecolor=ST[state][1],
                                      edgecolor="none"))

    for tag in ("P1", "P2", "P3"):
        base = ROW_Y[tag] - 190
        ax.plot([X0 - 1, X1 + 1], [base, base], color=V.PAL["ink"], lw=0.5)
        ax.plot([X0 - 1, X0 - 1], [base, base + BAR_H], color=V.PAL["ink"], lw=0.5)
        ax.text(18, base + BAR_H + 1, tag, fontsize=13, va="bottom",
                color=V.PAL["indigo"])
        ax.text(18, base + BAR_H - 5, "单点" if tag != "P3" else "双点", fontsize=10,
                va="bottom", color=V.PAL["dai"])

    for tag in ("P1", "P2", "P3"):
        base = ROW_Y[tag] - 190
        for v in (100, 1000, 10000, 100000):
            y = ypos(v, base)
            ax.plot([X0 - 2.6, X0 - 1], [y, y], color=V.PAL["grey"], lw=0.6)
            if tag == "P3":
                ax.text(X0 - 4, y, f"{v:,}", fontsize=10, ha="right", va="center",
                        color=V.PAL["dai"])


def plt_rect(x, y, w, h, **kw):
    from matplotlib.patches import Rectangle

    return Rectangle((x, y), w, h, **kw)


def brackets(fig, df):
    ax = V.mm_axes(fig, 0, 190, PAGE_W, 220)
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0, 220)
    n = len(df)
    pitch = (X1 - X0) / n
    y = 368 - 190
    groups = df.groupby("UniProt_ID").size()
    multi = groups[groups > 1]
    for idx, prot in enumerate(multi.index, 1):
        rows = df.index[df.UniProt_ID == prot]
        x0 = X0 + rows[0] * pitch
        x1 = X0 + (rows[-1] + 1) * pitch
        ax.plot([x0, x0, x1, x1], [y, y + 3.5, y + 3.5, y], color=V.PAL["dai"], lw=0.6)
        ax.text((x0 + x1) / 2, y + 4.5, str(idx), fontsize=10, ha="center",
                va="bottom", color=V.PAL["dai"])
    return list(multi.index)


def groups(fig, df):
    ax = V.mm_axes(fig, 0, 0, PAGE_W, 130)
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0, 130)
    ax.text(M, 100, "20 个多测定蛋白（共 47 个测定）在各轮的整组去向",
            fontsize=11, va="top")
    served = {}
    for prot, sub in df.groupby("UniProt_ID"):
        rounds = [r for r in ("P1", "P2", "P3") if (sub[r] == "test").any()]
        served[prot] = rounds[0] if rounds else "未检验"
    multi = df.groupby("UniProt_ID").size()
    multi = multi[multi > 1]
    for i, prot in enumerate(multi.index, 1):
        col = (i - 1) // 7
        row = (i - 1) % 7
        x = M + col * 186
        y = 88 - row * 9.0
        ax.text(x, y, str(i), fontsize=10, va="center", color=V.PAL["grey"])
        ax.text(x + 7, y, f"{prot} ×{multi[prot]}", fontsize=10, va="center")
        ax.text(x + 116, y, served[prot], fontsize=10, va="center",
                color=V.PAL["carmine"] if served[prot] != "未检验" else V.PAL["dai"])


def legend(fig, df):
    ax = V.mm_axes(fig, 0, 190, PAGE_W, 220)
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0, 220)
    x, y = 490, 208
    ax.text(x, y, "读法", fontsize=11, va="bottom")
    items = [("test", "检验：该轮统计读了它"),
             ("form", "成题：在形成集，未读"),
             ("seen", "已见：往轮检验过，本轮不用")]
    yy = y - 5
    for key, label in items:
        if key == "form":
            ax.add_patch(plt_rect(x, yy - 1.4, 5, 4.6, facecolor=V.PAL["paper"],
                                  edgecolor=V.PAL["indigo"], lw=0.6))
        else:
            ax.add_patch(plt_rect(x, yy - 1.4, 5, 4.6, facecolor=ST[key][1],
                                  edgecolor="none"))
        ax.text(x + 7.5, yy + 0.9, label, fontsize=10, va="center")
        yy -= 8
    ax.plot([x, x + 5], [yy + 0.9, yy + 0.9], color=V.PAL["grey"], lw=0.5)
    ax.text(x + 7.5, yy + 0.9, "无该轮数据：<32 个双点\n（P3 口径），或不在该轮池内",
            fontsize=10, va="top", linespacing=1.5)

    ty = 158
    ax.text(x, ty + 6, "各轮去向", fontsize=11, va="bottom")
    heads = ["轮", "检验", "成题", "已见", "无数据"]
    cols = [0, 9, 24, 39, 54]
    for c, h in zip(cols, heads):
        ax.text(x + c, ty, h, fontsize=10, va="bottom", color=V.PAL["dai"])
    yy = ty - 2
    for tag in ("P1", "P2", "P3"):
        c = df[tag].value_counts()
        vals = [tag, c.get("test", 0), c.get("form", 0), c.get("seen", 0),
                c.get("none", 0)]
        for cc, v in zip(cols, vals):
            ax.text(x + cc, yy, str(v), fontsize=10, va="top")
        yy -= 6

    ax.text(x, yy - 10,
            "检验集互斥：65 / 48 / 26，\n两两交集 0。20 个多测定蛋白\n（47 个测定）整组同侧，\n无一组在 P3 当检验集。",
            fontsize=10, va="top", linespacing=1.6)
    ax.text(x, 78, "两处执行修正\n（未改统计口径）", fontsize=11, va="top",
            linespacing=1.5)
    ax.text(x, 62, "① P3 首跑 CSV 误读：改读 A 段\n    产物后全量重跑。\n② 门禁误删：模块清单恢复后\n    全绿，未触数据与分析。",
            fontsize=10, va="top", linespacing=1.6)


def seal(fig):
    ev = {}
    for tag in ("P1", "P2", "P3"):
        d = f"results/{tag.lower()}"
        fr = json.loads(open(f"{d}/{tag}_FREEZE.json").read())
        ev[tag] = {
            "prop": dt.datetime.fromtimestamp(
                os.stat(f"{d}/{tag}_PROPOSITION.md").st_mtime, dt.timezone.utc),
            "freeze": dt.datetime.fromisoformat(fr["freeze_created_utc"]),
            "stats": dt.datetime.fromtimestamp(
                os.stat(f"{d}/{tag}_summary.json").st_mtime, dt.timezone.utc),
        }
    t0 = min(e["prop"] for e in ev.values()).replace(minute=35, second=0, microsecond=0)
    t1 = max(e["stats"] for e in ev.values()).replace(minute=30, second=0, microsecond=0)
    span = (t1 - t0).total_seconds() / 60
    W = 500

    def px(t):
        return 4 + (t - t0).total_seconds() / 60 / span * (W - 8)

    ax = V.mm_axes(fig, 70, 126, W, 96, True)
    ax.set_xlim(0, W)
    ax.set_ylim(0, 100)
    ax.set_yticks([])
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(V.PAL["ink"])
    ax.spines["bottom"].set_linewidth(0.6)

    ticks, t = [], t0.replace(minute=(t0.minute // 10) * 10, second=0, microsecond=0)
    while t <= t1:
        if t >= t0:
            ticks.append(t)
        t += dt.timedelta(minutes=10)
    ax.set_xticks([px(t) for t in ticks])
    ax.set_xticklabels([t.strftime("%H:%M") for t in ticks], fontsize=10)
    ax.tick_params(axis="x", length=3, pad=2, direction="out", color=V.PAL["ink"])

    for i, tag in enumerate(("P1", "P2", "P3")):
        e = ev[tag]
        y = 82 - i * 22
        ax.add_patch(plt_rect(px(e["freeze"]), y - 4.5,
                              px(e["stats"]) - px(e["freeze"]), 9,
                              facecolor=V.PAL["seal"], edgecolor="none", zorder=1))
        ax.plot([px(e["prop"]), px(e["stats"])], [y, y], color=V.PAL["indigo"], lw=0.8,
                zorder=2)
        ax.plot(px(e["prop"]), y, "o", ms=4.4, mfc=V.PAL["paper"], mec=V.PAL["ink"],
                mew=0.9, zorder=4)
        ax.plot(px(e["freeze"]), y, "D", ms=3.6, mfc=V.PAL["indigo"], mec=V.PAL["ink"],
                mew=0.7, zorder=4)
        ax.plot(px(e["stats"]), y, "s", ms=3.8, mfc=V.PAL["dai"], mec=V.PAL["ink"],
                mew=0.7, zorder=4)
        ax.text(px(e["prop"]), y + 5, "命题落盘", fontsize=10, ha="center", va="bottom")
        ax.text(px(e["freeze"]), y - 6, "冻结哈希", fontsize=10, ha="center", va="top")
        ax.text(px(e["stats"]), y + 5, "统计落盘", fontsize=10, ha="center", va="bottom")

    ax.set_yticks([82 - i * 22 for i in range(3)])
    ax.set_yticklabels(["P1", "P2", "P3"], fontsize=11)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.add_patch(plt_rect(6, 92, 7, 4, facecolor=V.PAL["seal"], edgecolor="none"))
    ax.text(15, 94, "浅色带 = 封存期：哈希之后、统计之前，不得读该轮检验集残差",
            fontsize=10, ha="left", va="center", color=V.PAL["dai"])


def main():
    V.apply_style()
    df, splits = load()
    df = order(df)
    fig = V.mm_figure(PAGE_W, PAGE_H)
    fig.text(
        M / PAGE_W,
        1 - 13 / PAGE_H,
        "三轮之间，191 个测定各自去了哪里",
        fontsize=18,
        va="top",
    )
    fig.text(
        M / PAGE_W,
        1 - 23 / PAGE_H,
        "每个测定在三轮里只当过一次检验集；同一 UniProt 的全部测定整组同侧。"
        "柱高为该轮可用突变量（对数轴）。",
        fontsize=10,
        va="top",
        color=V.PAL["dai"],
    )
    strip(fig, df)
    brackets(fig, df)
    legend(fig, df)
    seal(fig)
    groups(fig, df)
    fig.text(M / PAGE_W, 8 / PAGE_H,
             "数据：P*_FORMATION / P*_TEST.csv（去向）、master_tables（突变量与双点数）、"
             "各 FREEZE.json（封条时刻）。未改动任何结果文件。", fontsize=10, va="bottom",
             color=V.PAL["dai"])
    V.check_text(fig)
    V.save(fig, "ALLOC_A2")


if __name__ == "__main__":
    main()
