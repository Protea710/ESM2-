"""Build the four-figure plate (A2). Panel A: round-dependency graph."""

import json
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge

sys.path.insert(0, "figures")
import vizlib as V

AA = list("ACDEFGHIKLMNPQRSTVWY")
ROUNDS = [
    ("P1", "theta", "REFUTED", 0.10, "单点 · margin 四分位"),
    ("P2", "theta_norm", "INDETERMINATE", 0.05, "单点 · LLR/R_i 归一"),
    ("P3", "theta_3", "REFUTED", 0.10, "双点 · 非加性"),
]

PAGE_W, PAGE_H = 594, 420
MARGIN = 18
COL_L, COL_R = 259, 285
ROW_H, ROW_GAP = 174, 16
X_L = MARGIN
X_R = MARGIN + COL_L + 14
Y_BOT = MARGIN
Y_TOP = MARGIN + ROW_H + ROW_GAP

NODES = [
    # key, x, count (None = not tested), verdict, label inside/outside
    ("委托书", 26, 2411093, None),
    ("A 段", 72, 2411093, None),
    ("P1", 118, 223686, "REFUTED"),
    ("P2", 164, 145471, "INDETERMINATE"),
    ("P3", 210, 130340, "REFUTED"),
    ("P4", 250, None, "PROPOSED"),
]

VERDICT_FILL = {
    "REFUTED": V.PAL["carmine"],
    "INDETERMINATE": V.PAL["gamboge"],
    None: V.PAL["paper"],
    "PROPOSED": V.PAL["paper"],
}


def node_radius(n):
    return 0.5 * 20.0 * math.sqrt(n / 1e6)


def draw_node(ax, x, y, name, n, verdict, r=None):
    if r is None:
        r = node_radius(n) if n else 4.5
    if verdict == "INDETERMINATE":
        ax.add_patch(
            Circle(
                (x, y),
                r,
                facecolor=V.PAL["paper"],
                edgecolor=V.PAL["ink"],
                lw=1.0,
                zorder=3,
            )
        )
        ax.add_patch(
            Wedge(
                (x, y),
                r,
                90,
                270,
                facecolor=V.PAL["gamboge"],
                edgecolor="none",
                zorder=3.1,
            )
        )
        ax.add_patch(
            Circle(
                (x, y), r, facecolor="none", edgecolor=V.PAL["ink"], lw=1.0, zorder=3.2
            )
        )
    elif verdict == "PROPOSED":
        ax.add_patch(
            Circle(
                (x, y),
                r,
                facecolor=V.PAL["paper"],
                edgecolor=V.PAL["indigo"],
                lw=1.0,
                ls=(0, (3, 2)),
                zorder=3,
            )
        )
    else:
        ax.add_patch(
            Circle(
                (x, y),
                r,
                facecolor=VERDICT_FILL[verdict],
                edgecolor=V.PAL["ink"],
                lw=1.0,
                zorder=3,
            )
        )
        if verdict == "REFUTED":
            d = r * 0.52
            ax.plot(
                [x - d, x + d],
                [y - d, y + d],
                color=V.PAL["paper"],
                lw=2.0,
                solid_capstyle="round",
                zorder=3.3,
            )
            ax.plot(
                [x - d, x + d],
                [y + d, y - d],
                color=V.PAL["paper"],
                lw=2.0,
                solid_capstyle="round",
                zorder=3.3,
            )
    return r


def panel_a(fig, ax):
    W, H = COL_L, ROW_H
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    y0 = 120

    ax.text(0, 171, "A", fontsize=11, va="top", color=V.PAL["indigo"], weight="bold")
    ax.text(
        9, 171, "三轮：推翻 / 无法判断 / 推翻——每轮的命题由上轮的数字逼出来",
        fontsize=13, va="top",
    )

    for name, x, n, verdict in NODES:
        draw_node(ax, x, y0, name, n, verdict)

    for i in range(len(NODES) - 1):
        x0, x1 = NODES[i][1], NODES[i + 1][1]
        r0 = node_radius(NODES[i][2]) if NODES[i][2] else 4.5
        r1 = node_radius(NODES[i + 1][2]) if NODES[i + 1][2] else 4.5
        ax.add_patch(
            FancyArrowPatch(
                (x0 + r0 + 1.2, y0), (x1 - r1 - 1.2, y0),
                arrowstyle="-|>", mutation_scale=9, lw=0.9,
                color=V.PAL["indigo"], shrinkA=0, shrinkB=0, zorder=2,
            )
        )

    labels = [
        "清单与筛选规则冻结\nR2+R3 → 191 测定",
        "ρ 0.45042\n官方 0.45046\nmax 0.000506",
        "θ=−0.1390\n[−0.1767,−0.0839]\n48/61 负 → 归一",
        "θ_norm=+0.0699\n[+0.0012,+0.1289]\n未判定 → 转双点",
        "θ3=+0.1366\n[+0.0022,+0.2904]\n18 正 / 8 负",
    ]
    for i, lab in enumerate(labels):
        xm = (NODES[i][1] + NODES[i + 1][1]) / 2
        ax.plot([xm, xm], [y0 + 1.0, 141.0], color=V.PAL["grey"], lw=0.5, zorder=1)
        ax.text(xm, 142.0, lab, fontsize=10, ha="center", va="bottom", linespacing=1.35)

    ax.text(26, y0, "委托书", fontsize=10, ha="center", va="center")
    ax.text(72, y0, "A 段", fontsize=10, ha="center", va="center")
    ax.text(26, y0 - 17.0, "2,411,093 突变\n191 测定", fontsize=10,
            ha="center", va="top", linespacing=1.35)
    ax.text(72, y0 - 17.0, "2,411,093 突变\n全量重算", fontsize=10,
            ha="center", va="top", linespacing=1.35)
    ax.text(118, y0 - 7.4, "P1 · 223,686\n单点 · 推翻", fontsize=10, ha="center",
            va="top", linespacing=1.35)
    ax.text(164, y0 - 5.8, "P2 · 145,471\n单点 · 无法判断", fontsize=10, ha="center",
            va="top", linespacing=1.35)
    t_p3 = ax.text(210, y0 - 5.6, "P3 · 130,340\n双点 · 推翻", fontsize=10,
                   ha="center", va="top", linespacing=1.35)
    ax.text(258, y0 - 7.2, "P4 提议\n未检验", fontsize=10, ha="right", va="top",
            linespacing=1.35)

    ax.text(4, 76, "判定：颜色 + 形状", fontsize=10, va="top")
    for verdict, label, yy in [
        ("REFUTED", "推翻", 62),
        ("INDETERMINATE", "无法判断", 50),
        (None, "成立（本轮无）", 38),
    ]:
        draw_node(ax, 9, yy, "", None, verdict, r=3.5)
        ax.text(15, yy, label, fontsize=10, va="center")

    ax.text(72, 76, "突变量 n → 节点面积", fontsize=10, va="top")
    ax.add_patch(Circle((79, 62), node_radius(1e5), facecolor="none",
                        edgecolor=V.PAL["grey"], lw=0.8))
    ax.text(85.5, 62, "10 万", fontsize=10, va="center")
    ax.add_patch(Circle((84, 42), node_radius(1e6), facecolor="none",
                        edgecolor=V.PAL["grey"], lw=0.8))
    ax.text(98, 42, "100 万", fontsize=10, va="center")

    fig.canvas.draw()
    bb = t_p3.get_window_extent(fig.canvas.get_renderer()).transformed(
        ax.transData.inverted()
    )
    ax.plot([205, 205], [72, bb.y0 - 0.4], color=V.PAL["grey"], lw=0.6,
            ls=(0, (3, 2)), zorder=1)
    ax.text(150, 74, "中途执行修正（如实记入）", fontsize=10, va="top")
    ax.text(
        150, 67,
        "① CSV 误读（P3 首跑）：读错数据源，\n   改读 A 段产物后重跑，口径不变\n"
        "② 门禁误删（P3 重提交）：清单误删后恢复，\n   未触数据与分析",
        fontsize=10, va="top", linespacing=1.4,
    )



def load_rounds():
    out = {}
    for tag, key, verdict, thr, note in ROUNDS:
        d = json.loads(open(f"results/{tag.lower()}/{tag}_summary.json").read())
        out[tag] = {
            "theta": d[key], "lo": d["ci_low"], "hi": d["ci_high"],
            "verdict": d["verdict"], "thr": thr,
            "n_test": d["n_test_assays"], "n_ok": d["n_eligible"], "note": note,
        }
    return out


def load_heat(dms_id, singles):
    d = singles[singles.DMS_id == dms_id]
    piv = d.pivot_table(index="position", columns="mut_aa", values="r")
    piv = piv.reindex(columns=AA)
    return piv.to_numpy(float), len(piv)


def ramp():
    from matplotlib.colors import LinearSegmentedColormap

    base = LinearSegmentedColormap.from_list(
        "zresid", [V.PAL["indigo"], "#8F9794", V.PAL["paper"], "#C2AC97",
                   V.PAL["ochre"]], N=512
    )
    cols = [base(t) for t in np.linspace(0, 1, 15)]
    cmap = ListedColormap(cols)
    cmap.set_bad(V.PAL["paper"])
    return cmap


def verdict_marker(ax, x, y, verdict, r_pt=5.2):
    """Verdict glyph in data units, sized in points so it stays round at any dpi."""
    from matplotlib.patches import Ellipse, Polygon

    pos = ax.get_position()
    w_mm = pos.width * ax.figure.w_mm
    h_mm = pos.height * ax.figure.h_mm
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()
    sx = (xhi - xlo) / w_mm * 0.35278 * r_pt
    sy = (yhi - ylo) / h_mm * 0.35278 * r_pt
    ax.add_patch(Ellipse((x, y), 2 * sx, 2 * sy, facecolor=V.PAL["paper"],
                         edgecolor="none", zorder=6))
    if verdict == "INDETERMINATE":
        t = np.linspace(np.pi / 2, 3 * np.pi / 2, 40)
        verts = [(x, y)] + [(x + sx * np.cos(a), y + sy * np.sin(a)) for a in t]
        ax.add_patch(Polygon(verts, closed=True, facecolor=V.PAL["gamboge"],
                             edgecolor="none", zorder=6.1))
    elif verdict == "REFUTED":
        ax.add_patch(Ellipse((x, y), 2 * sx, 2 * sy, facecolor=V.PAL["carmine"],
                             edgecolor="none", zorder=6.1))
        d = 0.5
        ax.plot([x - d * sx, x + d * sx], [y - d * sy, y + d * sy],
                color=V.PAL["paper"], lw=1.5, solid_capstyle="round", zorder=6.3)
        ax.plot([x - d * sx, x + d * sx], [y + d * sy, y - d * sy],
                color=V.PAL["paper"], lw=1.5, solid_capstyle="round", zorder=6.3)
    ax.add_patch(Ellipse((x, y), 2 * sx, 2 * sy, facecolor="none",
                         edgecolor=V.PAL["ink"], lw=0.9, zorder=6.2))


def panel_b(fig):
    W, H = COL_R, ROW_H
    lay = V.mm_axes(fig, X_R, Y_TOP, W, H)
    lay.set_xlim(0, W)
    lay.set_ylim(0, H)

    lay.text(0, 171, "B", fontsize=11, va="top", color=V.PAL["indigo"], weight="bold")
    lay.text(9, 171, "证据图：残差结构（左）与三轮效应量及区间（右）",
             fontsize=13, va="top")

    singles = pd.read_csv("results/master_tables/singles.csv.gz",
                          usecols=["DMS_id", "position", "mut_aa", "r"])
    cmap = ramp()
    norm = BoundaryNorm(np.linspace(-3, 3, 16), cmap.N)

    mats = [
        ("SUMO1_HUMAN_Weile_2017",
         V.mm_axes(fig, X_R + 20, Y_TOP + 50, 84, 108, True), 84),
        ("TNKS2_HUMAN_Tsuboyama_2023_5JRT",
         V.mm_axes(fig, X_R + 108, Y_TOP + 50, 28, 108, True), 28),
        ("MLAC_ECOLI_MacRae_2023",
         V.mm_axes(fig, X_R + 140, Y_TOP + 50, 28, 108, True), 28),
    ]
    clipped = total = 0
    for dms, ax, w in mats:
        mat, npos = load_heat(dms, singles)
        total += int(np.isfinite(mat).sum())
        clipped += int((np.abs(mat[np.isfinite(mat)]) > 3).sum())
        ax.pcolormesh(np.arange(21), np.arange(npos + 1),
                      np.ma.masked_invalid(mat), cmap=cmap, norm=norm,
                      shading="flat", antialiased=False, rasterized=False)
        ax.set_xlim(0, 20)
        ax.set_ylim(npos, 0)
        ax.set_xticks(np.arange(20) + 0.5)
        ax.tick_params(axis="x", length=0, pad=0)
        if w == 84:
            ax.set_xticklabels(AA, fontsize=10)
            ticks = list(range(1, npos + 1, 20)) + [npos]
            ax.set_yticks(ticks)
            ax.set_yticklabels([str(t) for t in ticks], fontsize=10)
            ax.tick_params(axis="y", length=2, pad=2)
            ax.set_ylabel("位点序号", fontsize=11, labelpad=4)
            ax.set_xlabel("替换氨基酸", fontsize=11, labelpad=2)
        else:
            ax.set_xticklabels([])
            ax.set_yticks([])
        for side in ax.spines.values():
            side.set_linewidth(0.6)
            side.set_color(V.PAL["grey"])

    lay.text(20, 42, "SUMO1_HUMAN_Weile_2017 · 96 位点 · 1,700 单点",
             fontsize=10, va="top")
    lay.text(108, 160, "TNKS2", fontsize=10, va="bottom")
    lay.text(140, 160, "MLAC", fontsize=10, va="bottom")
    lay.text(20, 35, "小图：TNKS2_HUMAN_Tsuboyama_2023_5JRT 59 位点（1,118）；"
                     "MLAC_ECOLI_MacRae_2023 211 位点（4,007）", fontsize=10, va="top")

    cax = V.mm_axes(fig, X_R + 172, Y_TOP + 60, 4, 86, True)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      ticks=[-3, -2, -1, 0, 1, 2, 3])
    cb.ax.yaxis.set_ticks_position("right")
    cb.ax.tick_params(labelsize=10, length=2.5, pad=3, direction="out")
    cb.outline.set_linewidth(0.6)
    cb.outline.set_edgecolor(V.PAL["grey"])
    lay.text(172, 154.5, "残差 r", fontsize=10, va="bottom")
    lay.text(172, 150, "z 单位", fontsize=10, va="bottom")

    rnd = load_rounds()
    fax = V.mm_axes(fig, X_R + 214, Y_TOP + 54, 68, 94, True)
    fax.set_xlim(-0.25, 0.33)
    fax.set_ylim(0, 3)
    for side in ["top", "right", "left"]:
        fax.spines[side].set_visible(False)
    fax.spines["bottom"].set_linewidth(0.6)
    fax.spines["bottom"].set_color(V.PAL["ink"])
    fax.set_yticks([])
    fax.set_xticks([-0.2, 0, 0.2])
    fax.set_xticklabels(["−0.2", "0", "+0.2"], fontsize=10)
    fax.tick_params(axis="x", length=3, pad=2)
    fax.set_xlabel("效应量 θ（z 单位）", fontsize=11, labelpad=2)

    fax.axvline(0, color=V.PAL["grey"], lw=0.6, ls=(0, (4, 2)))
    for v in (-0.05, 0.05, 0.10):
        fax.axvline(v, color=V.PAL["grey"], lw=0.6, ls=(0, (1.5, 1.8)), ymax=0.97)

    for i, tag in enumerate(["P1", "P2", "P3"]):
        d = rnd[tag]
        y = 2.42 - i * 0.94
        fax.plot([d["lo"], d["hi"]], [y, y], color=V.PAL["ink"], lw=1.4,
                 solid_capstyle="butt", zorder=4)
        for e in (d["lo"], d["hi"]):
            fax.plot([e, e], [y - 0.05, y + 0.05], color=V.PAL["ink"], lw=1.2,
                     zorder=4)
        verdict_marker(fax, d["theta"], y, d["verdict"], r_pt=5.2)
        fax.text(d["theta"], y + 0.13, f"{d['theta']:+.4f}".replace("-", "−"),
                 fontsize=10, ha="center", va="bottom")
        fax.text((d["lo"] + d["hi"]) / 2, y - 0.18,
                 f"[{d['lo']:+.4f}, {d['hi']:+.4f}]".replace("-", "−"),
                 fontsize=10, ha="center", va="top")
        fax.text(-0.03, y, f"{tag} · {d['n_ok']}", fontsize=10, ha="right",
                 va="center", transform=fax.get_yaxis_transform())

    fax.text(-0.05, 2.90, "±0.05", fontsize=10, ha="right", va="top")
    fax.text(0.105, 2.90, "0.10", fontsize=10, ha="left", va="top")

    lay.text(0, 24,
             "残差量级 ÷ 效应量量级（同为 z 单位；残差量级 = 各测定（低组均值+高组均值）/2 的测定中位数）",
             fontsize=10, va="top")
    for xx, txt in [
        (0, "P1  n=61   0.79 ÷ 0.139 = 5.7"),
        (84, "P2  n=42   0.93 ÷ 0.070 = 13.3"),
        (168, "P3  n=26   1.02 ÷ 0.137 = 7.5"),
    ]:
        lay.text(xx, 17, txt, fontsize=10, va="top")
    print(f"heatmap cells beyond ±3 z: {clipped}/{total} ({clipped/total:.3%})")


def panel_c(fig):
    W, H = COL_L, ROW_H
    lay = V.mm_axes(fig, X_L, Y_BOT, W, H)
    lay.set_xlim(0, W)
    lay.set_ylim(0, H)
    lay.text(0, 171, "C", fontsize=11, va="top", color=V.PAL["indigo"], weight="bold")
    lay.text(9, 171, "校准：三测定的实验值与模型分（各图自有标尺）",
             fontsize=13, va="top")

    sing = pd.read_csv("results/master_tables/singles.csv.gz",
                       usecols=["DMS_id", "DMS_score", "model_score"])
    rho = pd.read_csv("results/assay_summary.csv").set_index("DMS_id")["rho_single"]

    specs = [
        ("SUMO1_HUMAN_Weile_2017", "SUMO1 · 1,700 单点", V.PAL["indigo"], 13, 87),
        ("TNKS2_HUMAN_Tsuboyama_2023_5JRT", "TNKS2 · 1,118 单点", V.PAL["ochre"], 97, 171),
        ("MLAC_ECOLI_MacRae_2023", "MLAC · 4,007 单点", V.PAL["dai"], 181, 255),
    ]
    kern = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], float)
    for i, (dms, name, col, x0, x1) in enumerate(specs):
        d = sing[sing.DMS_id == dms]
        x, y = d.DMS_score.to_numpy(float), d.model_score.to_numpy(float)
        ax = V.mm_axes(fig, X_L + x0, Y_BOT + 48, x1 - x0, 104, True)
        H2, xe, ye = np.histogram2d(x, y, bins=36)
        Hs = np.zeros_like(H2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                w = kern[dx + 1, dy + 1]
                Hs += w * np.roll(np.roll(H2, dx, 0), dy, 1)
        xc = 0.5 * (xe[:-1] + xe[1:])
        yc = 0.5 * (ye[:-1] + ye[1:])
        Hs[Hs < 1.5] = 0.0
        lv = np.array([0.08, 0.22, 0.45, 0.72]) * Hs.max()
        ax.contour(xc, yc, Hs, levels=lv, colors=[col] * 4,
                   linewidths=[0.5, 0.7, 0.9, 1.1])
        ax.set_xlim(xe[0], xe[-1])
        ax.set_ylim(ye[0], ye[-1])
        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)
        for side in ["bottom", "left"]:
            ax.spines[side].set_linewidth(0.6)
            ax.spines[side].set_color(V.PAL["grey"])
        ax.tick_params(labelsize=10, length=2.5, pad=2)
        from matplotlib.ticker import FuncFormatter

        fmt = FuncFormatter(lambda v, _: f"{v:g}".replace("-", "−"))
        ystep = 5 if (ye[-1] - ye[0]) > 15 else 2.5
        yt = np.arange(np.ceil(ye[0] / ystep) * ystep, ye[-1] + 1e-9, ystep)
        xt = np.arange(np.ceil(xe[0] * 2) / 2, xe[-1] + 1e-9, 0.5)
        ax.set_yticks(yt)
        ax.set_xticks(xt)
        ax.yaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_formatter(fmt)
        ax.set_title(f"{name.split(' · ')[0]} · ρ = {rho[dms]:.5f}", fontsize=10,
                     pad=4, color=col)
        if i == 0:
            ax.set_ylabel("模型分 LLR", fontsize=11, labelpad=3)

    lay.text(134, 40, "实验值 DMS_score", fontsize=11, va="top", ha="center")
    lay.text(4, 30,
             "A 段 191/191 测定（含多点）：mean ρ 0.45042，官方 0.45046，"
             "max|Δρ| 0.000506，无一超 0.02", fontsize=10, va="top")
    lay.text(4, 22,
             "图内 ρ 为该测定单点 Spearman（assay_summary.csv 列 rho_single）；"
             "等高线为 36×36 分箱 + 3×3 平滑后 4 条等值线（平滑计数 <1.5 的格不画）",
             fontsize=10, va="top")


def event_times():
    import datetime as dt
    import os

    out = {}
    for tag in ("P1", "P2", "P3"):
        d = f"results/{tag.lower()}"
        fr = json.loads(open(f"{d}/{tag}_FREEZE.json").read())
        out[tag] = {
            "prop": dt.datetime.fromtimestamp(
                os.stat(f"{d}/{tag}_PROPOSITION.md").st_mtime, dt.timezone.utc),
            "freeze": dt.datetime.fromisoformat(fr["freeze_created_utc"]),
            "stats": dt.datetime.fromtimestamp(
                os.stat(f"{d}/{tag}_summary.json").st_mtime, dt.timezone.utc),
        }
    return out


def panel_d(fig):
    W, H = COL_R, ROW_H
    lay = V.mm_axes(fig, X_R, Y_BOT, W, H)
    lay.set_xlim(0, W)
    lay.set_ylim(0, H)
    lay.text(0, 171, "D", fontsize=11, va="top", color=V.PAL["indigo"], weight="bold")
    lay.text(9, 171, "冻结早于统计：九个时刻，三轮同一顺序",
             fontsize=13, va="top")

    ev = event_times()
    t0 = min(e["prop"] for e in ev.values()).replace(minute=35, second=0, microsecond=0)
    t1 = max(e["stats"] for e in ev.values()).replace(minute=30, second=0, microsecond=0)
    span = (t1 - t0).total_seconds() / 60
    x0, x1 = 12, W - 6
    base = 60

    def px(t):
        return x0 + (t - t0).total_seconds() / 60 / span * (x1 - x0)

    ax = V.mm_axes(fig, X_R, Y_BOT + base, W, H - base, True)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H - base)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["bottom"].set_color(V.PAL["ink"])
    ax.set_yticks([])

    import datetime as dt

    ticks, t = [], t0.replace(minute=(t0.minute // 10) * 10,
                              second=0, microsecond=0)
    while t <= t1:
        if t >= t0:
            ticks.append(t)
        t += dt.timedelta(minutes=10)
    ax.set_xticks([px(t) for t in ticks])
    ax.set_xticklabels([t.strftime("%H:%M") for t in ticks], fontsize=10)
    ax.tick_params(axis="x", length=3, pad=2, color=V.PAL["ink"])
    for t in ticks:
        ax.plot([px(t), px(t)], [0, H - base - 4], color=V.PAL["grey"], lw=0.3,
                zorder=0)

    for tag, y in [("P1", 90), ("P2", 64), ("P3", 38)]:
        e = ev[tag]
        ax.text(4, y, tag, fontsize=10, va="center", ha="left",
                color=V.PAL["indigo"])
        ax.plot([px(e["prop"]), px(e["stats"])], [y, y], color=V.PAL["indigo"],
                lw=0.7, zorder=2)
        ax.plot(px(e["prop"]), y, "o", ms=4.2, mfc=V.PAL["paper"], mec=V.PAL["ink"],
                mew=0.9, zorder=4)
        ax.plot(px(e["freeze"]), y, "D", ms=3.4, mfc=V.PAL["indigo"],
                mec=V.PAL["ink"], mew=0.7, zorder=4)
        ax.plot(px(e["stats"]), y, "s", ms=3.6, mfc=V.PAL["dai"],
                mec=V.PAL["ink"], mew=0.7, zorder=4)
        ax.text(px(e["prop"]), y + 4.0, "命题落盘", fontsize=10, ha="center",
                va="bottom")
        ax.text(px(e["freeze"]), y - 4.4, "冻结哈希", fontsize=10, ha="center",
                va="top")
        ax.text(px(e["stats"]), y + 4.0, "统计落盘", fontsize=10, ha="center",
                va="bottom")
        xm = (px(e["freeze"]) + px(e["stats"])) / 2
        ax.annotate("", xy=(px(e["stats"]) - 1.6, y - 2.0),
                    xytext=(px(e["freeze"]) + 1.6, y - 2.0),
                    arrowprops=dict(arrowstyle="-|>", color=V.PAL["ochre"],
                                    lw=0.7, mutation_scale=7))
        ax.text(xm + 8, y - 2.3, "哈希在前", fontsize=10, ha="center", va="top",
                color=V.PAL["ochre"])

    lay.text(4, 50, "横轴为 UTC 时刻（刻度 10 分钟）；不标任何耗时。", fontsize=10,
             va="top")
    for i, tag in enumerate(["P1", "P2", "P3"]):
        e = ev[tag]
        lay.text(4, 40 - i * 7.5,
                 f"{tag}   命题 {e['prop']:%H:%M:%S}   哈希 {e['freeze']:%H:%M:%S}"
                 f"   统计 {e['stats']:%H:%M:%S}", fontsize=10, va="top")


def main():
    V.apply_style()
    fig = V.mm_figure(PAGE_W, PAGE_H)
    fig.text(
        MARGIN / PAGE_W,
        1 - 14 / PAGE_H,
        "模型知道自己错在哪吗：三轮预登记检验的执行链",
        fontsize=18,
        va="top",
    )
    fig.text(
        MARGIN / PAGE_W,
        1 - 25 / PAGE_H,
        "全部数字只读自 results/ 与 docs/，逐项出处见复核说明",
        fontsize=10,
        va="top",
        color=V.PAL["dai"],
    )
    ax = V.mm_axes(fig, X_L, Y_TOP, COL_L, ROW_H)
    panel_a(fig, ax)
    panel_b(fig)
    panel_c(fig)
    panel_d(fig)
    V.check_text(fig)
    V.save(fig, "FIG_PLATE_A2")


if __name__ == "__main__":
    main()
