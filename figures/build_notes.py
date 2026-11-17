"""A4 companion pages: 边界表 (p1) and 复核说明 (p2)."""

import sys

import matplotlib.pyplot as plt

sys.path.insert(0, "figures")
import vizlib as V

W, H = 210, 297
M = 16

DO = [
    "191 测定按 masked-marginals 口径重算零样本分",
    "逐测定 ρ 与官方 CSV 逐位对照，191/191 全对",
    "三轮命题先落盘再检验，哈希与时刻留档",
    "效应量 + 10,000 次测定级 bootstrap 区间",
    "三分类判定写死，按冻结规则报结论",
    "负结果如实入档（P1、P3），不换指标",
    "逐突变总表 2,411,093 行，含残差与模型内部特征",
    "单命令重跑，逐字节比对 A 段与三轮统计",
    "下一轮命题 + 值得测的理由（P4）",
]

DONT = [
    "不微调、不训练、不改动模型权重",
    "不排名，不与其他模型或测定比高低",
    "不写生物学结论，不给机制解释",
    "不放宽阈值、不换指标、不做事后挑选",
    "不读检验集残差之后再改命题",
    "不用 MSA、物种、结构、保守性作解释变量",
    "不省略或合并区间，不隐藏零结果",
    "不改动任何结果文件与冻结件",
    "不把结论外推到这 191 个测定之外",
]

ROWS = [
    (
        "A",
        "2,411,093 / 191 测定",
        "两节点同值同径",
        "MANIFEST.json total_rows_all_tables；assay_list_proposed.csv 191 行",
    ),
    (
        "A",
        "223,686（P1 单点）",
        "进入 P1 统计的行数",
        "singles.csv.gz 筛选 p1_set=test ∩ 61 合格测定 ∩ n_single_at_pos≥8",
    ),
    (
        "A",
        "145,471（P2 单点）",
        "同上口径",
        "singles.csv.gz 同筛选，p2_set=test ∩ 42 合格测定",
    ),
    (
        "A",
        "130,340（P3 双点）",
        "Σn_dbl",
        "P3_per_assay.csv 列 n_dbl 求和 = P3_FREEZE census.test_doubles",
    ),
    ("A", "θ=−0.1390 [−0.1767,−0.0839]", "P1 判定依据", "P1_summary.json"),
    ("A", "θ_norm=+0.0699 [+0.0012,+0.1289]", "P2 判定依据", "P2_summary.json"),
    ("A", "θ3=+0.1366 [+0.0022,+0.2904]", "P3 判定依据", "P3_summary.json"),
    ("A", "48/61 负向（P1）", "见差异 D1", "P1_per_assay.csv 按 delta_a<0 计数"),
    ("A", "18 正 / 8 负（P3）", "方向计数", "P3_per_assay.csv；P3_RESULTS.md"),
    (
        "A",
        "0.45042 / 0.45046 / 0.000506",
        "A 段均值与最大绝对差",
        "results/compare.csv 191 行；FINAL_REPORT §1",
    ),
    (
        "A",
        "3b8ae8f3 / ef0e814e / 618840c2",
        "命题哈希前 8 位",
        "FINAL_REPORT §1 = 各 P*_FREEZE.json artifacts[0].sha256",
    ),
    (
        "B",
        "1,700 / 1,118 / 4,007 单点",
        "热图三测定行数",
        "singles.csv.gz 该 DMS_id 行数",
    ),
    ("B", "96 / 59 / 211 位点", "热图纵轴长度", "同上，position 唯一值数"),
    ("B", "P1 61 / P2 42 / P3 26", "合格测定数", "各 P*_summary.json n_eligible"),
    (
        "B",
        "0.79 / 0.93 / 1.02 残差量级",
        "（低组均值+高组均值）/2 的测定中位数",
        "P*_per_assay.csv（P3 取 mean_absr_low/high）",
    ),
    ("B", "5.7 / 13.3 / 7.5", "残差量级 ÷ |θ|", "上两行之商，未做任何新统计"),
    (
        "B",
        "±3 z 色标",
        "42/6,825 格（0.62%）超端点按端点色",
        "singles.csv.gz 三测定 |r|>3 计数",
    ),
    (
        "C",
        "ρ = 0.50931 / 0.59304 / 0.38208",
        "该测定单点 Spearman",
        "assay_summary.csv 列 rho_single",
    ),
    ("C", "191/191，无一超 0.02", "|Δρ| 检查", "compare.csv 列 delta；FINAL_REPORT §1"),
    (
        "D",
        "九个时刻",
        "哈希=freeze_created_utc，其余=文件 mtime（UTC）",
        "P*_FREEZE.json；P*_PROPOSITION.md 与 P*_summary.json 的 mtime",
    ),
]

DIFF = [
    (
        "D1",
        "P1 方向计数",
        "P1_per_assay.csv 实测 48 负 / 13 正；此前 FINAL_REPORT §1 与 P1_RESULTS.md 误写 50/61、11 正，"
        "已按 ERRATA E1 更正为 48/61。图与表值一致。",
    ),
    (
        "D2",
        "P1 统计时刻",
        "P1_summary.json mtime = 17:51:38Z；P1_RESULTS.md 正文写「17:52Z 之后」。"
        "相差 ≤22 s，不影响「哈希早于统计」。图用 mtime。",
    ),
    (
        "D3",
        "P3 哈希指代",
        "P3_RESULTS.md 的 81fa6b81… 是 FREEZE.json 文件自身哈希；图为命题文件哈希 618840c2…，"
        "与 FINAL_REPORT §1 一致。二者不是同一对象，非冲突。",
    ),
]


def wrap(s, width_mm, size_pt):
    import re

    per = size_pt * 0.3528

    def cw(ch):
        return per if ord(ch) > 127 else per * 0.52

    toks = re.findall(r"[A-Za-z0-9_.,/%+\-\[\]()=:]+|.", s)
    lines, cur, w = [], "", 0.0
    for tok in toks:
        tw = sum(cw(c) for c in tok)
        if w + tw <= width_mm:
            cur += tok
            w += tw
            continue
        if cur.strip():
            lines.append(cur)
        cur, w = "", 0.0
        if tw <= width_mm:
            cur, w = tok, tw
        else:
            for c in tok:
                if w + cw(c) > width_mm and cur:
                    lines.append(cur)
                    cur, w = "", 0.0
                cur += c
                w += cw(c)
    if cur:
        lines.append(cur)
    return lines


def block(fig, x_mm, y_mm, s, width_mm, size_pt, colour=None, lead=1.45):
    lines = wrap(s, width_mm, size_pt)
    yy = y_mm
    for ln in lines:
        fig.text(x_mm / W, 1 - yy / H, ln, fontsize=size_pt, va="top",
                 color=colour or V.PAL["ink"])
        yy += size_pt * 0.3528 * lead
    return yy


def page1(fig):
    fig.text(M / W, 1 - M / H, "边界：做 / 不做", fontsize=18, va="top")
    block(fig, M, M + 10,
          "取自 docs/BAR.md、docs/FINAL_REPORT.md 与三轮 RESULT 文件的陈述边界；"
          "只描述这次交付。", 178, 10, V.PAL["dai"])

    for i, (head, items) in enumerate([("做", DO), ("不做", DONT)]):
        x = M + i * 92
        fig.text(x / W, 1 - (M + 28) / H, head, fontsize=14, va="top",
                 color=V.PAL["indigo"] if i == 0 else V.PAL["carmine"])
        y = M + 40
        for n, it in enumerate(items, 1):
            fig.text(x / W, 1 - y / H, f"{n:02d}", fontsize=9.5, va="top",
                     color=V.PAL["grey"])
            y = block(fig, x + 8, y, it, 76, 11.5) + 5.0

    yb = 172
    fig.add_artist(plt.Line2D([M / W, (M + 178) / W], [1 - yb / H, 1 - yb / H],
                              color=V.PAL["grey"], lw=0.5))
    block(fig, M, yb + 6,
          "「硬规则：命题被推翻算完成；只有「没跑完」「没说清为什么」算未完成；"
          "结论不超出这批测定。」——docs/BAR.md", 178, 10.5, V.PAL["dai"])
    block(fig, M, yb + 22,
          "这张边界页与四图一页、复核说明一页同属一次交付：图上一切数字在复核说明里能定位到"
          "文件与列；本页只划定这次没有做、也不打算做的事。", 178, 10.5)


def page2(fig):
    fig.text(M / W, 1 - M / H, "复核说明：图上每个数字的出处", fontsize=18, va="top")
    block(fig, M, M + 9,
          "全部为只读读取。未重算任何统计量；除行数计数与两处由既有列直接相除的比值外，"
          "图上数字都是文件里的原值。图列为四图的面板字母 A–D。", 178, 9, V.PAL["dai"])

    xs = [M, M + 11, M + 66, M + 106]
    widths = [10, 54, 38, 70]
    heads = ["图", "标注值", "含义", "出处（文件 / 列 / 行）"]
    y = M + 22
    for x, h in zip(xs, heads):
        fig.text(x / W, 1 - y / H, h, fontsize=9, va="top", color=V.PAL["dai"])
    y += 1.2
    fig.add_artist(plt.Line2D([M / W, (M + 178) / W], [1 - y / H, 1 - y / H],
                              color=V.PAL["grey"], lw=0.5))
    y += 4.4
    for tag, val, mean, src in ROWS:
        cells = [tag, val, mean, src]
        ends = []
        for x, w, c in zip(xs, widths, cells):
            ends.append(block(fig, x, y, c, w, 9, lead=1.4))
        y = max(ends) + 1.6

    y += 3
    fig.text(M / W, 1 - y / H, "三处需要说明的地方（如实列出，未改动任何文件）",
             fontsize=11, va="top", color=V.PAL["carmine"])
    y += 7
    for tag, what, txt in DIFF:
        fig.text(M / W, 1 - y / H, f"{tag}  {what}", fontsize=9, va="top")
        y = block(fig, M + 26, y, txt, 150, 8.5) + 2.2

    y += 2
    block(fig, M, y,
          "读法：「行数」由 results/master_tables/*.csv.gz 逐行计数得到，筛选条件即各轮冻结"
          "脚本的口径；未对任何数值做重算、平滑或合并。", 176, 8.5, V.PAL["dai"])


def main():
    from matplotlib.backends.backend_pdf import PdfPages

    V.apply_style()
    pages = []
    for builder in (page1, page2):
        fig = V.mm_figure(W, H)
        builder(fig)
        V.check_text(fig)
        pages.append(fig)
    with PdfPages("figures/NOTES_A4.pdf") as pdf:
        for fig in pages:
            pdf.savefig(fig, facecolor=V.PAL["paper"])
    for i, fig in enumerate(pages, 1):
        fig.savefig(f"figures/NOTES_A4_p{i}.png", dpi=150)
    print("saved figures/NOTES_A4.pdf (2 pages) + p1/p2 PNG")


if __name__ == "__main__":
    main()
