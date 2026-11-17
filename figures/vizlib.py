"""Shared palette, typography and layout checks for the verification-loop plate.

Read-only with respect to results/: this module only styles figures.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

PAL = {
    "indigo": "#2F5D6B",
    "ochre": "#9A6A4F",
    "carmine": "#9C4A54",
    "gamboge": "#C0A04A",
    "dai": "#46505C",
    "grey": "#8C8779",
    "ink": "#2B2B2B",
    "paper": "#FBF8F1",
    "seal": "#CBD6D6",
    "seen": "#ADA697",
}

MM = 1 / 25.4


def apply_style():
    plt.rcParams.update(
        {
            "font.family": ["STHeiti"],
            "font.size": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.facecolor": PAL["paper"],
            "savefig.facecolor": PAL["paper"],
            "axes.facecolor": PAL["paper"],
            "axes.edgecolor": PAL["ink"],
            "axes.linewidth": 0.6,
            "text.color": PAL["ink"],
            "axes.labelcolor": PAL["ink"],
            "xtick.color": PAL["ink"],
            "ytick.color": PAL["ink"],
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 10,
        }
    )


def mm_axes(fig, x, y, w, h, decorated=False):
    """Axes whose data coordinates are millimetres of the printed figure."""
    ax = fig.add_axes([x / fig.w_mm, y / fig.h_mm, w / fig.w_mm, h / fig.h_mm])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.set_facecolor("none")
    if not decorated:
        ax.set_axis_off()
    return ax


def mm_figure(w_mm, h_mm):
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM), dpi=72)
    fig.w_mm, fig.h_mm = w_mm, h_mm
    return fig


def text_boxes(fig, dpi=300):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    out = []
    for t in fig.findobj(Text):
        if not t.get_visible():
            continue
        s = t.get_text()
        if not s or not s.strip():
            continue
        bb = t.get_window_extent(r)
        if bb.width <= 0 or bb.height <= 0:
            continue
        out.append((s, bb))
    return out


def check_text(fig, pad=0.5, verbose=True):
    """Report text-text overlap and out-of-figure text, in points."""
    boxes = text_boxes(fig)
    W = fig.w_mm * 72 / 25.4
    H = fig.h_mm * 72 / 25.4
    issues = []
    for s, bb in boxes:
        if bb.x0 < -pad or bb.y0 < -pad or bb.x1 > W + pad or bb.y1 > H + pad:
            issues.append(("out", s, (bb.x0, bb.y0, bb.x1, bb.y1)))
    for i in range(len(boxes)):
        si, bi = boxes[i]
        for j in range(i + 1, len(boxes)):
            sj, bj = boxes[j]
            ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if ox > pad and oy > pad:
                issues.append(("overlap", (si, sj), (round(ox, 2), round(oy, 2))))
    if verbose:
        print(f"texts: {len(boxes)}  issues: {len(issues)}")
        for it in issues[:40]:
            print("  ", it)
    return issues


def save(fig, stem, png_dpi=150):
    fig.savefig(f"figures/{stem}.pdf", format="pdf")
    fig.savefig(f"figures/{stem}.png", dpi=png_dpi)
    print(f"saved figures/{stem}.pdf/.png")
