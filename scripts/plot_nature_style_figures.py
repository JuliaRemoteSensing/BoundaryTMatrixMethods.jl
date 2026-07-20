from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper_validation" / "data"
FIG = ROOT / "paper_validation" / "figures_nature"
REPORT = ROOT / "paper_validation" / "reports" / "nature_figure_manifest.md"
FIG.mkdir(parents=True, exist_ok=True)
REPORT.parent.mkdir(parents=True, exist_ok=True)


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7.6,
        "axes.titlesize": 8.4,
        "axes.labelsize": 7.7,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "figure.titlesize": 9.8,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "lines.linewidth": 1.25,
        "lines.markersize": 4.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "savefig.dpi": 450,
        "figure.dpi": 170,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


PALETTE = {
    "ink": "#1C1C1C",
    "muted": "#666666",
    "grid": "#DADADA",
    "pale_grid": "#ECECEC",
    "paper": "#FFFFFF",
    "wash": "#F7F7F7",
    "blue": "#3B6EA8",
    "teal": "#149C88",
    "orange": "#D96B27",
    "purple": "#8E6BBE",
    "rose": "#C76D9D",
    "grey": "#777777",
    "light_teal": "#DFF1EE",
    "light_orange": "#F7E8DC",
    "light_blue": "#E9F0F8",
    "light_purple": "#EFEAF7",
}

METHOD = {
    "EBCM": {"color": PALETTE["ink"], "linestyle": "-", "linewidth": 1.9},
    "SVM": {"color": PALETTE["blue"], "linestyle": (0, (3.2, 1.7)), "linewidth": 1.35},
    "PMM": {"color": PALETTE["teal"], "linestyle": "-", "linewidth": 1.35},
    "IITM coarse": {"color": PALETTE["orange"], "linestyle": (0, (3.0, 1.4, 1.0, 1.4)), "linewidth": 1.05},
    "IITM fine": {"color": PALETTE["grey"], "linestyle": ":", "linewidth": 1.4},
    "Auto selector": {"color": PALETTE["rose"], "linestyle": "-", "linewidth": 1.2},
}

ERR_CMAP = LinearSegmentedColormap.from_list(
    "nature_error",
    ["#F7F7F7", "#D7E7F0", "#8CBBD2", "#3F84A8", "#153B5A"],
)

CASE_COLORS = {
    "sphere": "#3B6EA8",
    "prolate_mild": "#149C88",
    "oblate_mild": "#7B9E3A",
    "chebyshev_n2": "#C54E4B",
    "prolate_strong": "#D96B27",
    "oblate_strong": "#A6611A",
    "chebyshev_n4_pos": "#8E6BBE",
    "chebyshev_n4_neg": "#C76D9D",
}

CASE_LABEL = {
    "sphere": "Sphere",
    "prolate_mild": "Mild prolate",
    "oblate_mild": "Mild oblate",
    "chebyshev_n2": "Cheb. n=2",
    "prolate_strong": "Strong prolate",
    "oblate_strong": "Strong oblate",
    "chebyshev_n4_pos": "Cheb. n=4+",
    "chebyshev_n4_neg": "Cheb. n=4-",
}

CASE_PARAMS = {
    "sphere": ("spheroid", 1.0, 1.0, 0, 0.0),
    "prolate_mild": ("spheroid", 0.8, 1.2, 0, 0.0),
    "oblate_mild": ("spheroid", 1.2, 0.8, 0, 0.0),
    "chebyshev_n2": ("chebyshev", 1.0, 1.0, 2, 0.08),
    "prolate_strong": ("spheroid", 0.65, 1.35, 0, 0.0),
    "oblate_strong": ("spheroid", 1.35, 0.65, 0, 0.0),
    "chebyshev_n4_pos": ("chebyshev", 1.0, 1.0, 4, 0.08),
    "chebyshev_n4_neg": ("chebyshev", 1.0, 1.0, 4, -0.08),
}

ELEMENTS = ["F11", "F12", "F22", "F33", "F34", "F44"]


def load_data():
    return (
        pd.read_csv(DATA / "validation_metrics.csv"),
        pd.read_csv(DATA / "validation_mueller.csv"),
        pd.read_csv(DATA / "validation_gate.csv"),
        pd.read_csv(DATA / "validation_selector.csv"),
        pd.read_csv(DATA / "validation_convergence.csv"),
    )


def save(fig, name):
    png = FIG / f"{name}.png"
    pdf = FIG / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [png, pdf]


def panel_label(ax, label, x=-0.12, y=1.05):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.8,
        fontweight="bold",
        color=PALETTE["ink"],
    )


def ordered_cases(frame, group=None):
    sub = frame if group is None else frame[frame["group"] == group]
    return list(dict.fromkeys(sub["case"].tolist()))


def shape_profile(case, n=240):
    theta = np.linspace(0, 2 * np.pi, n)
    kind, a, c, order, eps = CASE_PARAMS[case]
    if kind == "spheroid":
        r = np.ones_like(theta)
    else:
        r = 1.0 + eps * np.cos(order * theta)
    x = a * r * np.sin(theta)
    z = c * r * np.cos(theta)
    return x, z


def draw_shape(ax, case, accepted=True, lw=1.45):
    x, z = shape_profile(case)
    color = PALETTE["teal"] if accepted else PALETTE["orange"]
    fill = PALETTE["light_teal"] if accepted else PALETTE["light_orange"]
    ax.fill(x, z, facecolor=fill, edgecolor=color, lw=lw)
    ax.plot([0, 0], [min(z), max(z)], color="#CFCFCF", lw=0.45, zorder=0)
    ax.plot([min(x), max(x)], [0, 0], color="#CFCFCF", lw=0.45, zorder=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def metric(metrics, case, method, col):
    return float(metrics.loc[(metrics["case"] == case) & (metrics["method"] == method), col].iloc[0])


def method_legend(methods):
    handles = []
    for method in methods:
        style = METHOD[method]
        handles.append(Line2D([0], [0], label=method, **style))
    return handles


def light_grid(ax, axis="both"):
    ax.grid(True, axis=axis, color=PALETTE["pale_grid"], lw=0.55)
    ax.set_axisbelow(True)


def fig01_graphical_abstract(metrics, gate):
    fig = plt.figure(figsize=(7.2, 3.2))
    gs = fig.add_gridspec(2, 7, height_ratios=[1.0, 0.82], width_ratios=[1, 1, 0.18, 1.15, 0.18, 1, 1])

    ax_title = fig.add_subplot(gs[:, :])
    ax_title.set_axis_off()
    fig.suptitle(
        "Boundary-matching T matrices with a no-reference PMM gate",
        y=0.985,
        fontsize=10.4,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.925,
        "SVM/PMM -> TransitionMatrices.jl T matrix -> Mueller observables -> gate or IITM fallback",
        ha="center",
        va="center",
        fontsize=7.4,
        color=PALETTE["muted"],
    )

    shape_cases = ["sphere", "prolate_mild", "chebyshev_n2", "prolate_strong"]
    for i, case in enumerate(shape_cases[:2]):
        ax = fig.add_subplot(gs[0, i])
        accepted = bool(gate.loc[gate["case"] == case, "accepted"].iloc[0])
        draw_shape(ax, case, accepted)
        ax.set_title(CASE_LABEL[case], pad=1.5)
    for i, case in enumerate(shape_cases[2:]):
        ax = fig.add_subplot(gs[1, i])
        accepted = bool(gate.loc[gate["case"] == case, "accepted"].iloc[0])
        draw_shape(ax, case, accepted)
        ax.set_title(CASE_LABEL[case], pad=1.5)

    ax_mid = fig.add_subplot(gs[:, 3])
    ax_mid.set_axis_off()
    boxes = [
        (0.05, 0.64, 0.90, 0.22, "surface\nboundary system", PALETTE["light_blue"]),
        (0.05, 0.37, 0.40, 0.18, "SVM", "#EEF4FB"),
        (0.55, 0.37, 0.40, 0.18, "PMM", "#EAF5F1"),
        (0.05, 0.09, 0.90, 0.20, "T matrix\nobject", PALETTE["light_purple"]),
    ]
    for x, y, w, h, txt, fc in boxes:
        rect = patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.015,rounding_size=0.02",
            facecolor=fc,
            edgecolor=PALETTE["ink"],
            lw=0.8,
        )
        ax_mid.add_patch(rect)
        ax_mid.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.0)
    ax_mid.annotate("", xy=(0.25, 0.55), xytext=(0.45, 0.64), arrowprops=dict(arrowstyle="-|>", lw=0.8))
    ax_mid.annotate("", xy=(0.75, 0.55), xytext=(0.55, 0.64), arrowprops=dict(arrowstyle="-|>", lw=0.8))
    ax_mid.annotate("", xy=(0.43, 0.29), xytext=(0.25, 0.37), arrowprops=dict(arrowstyle="-|>", lw=0.8))
    ax_mid.annotate("", xy=(0.57, 0.29), xytext=(0.75, 0.37), arrowprops=dict(arrowstyle="-|>", lw=0.8))

    ax_curve = fig.add_subplot(gs[0, 5:])
    methods = ["EBCM", "PMM", "SVM"]
    mueller = pd.read_csv(DATA / "validation_mueller.csv")
    for method in methods:
        sub = mueller[(mueller["case"] == "prolate_mild") & (mueller["method"] == method)].sort_values("angle_deg")
        ax_curve.plot(sub["angle_deg"], sub["F11"], **METHOD[method])
    ax_curve.set_title("validated Mueller curve")
    ax_curve.set_xlabel("angle")
    ax_curve.set_ylabel("F11")
    ax_curve.set_xticks([0, 90, 180])
    light_grid(ax_curve)
    ax_curve.legend(handles=method_legend(methods), frameon=False, loc="upper right")

    ax_gate = fig.add_subplot(gs[1, 5:])
    cases = ordered_cases(gate)
    y = np.arange(len(cases))
    vals = gate.set_index("case").loc[cases, "convergence_error"].astype(float).values
    colors = [PALETTE["teal"] if bool(gate.set_index("case").loc[c, "accepted"]) else PALETTE["orange"] for c in cases]
    ax_gate.barh(y, vals, color=colors, height=0.56)
    ax_gate.axvline(5e-3, color=PALETTE["ink"], lw=0.8, ls="--")
    ax_gate.set_xscale("log")
    ax_gate.set_yticks([0, 3, 7])
    ax_gate.set_yticklabels([CASE_LABEL[cases[i]] for i in [0, 3, 7]])
    ax_gate.set_xlabel("PMM stability")
    ax_gate.set_title("gate prevents blind PMM")
    light_grid(ax_gate, "x")

    fig.subplots_adjust(left=0.03, right=0.985, bottom=0.10, top=0.80, wspace=0.55, hspace=0.55)
    return save(fig, "nature_fig01_graphical_abstract")


def fig02_case_library(metrics, gate, selector):
    cases = ordered_cases(metrics)
    fig = plt.figure(figsize=(7.2, 5.7))
    gs = fig.add_gridspec(4, 8, height_ratios=[0.95, 0.22, 1.25, 1.25], hspace=0.30, wspace=0.16)

    for i, case in enumerate(cases):
        ax = fig.add_subplot(gs[0, i])
        accepted = bool(gate.loc[gate["case"] == case, "accepted"].iloc[0])
        draw_shape(ax, case, accepted)
        ax.set_title(CASE_LABEL[case].replace(" ", "\n"), fontsize=6.6, pad=1.0)
        ax.text(
            0.5,
            -0.16,
            "pass" if accepted else "reject",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.4,
            color=PALETTE["teal"] if accepted else PALETTE["orange"],
        )

    ax = fig.add_subplot(gs[2, :])
    x = np.arange(len(cases))
    pmm = [metric(metrics, c, "PMM", "mueller_relerr_vs_ebcm") for c in cases]
    auto = [metric(metrics, c, "Auto selector", "mueller_relerr_vs_ebcm") for c in cases]
    svm = [metric(metrics, c, "SVM", "mueller_relerr_vs_ebcm") for c in cases]
    ax.plot(x, svm, color=PALETTE["blue"], marker="o", lw=1.0, label="SVM")
    ax.plot(x, pmm, color=PALETTE["teal"], marker="s", lw=1.0, label="PMM")
    ax.plot(x, auto, color=PALETTE["rose"], marker="D", lw=1.0, label="Auto")
    ax.set_yscale("log")
    ax.set_ylabel("Mueller error")
    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.set_title("Observable accuracy across the benchmark library", pad=5)
    light_grid(ax, "y")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[3, :])
    conv = gate.set_index("case").loc[cases, "convergence_error"].astype(float).values
    residual = gate.set_index("case").loc[cases, "boundary_residual"].astype(float).values
    ax.scatter(x - 0.08, conv, color=PALETTE["purple"], s=25, label="PMM stability")
    ax.scatter(x + 0.08, residual, color=PALETTE["orange"], s=25, label="boundary residual")
    ax.axhline(5e-3, color=PALETTE["purple"], lw=0.8, ls="--")
    ax.axhline(0.5, color=PALETTE["orange"], lw=0.8, ls=":")
    ax.set_yscale("log")
    ax.set_ylabel("gate metric")
    ax.set_xticks(x)
    ax.set_xticklabels([CASE_LABEL[c] for c in cases], rotation=25, ha="right")
    ax.set_title("No-reference diagnostics used by the selector", pad=5)
    light_grid(ax, "y")
    ax.legend(frameon=False, ncol=2, loc="lower left")
    panel_label(ax, "B")

    fig.suptitle("Benchmark particle library and diagnostic coverage", y=0.975, fontweight="bold")
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.11, top=0.87)
    return save(fig, "nature_fig02_case_library")


def fig03_accuracy_forest(metrics):
    cases = ordered_cases(metrics)
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    offsets = np.linspace(-0.25, 0.25, len(methods))
    markers = {"SVM": "o", "PMM": "s", "IITM coarse": "^", "IITM fine": "v", "Auto selector": "D"}

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    y = np.arange(len(cases))
    for i, method in enumerate(methods):
        vals = [metric(metrics, c, method, "mueller_relerr_vs_ebcm") for c in cases]
        ax.scatter(vals, y + offsets[i], color=METHOD[method]["color"], marker=markers[method], s=22, label=method, zorder=3)
        for v, yy in zip(vals, y + offsets[i]):
            ax.plot([1e-15, v], [yy, yy], color=METHOD[method]["color"], lw=0.35, alpha=0.35, zorder=1)

    for idx, case in enumerate(cases):
        if "strong" in case or "n4" in case:
            ax.axhspan(idx - 0.5, idx + 0.5, color=PALETTE["light_orange"], zorder=0, alpha=0.55)
    ax.set_xscale("log")
    ax.set_xlim(3e-16, 5e-2)
    ax.set_yticks(y)
    ax.set_yticklabels([CASE_LABEL[c] for c in cases])
    ax.invert_yaxis()
    ax.set_xlabel("relative Mueller-matrix error vs EBCM")
    ax.set_title("Accuracy forest plot: SVM/PMM T matrices and IITM fallback")
    light_grid(ax, "x")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    panel_label(ax, "A", -0.07, 1.02)
    return save(fig, "nature_fig03_accuracy_forest")


def fig04_mueller_accepted(mueller):
    cases = ordered_cases(mueller, "accepted")
    methods = ["EBCM", "SVM", "PMM", "IITM fine"]
    fig, axes = plt.subplots(len(cases), len(ELEMENTS), figsize=(7.4, 5.8), sharex=True)
    for r, case in enumerate(cases):
        for c, element in enumerate(ELEMENTS):
            ax = axes[r, c]
            for method in methods:
                sub = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values("angle_deg")
                z = 5 if method == "EBCM" else 3
                ax.plot(sub["angle_deg"], sub[element], zorder=z, **METHOD[method])
            if r == 0:
                ax.set_title(element, pad=3.0)
            if r < len(cases) - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xticks([0, 90, 180])
            ax.tick_params(axis="y", labelsize=5.7)
            ax.grid(True, color=PALETTE["pale_grid"], lw=0.45)
            ax.set_xlim(0, 180)
    fig.legend(handles=method_legend(methods), frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.54, 0.94))
    fig.text(0.51, 0.035, "scattering angle (deg)", ha="center")
    row_y = np.linspace(0.745, 0.205, len(cases))
    for y, case in zip(row_y, cases):
        fig.text(0.030, y, CASE_LABEL[case], ha="center", va="center", rotation=90, fontsize=7.3)
    fig.suptitle("Six-element Mueller validation for gate-accepted PMM/SVM T matrices", y=0.985, fontweight="bold")
    fig.subplots_adjust(left=0.115, right=0.995, bottom=0.09, top=0.855, wspace=0.28, hspace=0.28)
    return save(fig, "nature_fig04_mueller_accepted_atlas")


def element_errors(mueller, cases, method):
    arr = np.zeros((len(cases), len(ELEMENTS)))
    for i, case in enumerate(cases):
        ref = mueller[(mueller["case"] == case) & (mueller["method"] == "EBCM")].sort_values("angle_deg")
        cur = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values("angle_deg")
        for j, element in enumerate(ELEMENTS):
            num = np.linalg.norm(cur[element].to_numpy() - ref[element].to_numpy())
            den = max(np.linalg.norm(ref[element].to_numpy()), np.finfo(float).eps)
            arr[i, j] = num / den
    return arr


def heatmap_with_values(ax, data, row_labels, col_labels, title, cmap="mako"):
    log_data = np.log10(np.clip(data, 1e-16, None))
    im = ax.imshow(log_data, aspect="auto", cmap=ERR_CMAP, vmin=-6, vmax=-1)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=0)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_title(title)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            color = "white" if np.log10(max(val, 1e-16)) > -2.6 else PALETTE["ink"]
            ax.text(j, i, f"{val:.1e}", ha="center", va="center", fontsize=5.8, color=color)
    return im


def fig05_stress_deviation(mueller, metrics):
    stress = ordered_cases(mueller, "stress")
    labels = [CASE_LABEL[c] for c in stress]
    pmm_err = element_errors(mueller, stress, "PMM")
    auto_err = element_errors(mueller, stress, "Auto selector")

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], width_ratios=[1, 1], hspace=0.45, wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    im = heatmap_with_values(ax1, pmm_err, labels, ELEMENTS, "PMM deviation by Mueller element", "viridis")
    panel_label(ax1, "A", -0.20, 1.08)
    ax2 = fig.add_subplot(gs[0, 1])
    heatmap_with_values(ax2, auto_err, labels, ELEMENTS, "Auto selector deviation", "viridis")
    ax2.set_yticklabels([])
    panel_label(ax2, "B", -0.12, 1.08)
    cax = fig.add_axes([0.92, 0.57, 0.012, 0.25])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log10(error)")

    ax3 = fig.add_subplot(gs[1, :])
    x = np.arange(len(stress))
    pmm = [metric(metrics, c, "PMM", "mueller_relerr_vs_ebcm") for c in stress]
    auto = [metric(metrics, c, "Auto selector", "mueller_relerr_vs_ebcm") for c in stress]
    iitm = [metric(metrics, c, "IITM coarse", "mueller_relerr_vs_ebcm") for c in stress]
    ax3.plot(x, pmm, color=PALETTE["teal"], marker="s", label="PMM")
    ax3.plot(x, auto, color=PALETTE["rose"], marker="D", label="Auto selector")
    ax3.plot(x, iitm, color=PALETTE["orange"], marker="^", label="IITM fallback")
    for i in range(len(stress)):
        ax3.plot([i, i], [min(pmm[i], auto[i]), max(pmm[i], auto[i])], color="#BBBBBB", lw=0.55, zorder=0)
    ax3.set_yscale("log")
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=18, ha="right")
    ax3.set_ylabel("whole Mueller error")
    ax3.set_title("Rejected PMM candidates are not used blindly")
    light_grid(ax3, "y")
    ax3.legend(frameon=False, ncol=3, loc="upper right")
    panel_label(ax3, "C", -0.08, 1.07)

    fig.suptitle("Stress-case failure modes and IITM fallback", y=0.985, fontweight="bold")
    fig.subplots_adjust(left=0.12, right=0.90, bottom=0.12, top=0.88)
    return save(fig, "nature_fig05_stress_deviation_atlas")


def fig06_gate_diagnostics(metrics, gate, selector):
    cases = ordered_cases(gate)
    base = gate.set_index("case")
    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], width_ratios=[1.1, 1])

    ax = fig.add_subplot(gs[0, 0])
    label_cases = {"sphere", "chebyshev_n2", "prolate_strong", "oblate_strong"}
    offsets = {
        "sphere": (3, -2),
        "chebyshev_n2": (4, -8),
        "prolate_strong": (4, 9),
        "oblate_strong": (4, -10),
    }
    for case in cases:
        accepted = bool(base.loc[case, "accepted"])
        color = PALETTE["teal"] if accepted else PALETTE["orange"]
        x = float(base.loc[case, "convergence_error"])
        y = metric(metrics, case, "PMM", "mueller_relerr_vs_ebcm")
        s = 35 + 120 * min(float(base.loc[case, "boundary_residual"]) / 0.5, 1.0)
        ax.scatter(x, y, s=s, color=color, edgecolor=PALETTE["ink"], lw=0.45, zorder=3)
        if case in label_cases:
            ax.annotate(
                CASE_LABEL[case],
                xy=(x, y),
                xytext=offsets.get(case, (4, 4)),
                textcoords="offset points",
                fontsize=5.9,
                color=PALETTE["ink"],
                arrowprops=dict(arrowstyle="-", lw=0.35, color="#8A8A8A", shrinkA=2, shrinkB=2),
            )
    ax.axvline(5e-3, color=PALETTE["ink"], ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2e-16, 7e-2)
    ax.set_ylim(2e-16, 5e-2)
    ax.set_xlabel("refined-grid PMM change")
    ax.set_ylabel("PMM Mueller error")
    ax.set_title("PMM stability predicts usable regime")
    light_grid(ax)
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    crit = np.zeros((len(cases), 4))
    crit[:, 0] = base.loc[cases, "convergence_error"].astype(float).values <= 5e-3
    crit[:, 1] = base.loc[cases, "boundary_residual"].astype(float).values <= 0.5
    crit[:, 2] = base.loc[cases, "projected_boundary_residual"].astype(float).values <= 0.5
    crit[:, 3] = base.loc[cases, "energy_violation"].astype(float).values <= 1e-10
    cmap = mpl.colors.ListedColormap([PALETTE["orange"], PALETTE["teal"]])
    ax.imshow(crit, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(["stable", "resid.", "proj.", "passive"], rotation=28, ha="right")
    ax.set_yticks(np.arange(len(cases)))
    ax.set_yticklabels([CASE_LABEL[c] for c in cases])
    ax.set_title("Gate criteria")
    for i in range(len(cases)):
        for j in range(4):
            ax.text(j, i, "yes" if crit[i, j] else "no", ha="center", va="center", color="white", fontsize=6.3)
    panel_label(ax, "B", -0.18)

    ax = fig.add_subplot(gs[1, :])
    y = np.arange(len(cases))
    selected = selector.set_index("case").loc[cases, "selected_method"].tolist()
    colors = [PALETTE["teal"] if s == "pmm" else PALETTE["orange"] for s in selected]
    ax.barh(y, np.ones(len(cases)), color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([CASE_LABEL[c] for c in cases])
    ax.set_xticks([])
    ax.invert_yaxis()
    for i, s in enumerate(selected):
        ax.text(0.5, i, "PMM accepted" if s == "pmm" else "IITM fallback", ha="center", va="center", color="white", fontsize=7.4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Automatic method selection")
    panel_label(ax, "C", -0.07)

    fig.suptitle("Residual diagnostic gate: why PMM is accepted or rejected", y=0.985, fontweight="bold")
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.08, top=0.88, wspace=0.42, hspace=0.55)
    return save(fig, "nature_fig06_gate_diagnostics")


def fig07_convergence(convergence):
    cases = ordered_cases(convergence)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.7), sharex=True)
    axes = axes.ravel()
    for idx, (ax, case) in enumerate(zip(axes, cases)):
        for method in ["SVM", "PMM"]:
            sub = convergence[(convergence["case"] == case) & (convergence["method"] == method)].sort_values("nmax")
            ax.plot(
                sub["nmax"],
                sub["mueller_relerr_vs_ebcm"],
                marker="o",
                color=METHOD[method]["color"],
                linestyle=METHOD[method]["linestyle"],
                label=method,
            )
        ax.set_yscale("log")
        ax.set_title(CASE_LABEL[case])
        ax.set_xlabel("nmax")
        ax.set_ylabel("Mueller error")
        light_grid(ax)
        if idx == 0:
            ax.legend(frameon=False)
            panel_label(ax, "A")
    fig.suptitle("Convergence of SVM/PMM boundary-matching T matrices", y=0.99, fontweight="bold")
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.10, top=0.88, wspace=0.30, hspace=0.45)
    return save(fig, "nature_fig07_convergence")


def fig08_observable_decomposition(metrics):
    accepted = ordered_cases(metrics, "accepted")
    obs = [
        ("mueller_relerr_vs_ebcm", "Mueller"),
        ("F11_relerr_vs_ebcm", "F11"),
        ("tmatrix_relerr_vs_ebcm", "T matrix"),
        ("Csca_relerr_vs_ebcm", "Csca"),
        ("Cext_relerr_vs_ebcm", "Cext"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), gridspec_kw={"width_ratios": [1.35, 1.0]})

    ax = axes[0]
    x = np.arange(len(obs))
    for i, case in enumerate(accepted):
        vals = [max(metric(metrics, case, "PMM", col), 1e-16) for col, _ in obs]
        ax.plot(x, vals, marker="o", lw=1.0, label=CASE_LABEL[case], alpha=0.95, color=CASE_COLORS[case])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([name for _, name in obs], rotation=25, ha="right")
    ax.set_ylabel("relative error vs EBCM")
    ax.set_title("PMM error decomposed by observable")
    light_grid(ax, "y")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    panel_label(ax, "A")

    ax = axes[1]
    methods = ["SVM", "PMM", "IITM fine"]
    vals = []
    for method in methods:
        vals.append([metric(metrics, case, method, "Cext_relerr_vs_ebcm") for case in accepted])
    vals = np.array(vals)
    im = ax.imshow(np.log10(np.clip(vals, 1e-16, None)), aspect="auto", cmap=ERR_CMAP, vmin=-15, vmax=-2)
    ax.set_xticks(np.arange(len(accepted)))
    ax.set_xticklabels([CASE_LABEL[c] for c in accepted], rotation=28, ha="right")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    ax.set_title("Extinction cross-section error")
    cb = fig.colorbar(im, ax=ax, shrink=0.72)
    cb.set_label("log10(error)")
    panel_label(ax, "B")

    fig.suptitle("T-matrix validation should include observables and cross sections", y=0.99, fontweight="bold")
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.22, top=0.82, wspace=0.42)
    return save(fig, "nature_fig08_observable_decomposition")


def fig09_api_reproducibility():
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    columns = [
        (0.05, "Solvers", ["svm_tmatrix", "pmm_tmatrix", "residual_corrected_iitm"]),
        (0.29, "Diagnostics", ["boundary_residual", "projected_residual", "diagnose_pmm_applicability"]),
        (0.56, "Selector", ["auto_tmatrix", "PMM if safe", "IITM fallback"]),
        (0.79, "Outputs", ["AxisymmetricTransitionMatrix", "Mueller curves", "CSV + figures"]),
    ]
    for x, title, items in columns:
        rect = patches.FancyBboxPatch(
            (x, 0.22),
            0.17,
            0.55,
            boxstyle="round,pad=0.018,rounding_size=0.025",
            facecolor="#F8F8F8",
            edgecolor=PALETTE["ink"],
            lw=0.75,
        )
        ax.add_patch(rect)
        ax.text(x + 0.085, 0.70, title, ha="center", va="center", fontsize=8.7, fontweight="bold")
        for k, item in enumerate(items):
            ax.text(x + 0.085, 0.58 - 0.12 * k, item, ha="center", va="center", fontsize=7.0, color=PALETTE["muted"])
    for x0, x1 in [(0.22, 0.29), (0.46, 0.56), (0.73, 0.79)]:
        ax.annotate("", xy=(x1, 0.50), xytext=(x0, 0.50), arrowprops=dict(arrowstyle="-|>", lw=0.85, color=PALETTE["ink"]))

    code = (
        "T = auto_tmatrix(shape, 2pi, nmax; pmm_ngauss=120)\n"
        "F = scattering_matrix(T.T, 2pi, angles)\n"
        "gate = diagnose_pmm_applicability(shape, 2pi, nmax)"
    )
    rect = patches.FancyBboxPatch(
        (0.18, 0.03),
        0.64,
        0.12,
        boxstyle="round,pad=0.014,rounding_size=0.018",
        facecolor=PALETTE["light_blue"],
        edgecolor="#B8C6D8",
        lw=0.65,
    )
    ax.add_patch(rect)
    ax.text(0.50, 0.09, code, ha="center", va="center", family="monospace", fontsize=6.7)
    ax.text(0.50, 0.91, "Package-level reproducibility and API design", ha="center", va="center", fontsize=10.3, fontweight="bold")
    ax.text(
        0.50,
        0.84,
        "All validation data and figures are regenerated from scripts, not manual plotting.",
        ha="center",
        va="center",
        fontsize=7.6,
        color=PALETTE["muted"],
    )
    return save(fig, "nature_fig09_api_reproducibility")


def fig10_manuscript_roadmap():
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    rows = [
        ("1. Motivation", "T matrices need reusable software interfaces", "Fig. 1"),
        ("2. Methods", "SVM/PMM boundary matching -> TransitionMatrices object", "Fig. 1, 9"),
        ("3. Diagnostics", "PMM is gated by stability, residual and passivity", "Fig. 6"),
        ("4. Validation", "sphere, spheroids and Chebyshev particles", "Fig. 2"),
        ("5. Accuracy", "six Mueller elements and cross sections", "Fig. 3, 4, 8"),
        ("6. Failure modes", "stress cases trigger IITM fallback", "Fig. 5, 6"),
        ("7. Reproducibility", "CSV data, scripts, tests and examples", "Fig. 9"),
    ]
    y0 = 0.84
    dy = 0.105
    for i, (section, claim, figs) in enumerate(rows):
        y = y0 - i * dy
        color = PALETTE["light_teal"] if i in [2, 4] else PALETTE["light_blue"] if i in [1, 6] else "#F8F8F8"
        rect = patches.FancyBboxPatch(
            (0.05, y - 0.045),
            0.90,
            0.072,
            boxstyle="round,pad=0.01,rounding_size=0.018",
            facecolor=color,
            edgecolor="#D0D0D0",
            lw=0.55,
        )
        ax.add_patch(rect)
        ax.text(0.08, y - 0.009, section, ha="left", va="center", fontsize=7.6, fontweight="bold")
        ax.text(0.31, y - 0.009, claim, ha="left", va="center", fontsize=7.2, color=PALETTE["muted"])
        ax.text(0.91, y - 0.009, figs, ha="right", va="center", fontsize=7.0, color=PALETTE["ink"])
    ax.text(0.5, 0.95, "Manuscript story map", ha="center", va="center", fontsize=10.5, fontweight="bold")
    ax.text(
        0.5,
        0.905,
        "The CPC paper should be written as a software and reliability workflow, not as a claim of a new SVM/PMM algorithm.",
        ha="center",
        va="center",
        fontsize=7.2,
        color=PALETTE["muted"],
    )
    return save(fig, "nature_fig10_manuscript_roadmap")


def supp_case_mueller(mueller):
    paths = []
    methods = ["EBCM", "SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    for case in ordered_cases(mueller):
        group = mueller.loc[mueller["case"] == case, "group"].iloc[0]
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.2), sharex=True)
        axes = axes.ravel()
        for ax, element in zip(axes, ELEMENTS):
            for method in methods:
                sub = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values("angle_deg")
                style = METHOD[method].copy()
                if method in ["IITM coarse", "IITM fine", "Auto selector"]:
                    style["alpha"] = 0.72
                ax.plot(sub["angle_deg"], sub[element], **style)
            ax.set_title(element)
            ax.set_xlim(0, 180)
            ax.set_xticks([0, 90, 180])
            ax.grid(True, color=PALETTE["pale_grid"], lw=0.45)
        for ax in axes[3:]:
            ax.set_xlabel("angle (deg)")
        axes[0].set_ylabel("value")
        axes[3].set_ylabel("value")
        fig.legend(handles=method_legend(methods), frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.52, 0.94))
        status = "gate accepted" if group == "accepted" else "stress / PMM rejected"
        fig.suptitle(f"{CASE_LABEL[case]}: six Mueller elements ({status})", y=0.995, fontweight="bold")
        fig.subplots_adjust(left=0.08, right=0.99, bottom=0.10, top=0.82, wspace=0.30, hspace=0.38)
        paths.extend(save(fig, f"nature_supp_six_elements_{case}"))
    return paths


def write_manifest(paths):
    main = [
        ("nature_fig01_graphical_abstract", "Graphical abstract and package workflow."),
        ("nature_fig02_case_library", "Benchmark shapes, gate statuses, and diagnostic coverage."),
        ("nature_fig03_accuracy_forest", "Forest plot of Mueller-matrix errors across methods and cases."),
        ("nature_fig04_mueller_accepted_atlas", "Six-element Mueller validation for gate-accepted cases."),
        ("nature_fig05_stress_deviation_atlas", "Stress-case element deviations and IITM fallback."),
        ("nature_fig06_gate_diagnostics", "No-reference gate criteria and selector decisions."),
        ("nature_fig07_convergence", "Convergence of SVM/PMM against EBCM."),
        ("nature_fig08_observable_decomposition", "Observable-level error decomposition."),
        ("nature_fig09_api_reproducibility", "Package API and reproducibility workflow."),
        ("nature_fig10_manuscript_roadmap", "Section-by-section manuscript story map."),
    ]
    lines = [
        "# Nature-Style Figure Manifest",
        "",
        "Figures are generated from `paper_validation/data/*.csv` by",
        "`scripts/plot_nature_style_figures.py`.",
        "",
        "## Main Figure Set",
        "",
    ]
    for stem, desc in main:
        lines.append(f"- `paper_validation/figures_nature/{stem}.png` - {desc}")
        lines.append(f"- `paper_validation/figures_nature/{stem}.pdf`")
    lines.extend(["", "## Supplementary Six-Element Figures", ""])
    for path in sorted(FIG.glob("nature_supp_six_elements_*.png")):
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    metrics, mueller, gate, selector, convergence = load_data()
    paths = []
    paths.extend(fig01_graphical_abstract(metrics, gate))
    paths.extend(fig02_case_library(metrics, gate, selector))
    paths.extend(fig03_accuracy_forest(metrics))
    paths.extend(fig04_mueller_accepted(mueller))
    paths.extend(fig05_stress_deviation(mueller, metrics))
    paths.extend(fig06_gate_diagnostics(metrics, gate, selector))
    paths.extend(fig07_convergence(convergence))
    paths.extend(fig08_observable_decomposition(metrics))
    paths.extend(fig09_api_reproducibility())
    paths.extend(fig10_manuscript_roadmap())
    paths.extend(supp_case_mueller(mueller))
    write_manifest(paths)
    for path in paths:
        print(f"Wrote {path}")
    print(f"Updated {REPORT}")


if __name__ == "__main__":
    main()
