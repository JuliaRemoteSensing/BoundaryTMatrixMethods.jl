from pathlib import Path
from textwrap import fill

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper_validation" / "data"
FIG = ROOT / "paper_validation" / "figures_publication"
REPORT = ROOT / "paper_validation" / "reports" / "publication_figure_manifest.md"
FIG.mkdir(parents=True, exist_ok=True)
REPORT.parent.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "figure.titlesize": 13,
        "figure.dpi": 160,
        "savefig.dpi": 360,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.7,
        "lines.markersize": 4.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

COLORS = {
    "EBCM": "#111111",
    "SVM": "#0072B2",
    "PMM": "#009E73",
    "IITM coarse": "#D55E00",
    "IITM fine": "#6F6F6F",
    "Auto selector": "#CC79A7",
    "accepted": "#009E73",
    "stress": "#D55E00",
    "background": "#F7F7F2",
    "ink": "#202020",
}

METHOD_STYLE = {
    "EBCM": dict(color=COLORS["EBCM"], linestyle="-", linewidth=2.2, zorder=5),
    "SVM": dict(color=COLORS["SVM"], linestyle=(0, (4, 2)), linewidth=1.8),
    "PMM": dict(color=COLORS["PMM"], linestyle="-", linewidth=1.8),
    "IITM coarse": dict(color=COLORS["IITM coarse"], linestyle=(0, (4, 2, 1, 2)), linewidth=1.4),
    "IITM fine": dict(color=COLORS["IITM fine"], linestyle=":", linewidth=1.7),
    "Auto selector": dict(color=COLORS["Auto selector"], linestyle=(0, (2, 1)), linewidth=1.8),
}

SHORT_LABELS = {
    "sphere": "Sphere",
    "prolate_mild": "Mild\nprolate",
    "oblate_mild": "Mild\noblate",
    "chebyshev_n2": "Cheb.\nn=2",
    "prolate_strong": "Strong\nprolate",
    "oblate_strong": "Strong\noblate",
    "chebyshev_n4_pos": "Cheb.\nn=4+",
    "chebyshev_n4_neg": "Cheb.\nn=4-",
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


def load_data():
    metrics = pd.read_csv(DATA / "validation_metrics.csv")
    mueller = pd.read_csv(DATA / "validation_mueller.csv")
    gate = pd.read_csv(DATA / "validation_gate.csv")
    selector = pd.read_csv(DATA / "validation_selector.csv")
    convergence = pd.read_csv(DATA / "validation_convergence.csv")
    return metrics, mueller, gate, selector, convergence


def save(fig, stem):
    png = FIG / f"{stem}.png"
    pdf = FIG / f"{stem}.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


def ordered_cases(frame, group=None):
    sub = frame if group is None else frame[frame["group"] == group]
    return list(dict.fromkeys(sub["case"]))


def setup_log_axis(ax, which="y"):
    ax.grid(True, axis=which, color="#D9D9D9", linewidth=0.65, alpha=0.75)
    ax.set_axisbelow(True)


def shape_profile(case, n=220):
    kind, a, c, order, eps = CASE_PARAMS[case]
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    if kind == "spheroid":
        x = a * np.sin(theta)
        z = c * np.cos(theta)
    else:
        r = 1.0 + eps * np.cos(order * theta)
        x = r * np.sin(theta)
        z = r * np.cos(theta)
    return x, z


def add_panel_label(ax, label):
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def fig01_workflow():
    fig, ax = plt.subplots(figsize=(12.4, 6.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x, y, w, h, text, color, lw=1.4):
        rect = patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            facecolor=color,
            edgecolor="#2B2B2B",
            linewidth=lw,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="#1D1D1D")
        return rect

    def arrow(xy0, xy1, color="#333333", lw=1.4):
        ax.annotate(
            "",
            xy=xy1,
            xytext=xy0,
            arrowprops=dict(arrowstyle="-|>", lw=lw, color=color, shrinkA=4, shrinkB=4),
        )

    b0 = box(0.04, 0.64, 0.18, 0.17, "Particle shape\nand material", "#F0F3F7")
    b1 = box(0.29, 0.64, 0.23, 0.17, "Axisymmetric\nboundary system", "#EFF7F3")
    b2 = box(0.60, 0.72, 0.18, 0.13, "SVM\nsurface solve", "#EAF3FB")
    b3 = box(0.60, 0.52, 0.18, 0.13, "PMM\nleast squares", "#EAF8F1")
    b4 = box(0.83, 0.62, 0.13, 0.16, "T matrix\nobject", "#F5F0FA")
    b5 = box(0.32, 0.28, 0.25, 0.15, "PMM residual gate\nstability + residual + passivity", "#FFF4E5")
    b6 = box(0.66, 0.28, 0.24, 0.15, "Auto selector\nPMM if safe, IITM otherwise", "#F7F1EE")
    b7 = box(0.38, 0.06, 0.30, 0.13, "Mueller curves,\ncross sections, validation data", "#F4F4F4")

    arrow((0.22, 0.725), (0.29, 0.725))
    arrow((0.52, 0.73), (0.60, 0.78))
    arrow((0.52, 0.69), (0.60, 0.58))
    arrow((0.78, 0.785), (0.83, 0.70))
    arrow((0.78, 0.585), (0.83, 0.70))
    arrow((0.69, 0.52), (0.48, 0.43), COLORS["PMM"])
    arrow((0.57, 0.355), (0.66, 0.355))
    arrow((0.78, 0.43), (0.86, 0.62))
    arrow((0.49, 0.28), (0.53, 0.19))
    arrow((0.78, 0.28), (0.60, 0.19))
    arrow((0.895, 0.62), (0.62, 0.19))

    ax.text(
        0.50,
        0.93,
        "BoundaryTMatrixMethods.jl validation pipeline",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        0.50,
        0.885,
        "SVM and PMM generate TransitionMatrices-compatible T matrices; a no-reference gate prevents unsafe PMM use.",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#555555",
    )
    return save(fig, "fig01_workflow_schematic")


def fig02_benchmark_landscape(metrics, gate):
    cases = ordered_cases(metrics)
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    fig = plt.figure(figsize=(13.2, 8.4), constrained_layout=True)
    gs = fig.add_gridspec(3, 8, height_ratios=[1.0, 0.08, 2.1], hspace=0.08, wspace=0.12)

    for i, case in enumerate(cases):
        ax = fig.add_subplot(gs[0, i])
        x, z = shape_profile(case)
        group = gate.loc[gate["case"] == case, "group"].iloc[0]
        accepted = str(gate.loc[gate["case"] == case, "accepted"].iloc[0]).lower() == "true"
        edge = COLORS["accepted"] if accepted else COLORS["stress"]
        ax.fill(x, z, facecolor="#F2F2F2", edgecolor=edge, linewidth=2.0)
        ax.axhline(0, color="#D7D7D7", lw=0.5)
        ax.axvline(0, color="#D7D7D7", lw=0.5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(SHORT_LABELS[case], fontsize=9.2, pad=2)
        if i == 0:
            add_panel_label(ax, "A")
        ax.text(
            0.5,
            -0.10,
            "gate pass" if accepted else "gate reject",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=7.8,
            color=edge,
        )

    ax = fig.add_subplot(gs[2, :])
    x = np.arange(len(cases))
    width = 0.145
    for i, method in enumerate(methods):
        vals = [
            float(
                metrics.loc[
                    (metrics["case"] == case) & (metrics["method"] == method),
                    "mueller_relerr_vs_ebcm",
                ].iloc[0]
            )
            for case in cases
        ]
        ax.bar(
            x + (i - 2) * width,
            vals,
            width=width,
            color=COLORS[method],
            label=method,
            alpha=0.94,
        )

    for i, case in enumerate(cases):
        group = gate.loc[gate["case"] == case, "group"].iloc[0]
        if group == "stress":
            ax.axvspan(i - 0.5, i + 0.5, color="#FFF0E7", zorder=-5)
    ax.set_yscale("log")
    ax.set_ylim(5e-16, 6e-2)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_LABELS[c].replace("\n", " ") for c in cases], rotation=18, ha="right")
    ax.set_ylabel("Relative Mueller-matrix error vs EBCM")
    ax.set_title("Validation landscape across accepted and stress cases", pad=10)
    setup_log_axis(ax)
    ax.legend(ncol=5, loc="upper left", frameon=False)
    ax.text(
        4.55,
        2.6e-15,
        "orange background: stress cases\nselector falls back to IITM",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#6A3A19",
    )
    add_panel_label(ax, "B")
    return save(fig, "fig02_validation_landscape")


def fig03_accepted_mueller_atlas(mueller):
    cases = ordered_cases(mueller, "accepted")
    elements = ["F11", "F12", "F22", "F33", "F34", "F44"]
    methods = ["EBCM", "SVM", "PMM", "IITM fine"]
    fig, axes = plt.subplots(len(cases), len(elements), figsize=(14.4, 9.4), sharex=True)
    for r, case in enumerate(cases):
        for c, element in enumerate(elements):
            ax = axes[r, c]
            for method in methods:
                sub = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values(
                    "angle_deg"
                )
                ax.plot(sub["angle_deg"], sub[element], label=method, **METHOD_STYLE[method])
            if r == 0:
                ax.set_title(element)
            if c == 0:
                ax.annotate(
                    SHORT_LABELS[case].replace("\n", " "),
                    xy=(-0.30, 0.5),
                    xycoords="axes fraction",
                    ha="center",
                    va="center",
                    rotation=90,
                    fontsize=10.0,
                    color=COLORS["ink"],
                )
            ax.grid(True, color="#E2E2E2", lw=0.55, alpha=0.8)
            ax.set_xlim(0, 180)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="upper center", frameon=False, bbox_to_anchor=(0.5, 0.945))
    fig.suptitle("Six independent Mueller-matrix elements for gate-accepted cases", y=0.99)
    fig.text(0.50, 0.025, "Scattering angle (deg)", ha="center")
    fig.text(0.008, 0.50, "Mueller element value", va="center", rotation=90)
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.075, top=0.875, wspace=0.25, hspace=0.22)
    return save(fig, "fig03_accepted_mueller_atlas")


def fig04_stress_gate_fallback(metrics, gate, selector):
    stress = ordered_cases(metrics, "stress")
    all_cases = ordered_cases(metrics)
    fig = plt.figure(figsize=(13.2, 8.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(stress))
    width = 0.22
    for i, method in enumerate(["PMM", "IITM coarse", "Auto selector"]):
        vals = [
            float(
                metrics.loc[
                    (metrics["case"] == case) & (metrics["method"] == method),
                    "mueller_relerr_vs_ebcm",
                ].iloc[0]
            )
            for case in stress
        ]
        ax.bar(x + (i - 1) * width, vals, width=width, color=COLORS[method], label=method)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_LABELS[c].replace("\n", " ") for c in stress], rotation=15, ha="right")
    ax.set_ylabel("Relative Mueller error vs EBCM")
    ax.set_title("Blind PMM is unsafe for stress cases")
    setup_log_axis(ax)
    ax.legend(frameon=False, ncol=3)
    add_panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    conv = gate.set_index("case").loc[all_cases, "convergence_error"].astype(float)
    colors = [COLORS["accepted"] if c in ordered_cases(metrics, "accepted") else COLORS["stress"] for c in all_cases]
    ax.barh(np.arange(len(all_cases)), conv.values, color=colors, alpha=0.95)
    ax.axvline(5e-3, color="#202020", linestyle="--", linewidth=1.25, label="threshold")
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(all_cases)))
    ax.set_yticklabels([SHORT_LABELS[c].replace("\n", " ") for c in all_cases])
    ax.invert_yaxis()
    ax.set_xlabel("PMM refined-grid change")
    ax.set_title("No-reference gate criterion")
    ax.grid(True, axis="x", color="#D9D9D9", linewidth=0.65)
    ax.legend(frameon=False)
    add_panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    base = gate.set_index("case")
    offsets = {
        "sphere": (8, -2),
        "prolate_mild": (-78, 8),
        "oblate_mild": (-28, 12),
        "chebyshev_n2": (7, -8),
        "prolate_strong": (8, 11),
        "oblate_strong": (8, -10),
        "chebyshev_n4_pos": (-72, -2),
        "chebyshev_n4_neg": (-72, 12),
    }
    for case in all_cases:
        row = base.loc[case]
        accepted = str(row["accepted"]).lower() == "true"
        color = COLORS["accepted"] if accepted else COLORS["stress"]
        pmm_err = float(
            metrics.loc[
                (metrics["case"] == case) & (metrics["method"] == "PMM"),
                "mueller_relerr_vs_ebcm",
            ].iloc[0]
        )
        ax.scatter(
            float(row["convergence_error"]),
            float(row["boundary_residual"]),
            s=80 + 550 * min(pmm_err / 0.02, 1.0),
            color=color,
            edgecolor="#222222",
            linewidth=0.65,
            alpha=0.88,
        )
        ax.annotate(
            SHORT_LABELS[case].replace("\n", " "),
            xy=(float(row["convergence_error"]), float(row["boundary_residual"])),
            xytext=offsets.get(case, (6, 6)),
            textcoords="offset points",
            fontsize=7.5,
            color="#303030",
            arrowprops=dict(arrowstyle="-", color="#777777", lw=0.55, shrinkA=2, shrinkB=2),
        )
    ax.axvline(5e-3, color="#202020", linestyle="--", linewidth=1.15)
    ax.axhline(0.5, color="#202020", linestyle=":", linewidth=1.15)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("PMM refined-grid change")
    ax.set_ylabel("Boundary residual")
    ax.set_title("Gate map: stable PMM candidates stay left")
    ax.grid(True, color="#D9D9D9", linewidth=0.65, alpha=0.75)
    add_panel_label(ax, "C")

    ax = fig.add_subplot(gs[1, 1])
    selected = selector.set_index("case").loc[all_cases, "selected_method"]
    y = np.arange(len(all_cases))
    vals = [1 if selected[c] == "pmm" else 0 for c in all_cases]
    colors = [COLORS["PMM"] if v == 1 else COLORS["IITM coarse"] for v in vals]
    ax.barh(y, np.ones(len(all_cases)), color=colors, alpha=0.92)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT_LABELS[c].replace("\n", " ") for c in all_cases])
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.invert_yaxis()
    for i, v in enumerate(vals):
        ax.text(0.5, i, "PMM" if v else "IITM fallback", ha="center", va="center", color="white")
    ax.set_title("Automatic selector decision")
    for spine in ax.spines.values():
        spine.set_visible(False)
    add_panel_label(ax, "D")
    return save(fig, "fig04_stress_gate_fallback")


def fig05_error_decomposition(metrics):
    cases = ordered_cases(metrics)
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    columns = [
        ("mueller_relerr_vs_ebcm", "Mueller"),
        ("tmatrix_relerr_vs_ebcm", "T matrix"),
        ("Csca_relerr_vs_ebcm", "Csca"),
        ("Cext_relerr_vs_ebcm", "Cext"),
        ("boundary_residual", "Boundary residual"),
    ]

    matrix = []
    ylabels = []
    for method in methods:
        vals = []
        for case in cases:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == method)].iloc[0]
            vals.append(float(max(row["mueller_relerr_vs_ebcm"], 1e-16)))
        matrix.append(vals)
        ylabels.append(method)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6), constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(np.log10(np.asarray(matrix)), aspect="auto", cmap="viridis_r", vmin=-15, vmax=-1)
    ax.set_xticks(np.arange(len(cases)))
    ax.set_xticklabels([SHORT_LABELS[c].replace("\n", " ") for c in cases], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(ylabels)
    ax.set_title("Mueller error heat map")
    cb = fig.colorbar(im, ax=ax, shrink=0.92)
    cb.set_label("log10(relative error)")
    add_panel_label(ax, "A")

    ax = axes[1]
    accepted = ordered_cases(metrics, "accepted")
    x = np.arange(len(accepted))
    width = 0.18
    for i, (col, label) in enumerate(columns[1:4]):
        vals = []
        for case in accepted:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == "PMM")].iloc[0]
            vals.append(float(max(row[col], 1e-16)))
        ax.bar(x + (i - 1) * width, vals, width=width, label=label, alpha=0.92)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_LABELS[c].replace("\n", " ") for c in accepted], rotation=20, ha="right")
    ax.set_ylabel("Relative value vs EBCM")
    ax.set_title("Accepted PMM: errors by observable")
    setup_log_axis(ax)
    ax.legend(frameon=False, ncol=3)
    add_panel_label(ax, "B")
    return save(fig, "fig05_error_decomposition")


def fig06_convergence_runtime(metrics, convergence):
    cases = ordered_cases(convergence)
    fig = plt.figure(figsize=(13.2, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.05])

    for i, case in enumerate(cases):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        for method in ["SVM", "PMM"]:
            sub = convergence[(convergence["case"] == case) & (convergence["method"] == method)]
            ax.semilogy(
                sub["nmax"],
                sub["mueller_relerr_vs_ebcm"],
                marker="o",
                label=method,
                color=COLORS[method],
                linestyle=METHOD_STYLE[method]["linestyle"],
            )
        ax.set_title(SHORT_LABELS[case].replace("\n", " "))
        ax.set_xlabel("nmax")
        ax.set_ylabel("Mueller error")
        ax.grid(True, color="#D9D9D9", linewidth=0.65, alpha=0.75)
        if i == 0:
            ax.legend(frameon=False)
            add_panel_label(ax, "A")

    ax = fig.add_subplot(gs[:, 2])
    accepted = ordered_cases(metrics, "accepted")
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine"]
    x = np.arange(len(accepted))
    width = 0.18
    for i, method in enumerate(methods):
        vals = []
        for case in accepted:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == method)].iloc[0]
            vals.append(float(row["runtime_s"]))
        ax.bar(x + (i - 1.5) * width, vals, width=width, color=COLORS[method], label=method)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_LABELS[c].replace("\n", " ") for c in accepted], rotation=25, ha="right")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("Runtime snapshot")
    setup_log_axis(ax)
    ax.legend(frameon=False)
    add_panel_label(ax, "B")
    return save(fig, "fig06_convergence_runtime")


def fig07_case_closeups(mueller):
    saved = []
    elements = ["F11", "F12", "F22", "F33", "F34", "F44"]
    methods = ["EBCM", "SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    for case in ordered_cases(mueller):
        label = SHORT_LABELS[case].replace("\n", " ")
        group = mueller.loc[mueller["case"] == case, "group"].iloc[0]
        fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8), sharex=True, constrained_layout=True)
        axes = axes.ravel()
        for ax, element in zip(axes, elements):
            for method in methods:
                sub = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values(
                    "angle_deg"
                )
                alpha = 0.55 if method in ["IITM coarse", "IITM fine", "Auto selector"] else 1.0
                style = METHOD_STYLE[method].copy()
                style["alpha"] = alpha
                ax.plot(sub["angle_deg"], sub[element], label=method, **style)
            ax.set_title(element)
            ax.grid(True, color="#E2E2E2", lw=0.55, alpha=0.8)
            ax.set_xlim(0, 180)
        for ax in axes[3:]:
            ax.set_xlabel("Angle (deg)")
        axes[0].set_ylabel("Value")
        axes[3].set_ylabel("Value")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, ncol=6, loc="upper center", frameon=False, bbox_to_anchor=(0.5, 1.03))
        status = "gate accepted" if group == "accepted" else "stress / gate rejected"
        fig.suptitle(f"{label}: six Mueller-matrix elements ({status})", y=1.08)
        saved.extend(save(fig, f"supp_six_elements_{case}"))
    return saved


def write_manifest(paths):
    lines = [
        "# Publication Figure Manifest",
        "",
        "These figures are generated from the machine-readable validation CSV files.",
        "Each figure is saved as both PNG and PDF in `paper_validation/figures_publication`.",
        "",
        "## Main Figures",
        "",
    ]
    descriptions = {
        "fig01_workflow_schematic": "Workflow of SVM/PMM T-matrix generation, residual gate, and IITM fallback.",
        "fig02_validation_landscape": "Shape gallery and relative Mueller error overview across all cases.",
        "fig03_accepted_mueller_atlas": "Six independent Mueller elements for all gate-accepted validation cases.",
        "fig04_stress_gate_fallback": "Stress-case failure analysis and automatic selector decisions.",
        "fig05_error_decomposition": "Error heat map and accepted-case PMM observable decomposition.",
        "fig06_convergence_runtime": "Convergence behavior and runtime snapshot.",
    }
    for stem, desc in descriptions.items():
        lines.append(f"- `{(FIG / (stem + '.png')).relative_to(ROOT).as_posix()}` - {desc}")
        lines.append(f"- `{(FIG / (stem + '.pdf')).relative_to(ROOT).as_posix()}`")
    lines.extend(["", "## Supplementary Mueller Closeups", ""])
    for path in sorted(FIG.glob("supp_six_elements_*.png")):
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    metrics, mueller, gate, selector, convergence = load_data()
    paths = []
    paths.extend(fig01_workflow())
    paths.extend(fig02_benchmark_landscape(metrics, gate))
    paths.extend(fig03_accepted_mueller_atlas(mueller))
    paths.extend(fig04_stress_gate_fallback(metrics, gate, selector))
    paths.extend(fig05_error_decomposition(metrics))
    paths.extend(fig06_convergence_runtime(metrics, convergence))
    paths.extend(fig07_case_closeups(mueller))
    write_manifest(paths)
    for path in paths:
        print(f"Wrote {path}")
    print(f"Updated {REPORT}")


if __name__ == "__main__":
    main()
