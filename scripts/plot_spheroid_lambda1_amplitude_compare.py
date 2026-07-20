from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "paper_validation" / "spheroid_lambda1"
DATA = BASE / "data"
FIG = BASE / "figures"
REPORT = BASE / "reports" / "spheroid_lambda1_figures.md"
FIG.mkdir(parents=True, exist_ok=True)
REPORT.parent.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.0,
        "axes.titlesize": 9.2,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "legend.fontsize": 7.0,
        "figure.titlesize": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "savefig.dpi": 420,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

COLORS = {
    "EBCM n18": "#1C1C1C",
    "EBCM n16": "#777777",
    "EBCM n14": "#A9A9A9",
    "PMM n16": "#149C88",
    "SVM n16": "#3B6EA8",
}

STYLES = {
    "EBCM n18": dict(color=COLORS["EBCM n18"], lw=1.8, ls="-"),
    "EBCM n16": dict(color=COLORS["EBCM n16"], lw=1.2, ls=":"),
    "EBCM n14": dict(color=COLORS["EBCM n14"], lw=1.0, ls=(0, (3, 1.5))),
    "PMM n16": dict(color=COLORS["PMM n16"], lw=1.25, ls="-"),
    "SVM n16": dict(color=COLORS["SVM n16"], lw=1.25, ls=(0, (3, 1.7))),
}


def panel(ax, label):
    ax.text(-0.13, 1.05, label, transform=ax.transAxes, ha="left",
            va="top", fontsize=10, fontweight="bold")


def light_grid(ax, axis="both"):
    ax.grid(True, axis=axis, color="#ECECEC", lw=0.55)
    ax.set_axisbelow(True)


def save(fig, name):
    png = FIG / f"{name}.png"
    pdf = FIG / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


def complex_col(frame, stem):
    return frame[f"{stem}_real"].to_numpy() + 1j * frame[f"{stem}_imag"].to_numpy()


def plot_elements(amp):
    methods = ["EBCM n18", "EBCM n16", "PMM n16", "SVM n16"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), sharex=True)
    specs = [
        ("S11", np.real, "Re S11"),
        ("S11", np.imag, "Im S11"),
        ("S22", np.real, "Re S22"),
        ("S22", np.imag, "Im S22"),
    ]
    for ax, (stem, op, title) in zip(axes.ravel(), specs):
        for method in methods:
            sub = amp[amp["method"] == method].sort_values("angle_deg")
            ax.plot(sub["angle_deg"], op(complex_col(sub, stem)),
                    label=method, **STYLES[method])
        ax.set_title(title)
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
        ax.set_ylabel("amplitude")
        light_grid(ax)
    for ax in axes[1, :]:
        ax.set_xlabel("scattering angle (deg)")
    axes[0, 0].legend(frameon=False, ncol=2, loc="upper right")
    for label, ax in zip("ABCD", axes.ravel()):
        panel(ax, label)
    fig.suptitle("Spheroid at lambda=1: complex amplitude-matrix elements", y=0.995)
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.11, top=0.86, wspace=0.24, hspace=0.38)
    return save(fig, "spheroid_lambda1_amplitude_elements")


def plot_diagnostics(amp, summary):
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.15, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    for method in ["EBCM n14", "EBCM n16", "PMM n16", "SVM n16"]:
        sub = amp[amp["method"] == method].sort_values("angle_deg")
        ax.plot(sub["angle_deg"], np.clip(sub["local_relerr_vs_ebcm18"], 1e-16, None),
                label=method, **STYLES[method])
    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_ylim(1e-10, 1)
    ax.set_xlabel("scattering angle (deg)")
    ax.set_ylabel("local relative error")
    ax.set_title("Amplitude matrix error vs EBCM n18")
    light_grid(ax)
    ax.legend(frameon=False, ncol=2, loc="upper left")
    panel(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    methods = ["EBCM n14", "EBCM n16", "EBCM n18", "PMM n16", "SVM n16"]
    vals = [float(summary.loc[summary["method"] == m, "amplitude_relerr_vs_ebcm18"].iloc[0])
            for m in methods]
    ax.barh(np.arange(len(methods)), np.clip(vals, 1e-16, None),
            color=[COLORS[m] for m in methods], height=0.58)
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    ax.invert_yaxis()
    ax.set_xlabel("global relative error")
    ax.set_title("Full-angle Frobenius error")
    light_grid(ax, "x")
    panel(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    x = np.arange(len(methods))
    csca = [float(summary.loc[summary["method"] == m, "Csca_relerr_vs_ebcm18"].iloc[0])
            for m in methods]
    cext = [float(summary.loc[summary["method"] == m, "Cext_relerr_vs_ebcm18"].iloc[0])
            for m in methods]
    width = 0.36
    ax.bar(x - width / 2, np.clip(csca, 1e-16, None), width, color="#3B6EA8", label="Csca err")
    ax.bar(x + width / 2, np.clip(cext, 1e-16, None), width, color="#D96B27", label="Cext err")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=30, ha="right")
    ax.set_ylabel("relative error")
    ax.set_title("Cross-section errors")
    light_grid(ax, "y")
    ax.legend(frameon=False)
    panel(ax, "C")

    ax = fig.add_subplot(gs[1, 1])
    ebcm = amp[amp["method"] == "EBCM n18"].sort_values("angle_deg")
    ax.plot(ebcm["angle_deg"], np.abs(complex_col(ebcm, "S11")),
            color="#1C1C1C", lw=1.6, label="|S11|")
    ax.plot(ebcm["angle_deg"], np.abs(complex_col(ebcm, "S22")),
            color="#149C88", lw=1.4, label="|S22|")
    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_xlabel("scattering angle (deg)")
    ax.set_ylabel("EBCM amplitude magnitude")
    ax.set_title("Reference diagonal amplitudes")
    light_grid(ax)
    ax.legend(frameon=False)
    panel(ax, "D")

    fig.suptitle("Lambda=1 EBCM comparison and PMM/SVM diagnostics", y=0.995)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.15, top=0.87, wspace=0.42, hspace=0.52)
    return save(fig, "spheroid_lambda1_ebcm_diagnostics")


def write_report(paths):
    lines = ["# Spheroid Lambda=1 Amplitude Figures", ""]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    amp = pd.read_csv(DATA / "amplitude_matrix.csv")
    summary = pd.read_csv(DATA / "summary.csv")
    paths = []
    paths.extend(plot_elements(amp))
    paths.extend(plot_diagnostics(amp, summary))
    write_report(paths)
    for path in paths:
        print(f"Wrote {path}")
    print(f"Updated {REPORT}")


if __name__ == "__main__":
    main()
