from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "paper_validation" / "large_spheroid_lambda01"
DATA = BASE / "data"
FIG = BASE / "figures"
REPORT = BASE / "reports" / "large_spheroid_lambda01_figures.md"
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
    "EBCM n80": "#1C1C1C",
    "EBCM n70": "#777777",
    "PMM n60": "#149C88",
    "SVM n60": "#3B6EA8",
    "PMM n80": "#D96B27",
    "SVM n80": "#C76D9D",
}

STYLES = {
    "EBCM n80": dict(color=COLORS["EBCM n80"], lw=1.8, ls="-"),
    "EBCM n70": dict(color=COLORS["EBCM n70"], lw=1.2, ls=":"),
    "PMM n60": dict(color=COLORS["PMM n60"], lw=1.25, ls="-"),
    "SVM n60": dict(color=COLORS["SVM n60"], lw=1.25, ls=(0, (3, 1.7))),
    "PMM n80": dict(color=COLORS["PMM n80"], lw=1.0, ls=(0, (3, 1.5, 1, 1.5))),
    "SVM n80": dict(color=COLORS["SVM n80"], lw=1.0, ls=(0, (2, 1.5))),
}


def panel(ax, label):
    ax.text(
        -0.13,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
    )


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


def load():
    amp = pd.read_csv(DATA / "amplitude_matrix.csv")
    summary = pd.read_csv(DATA / "summary.csv")
    return amp, summary


def complex_col(frame, stem):
    return frame[f"{stem}_real"].to_numpy() + 1j * frame[f"{stem}_imag"].to_numpy()


def plot_amplitude_elements(amp):
    methods = ["EBCM n80", "EBCM n70", "PMM n60", "SVM n60"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), sharex=True)
    specs = [
        ("S11", "real", np.real, "Re S11"),
        ("S11", "imag", np.imag, "Im S11"),
        ("S22", "real", np.real, "Re S22"),
        ("S22", "imag", np.imag, "Im S22"),
    ]
    for ax, (stem, _, op, title) in zip(axes.ravel(), specs):
        for method in methods:
            sub = amp[amp["method"] == method].sort_values("angle_deg")
            z = complex_col(sub, stem)
            ax.plot(sub["angle_deg"], op(z), label=method, **STYLES[method])
        ax.set_title(title)
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
        ax.set_ylabel("amplitude")
        light_grid(ax)
    for ax in axes[1, :]:
        ax.set_xlabel("scattering angle (deg)")
    axes[0, 0].legend(frameon=False, ncol=2, loc="upper right")
    panel(axes[0, 0], "A")
    panel(axes[0, 1], "B")
    panel(axes[1, 0], "C")
    panel(axes[1, 1], "D")
    fig.suptitle("Large prolate spheroid: complex amplitude-matrix elements", y=0.995)
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.11, top=0.86, wspace=0.24, hspace=0.38)
    return save(fig, "large_spheroid_amplitude_elements")


def plot_comparison_diagnostics(amp, summary):
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.15, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    for method in ["EBCM n70", "PMM n60", "SVM n60", "PMM n80", "SVM n80"]:
        sub = amp[amp["method"] == method].sort_values("angle_deg")
        ax.plot(
            sub["angle_deg"],
            np.clip(sub["local_relerr_vs_ebcm80"], 1e-16, None),
            label=method,
            **STYLES[method],
        )
    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_ylim(1e-7, 1e2)
    ax.set_xlabel("scattering angle (deg)")
    ax.set_ylabel("local relative error")
    ax.set_title("Amplitude matrix error vs EBCM n80")
    light_grid(ax)
    ax.legend(frameon=False, ncol=2, loc="upper left")
    panel(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    methods = ["EBCM n70", "EBCM n80", "PMM n60", "SVM n60", "PMM n80", "SVM n80"]
    vals = [float(summary.loc[summary["method"] == m, "amplitude_relerr_vs_ebcm80"].iloc[0]) for m in methods]
    colors = [COLORS[m] for m in methods]
    ax.barh(np.arange(len(methods)), np.clip(vals, 1e-16, None), color=colors, height=0.58)
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    ax.invert_yaxis()
    ax.set_xlabel("global relative error")
    ax.set_title("Full-angle Frobenius error")
    light_grid(ax, "x")
    panel(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    methods2 = ["EBCM n70", "EBCM n80", "PMM n60", "SVM n60", "PMM n80", "SVM n80"]
    x = np.arange(len(methods2))
    csca = [float(summary.loc[summary["method"] == m, "Csca"].iloc[0]) for m in methods2]
    cext = [float(summary.loc[summary["method"] == m, "Cext"].iloc[0]) for m in methods2]
    width = 0.36
    ax.bar(x - width / 2, np.clip(csca, 1e-30, None), width, color="#3B6EA8", label="Csca")
    ax.bar(x + width / 2, np.clip(cext, 1e-30, None), width, color="#D96B27", label="Cext")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods2, rotation=30, ha="right")
    ax.set_ylabel("cross section")
    ax.set_title("Energy/passivity check")
    light_grid(ax, "y")
    ax.legend(frameon=False, ncol=2)
    panel(ax, "C")

    ax = fig.add_subplot(gs[1, 1])
    ebcm = amp[amp["method"] == "EBCM n80"].sort_values("angle_deg")
    s11 = complex_col(ebcm, "S11")
    s22 = complex_col(ebcm, "S22")
    ax.plot(ebcm["angle_deg"], np.abs(s11), color="#1C1C1C", lw=1.6, label="|S11|")
    ax.plot(ebcm["angle_deg"], np.abs(s22), color="#149C88", lw=1.4, label="|S22|")
    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_xlabel("scattering angle (deg)")
    ax.set_ylabel("EBCM amplitude magnitude")
    ax.set_title("Reference diagonal amplitudes")
    light_grid(ax)
    ax.legend(frameon=False)
    panel(ax, "D")

    fig.suptitle("EBCM comparison and PMM/SVM stability diagnostics", y=0.995)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.15, top=0.87, wspace=0.42, hspace=0.52)
    return save(fig, "large_spheroid_ebcm_diagnostics")


def write_report(paths):
    lines = [
        "# Large Spheroid Amplitude Figures",
        "",
        "Files generated from `large_spheroid_lambda01/data/*.csv`.",
        "",
    ]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    amp, summary = load()
    paths = []
    paths.extend(plot_amplitude_elements(amp))
    paths.extend(plot_comparison_diagnostics(amp, summary))
    write_report(paths)
    for path in paths:
        print(f"Wrote {path}")
    print(f"Updated {REPORT}")


if __name__ == "__main__":
    main()
