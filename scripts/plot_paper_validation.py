from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper_validation" / "data"
FIG = ROOT / "paper_validation" / "figures"
REPORT = ROOT / "paper_validation" / "reports" / "paper_validation_report.md"
FIG.mkdir(parents=True, exist_ok=True)

METHOD_STYLE = {
    "EBCM": ("#111111", "-", "EBCM"),
    "SVM": ("#1f77b4", "--", "SVM"),
    "PMM": ("#2ca02c", "-", "PMM"),
    "IITM coarse": ("#d95f02", "-.", "IITM coarse"),
    "IITM fine": ("#6c6c6c", ":", "IITM fine"),
    "Auto selector": ("#9467bd", (0, (4, 1, 1, 1)), "Auto selector"),
}


def ordered_cases(frame, group=None):
    sub = frame if group is None else frame[frame["group"] == group]
    return list(dict.fromkeys(sub["case"]))


def labels_for_cases(frame, cases):
    out = []
    for case in cases:
        out.append(str(frame[frame["case"] == case]["label"].iloc[0]))
    return out


def plot_mueller_error_bars(metrics):
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    cases = ordered_cases(metrics, "accepted")
    labels = labels_for_cases(metrics, cases)
    x = np.arange(len(cases))
    width = 0.15

    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    for i, method in enumerate(methods):
        vals = []
        for case in cases:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == method)]
            vals.append(float(row["mueller_relerr_vs_ebcm"].iloc[0]))
        color, _, label = METHOD_STYLE[method]
        ax.bar(x + (i - 2) * width, vals, width=width, color=color, label=label)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Relative Mueller error vs EBCM")
    ax.set_title("Accepted Validation Cases: SVM/PMM T-Matrix Accuracy")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(ncols=3, fontsize=8)
    path = FIG / "accepted_mueller_error_bars.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_stress_gate_bars(metrics, gate):
    cases = ordered_cases(metrics, "stress")
    labels = labels_for_cases(metrics, cases)
    x = np.arange(len(cases))
    width = 0.26
    methods = ["PMM", "IITM coarse", "Auto selector"]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax = axes[0]
    for i, method in enumerate(methods):
        vals = []
        for case in cases:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == method)]
            vals.append(float(row["mueller_relerr_vs_ebcm"].iloc[0]))
        color, _, label = METHOD_STYLE[method]
        ax.bar(x + (i - 1) * width, vals, width=width, color=color, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Relative Mueller error vs EBCM")
    ax.set_title("Stress Cases: Gate Prevents Blind PMM Use")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1]
    conv = [float(gate[gate["case"] == case]["convergence_error"].iloc[0]) for case in cases]
    colors = ["#c44e52" for _ in cases]
    ax.bar(x, conv, color=colors, width=0.55)
    ax.axhline(5e-3, color="#111111", linestyle="--", lw=1.2, label="gate threshold")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("PMM refined-grid convergence error")
    ax.set_title("Rejected by No-Reference PMM Gate")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    path = FIG / "stress_gate_summary.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_error_decomposition(metrics):
    cases = ordered_cases(metrics, "accepted")
    labels = labels_for_cases(metrics, cases)
    methods = ["SVM", "PMM"]
    columns = [
        ("tmatrix_relerr_vs_ebcm", "T-matrix relative error"),
        ("Csca_relerr_vs_ebcm", "Csca relative error"),
        ("Cext_relerr_vs_ebcm", "Cext relative error"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), sharex=True)
    x = np.arange(len(cases))
    width = 0.32
    for ax, (column, title) in zip(axes, columns):
        for i, method in enumerate(methods):
            vals = []
            for case in cases:
                row = metrics[(metrics["case"] == case) & (metrics["method"] == method)]
                vals.append(float(row[column].iloc[0]))
            color, _, label = METHOD_STYLE[method]
            ax.bar(x + (i - 0.5) * width, vals, width=width, color=color, label=label)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=22, ha="right")
    axes[0].set_ylabel("Relative error vs EBCM")
    axes[0].legend(fontsize=8)
    path = FIG / "accepted_error_decomposition.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_convergence(convergence):
    cases = ordered_cases(convergence)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.8), sharex=True)
    axes = axes.ravel()
    for ax, case in zip(axes, cases):
        label = convergence[convergence["case"] == case]["label"].iloc[0]
        for method in ["SVM", "PMM"]:
            sub = convergence[(convergence["case"] == case) & (convergence["method"] == method)].sort_values("nmax")
            color, linestyle, display = METHOD_STYLE[method]
            ax.semilogy(sub["nmax"], sub["mueller_relerr_vs_ebcm"],
                        color=color, linestyle=linestyle, marker="o", lw=1.5,
                        label=display)
        ax.set_title(label)
        ax.set_ylabel("Mueller error")
        ax.grid(True, alpha=0.25)
    for ax in axes[2:]:
        ax.set_xlabel("nmax")
    axes[0].legend(fontsize=8)
    fig.suptitle("SVM/PMM Convergence Against EBCM", y=0.995)
    path = FIG / "svm_pmm_convergence.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_residual_vs_error(metrics):
    keep = metrics[(metrics["method"].isin(["SVM", "PMM", "IITM coarse", "IITM fine"])) &
                   np.isfinite(metrics["boundary_residual"])]
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    markers = {
        "SVM": "o",
        "PMM": "s",
        "IITM coarse": "^",
        "IITM fine": "D",
    }
    for method, sub in keep.groupby("method"):
        color, _, label = METHOD_STYLE[method]
        ax.scatter(sub["boundary_residual"], sub["mueller_relerr_vs_ebcm"],
                   marker=markers[method], color=color, s=48, alpha=0.82,
                   label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("PMM boundary residual")
    ax.set_ylabel("Relative Mueller error vs EBCM")
    ax.set_title("Boundary Residual as a Diagnostic, Not a Standalone Error Estimator")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    path = FIG / "residual_vs_mueller_error.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_runtime(metrics):
    cases = ordered_cases(metrics, "accepted")
    labels = labels_for_cases(metrics, cases)
    methods = ["SVM", "PMM", "IITM coarse", "IITM fine"]
    x = np.arange(len(cases))
    width = 0.18
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    for i, method in enumerate(methods):
        vals = []
        for case in cases:
            row = metrics[(metrics["case"] == case) & (metrics["method"] == method)]
            vals.append(float(row["runtime_s"].iloc[0]))
        color, _, label = METHOD_STYLE[method]
        ax.bar(x + (i - 1.5) * width, vals, width=width, color=color, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("Runtime Snapshot for Accepted Cases")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(ncols=2, fontsize=8)
    path = FIG / "accepted_runtime_snapshot.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_six_elements(mueller, case):
    elements = ["F11", "F12", "F22", "F33", "F34", "F44"]
    methods = ["EBCM", "SVM", "PMM", "IITM coarse", "IITM fine", "Auto selector"]
    label = mueller[mueller["case"] == case]["label"].iloc[0]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8), sharex=True)
    axes = axes.ravel()
    for ax, element in zip(axes, elements):
        for method in methods:
            sub = mueller[(mueller["case"] == case) & (mueller["method"] == method)].sort_values("angle_deg")
            if sub.empty:
                continue
            color, linestyle, display = METHOD_STYLE[method]
            lw = 1.8 if method == "EBCM" else 1.25
            alpha = 1.0 if method in ("EBCM", "SVM", "PMM") else 0.85
            ax.plot(sub["angle_deg"], sub[element], color=color, linestyle=linestyle,
                    lw=lw, alpha=alpha, label=display)
        ax.set_ylabel(element)
        ax.grid(True, alpha=0.25)
    axes[0].legend(fontsize=7)
    for ax in axes[3:]:
        ax.set_xlabel("Scattering angle (deg)")
    fig.suptitle(f"{label}: Six Independent Mueller Matrix Elements", y=0.995)
    path = FIG / f"six_elements_{case}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_gate_scorecard(gate):
    cases = list(gate["case"])
    labels = list(gate["label"])
    accepted = gate["accepted"].astype(str).str.lower().eq("true").to_numpy()
    conv = gate["convergence_error"].to_numpy(dtype=float)
    residual = gate["boundary_residual"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    x = np.arange(len(cases))
    colors = ["#2ca02c" if ok else "#c44e52" for ok in accepted]
    axes[0].bar(x, conv, color=colors, width=0.58)
    axes[0].axhline(5e-3, color="#111111", linestyle="--", lw=1.2)
    axes[0].set_yscale("log")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=25, ha="right")
    axes[0].set_ylabel("PMM convergence error")
    axes[0].set_title("Gate Criterion 1: Refined PMM Stability")
    axes[0].grid(True, axis="y", alpha=0.25)

    axes[1].bar(x, residual, color=colors, width=0.58)
    axes[1].axhline(0.5, color="#111111", linestyle="--", lw=1.2)
    axes[1].set_yscale("log")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1].set_ylabel("Boundary residual")
    axes[1].set_title("Gate Criterion 2: Boundary Residual")
    axes[1].grid(True, axis="y", alpha=0.25)

    path = FIG / "pmm_gate_scorecard.png"
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fmt(x):
    return f"{float(x):.3e}"


def write_report(metrics, gate, selector, paths, six_paths):
    accepted_cases = ordered_cases(metrics, "accepted")
    stress_cases = ordered_cases(metrics, "stress")
    lines = [
        "# Paper Validation Report",
        "",
        "This report validates `BoundaryTMatrixMethods.jl` against `TransitionMatrices.jl` EBCM/IITM references.",
        "SVM and PMM both return `AxisymmetricTransitionMatrix` objects and can be used directly with Mueller-matrix and cross-section routines.",
        "",
        "## Main Result",
        "",
        "For gate-accepted cases, SVM/PMM reproduce EBCM Mueller curves with small relative errors; for stress cases, the PMM gate rejects unstable PMM targets and the auto selector falls back to IITM.",
        "",
        "## Accepted Cases",
        "",
        "| case | SVM Mueller err | PMM Mueller err | PMM T err | PMM Cext err | gate |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for case in accepted_cases:
        label = metrics[metrics["case"] == case]["label"].iloc[0]
        svm = metrics[(metrics["case"] == case) & (metrics["method"] == "SVM")].iloc[0]
        pmm = metrics[(metrics["case"] == case) & (metrics["method"] == "PMM")].iloc[0]
        g = gate[gate["case"] == case].iloc[0]
        lines.append(
            f"| {label} | {fmt(svm['mueller_relerr_vs_ebcm'])} | "
            f"{fmt(pmm['mueller_relerr_vs_ebcm'])} | "
            f"{fmt(pmm['tmatrix_relerr_vs_ebcm'])} | "
            f"{fmt(pmm['Cext_relerr_vs_ebcm'])} | {g['status']} |"
        )

    lines.extend([
        "",
        "## Stress Cases",
        "",
        "| case | PMM Mueller err | IITM coarse err | auto method | gate reason |",
        "|---|---:|---:|---|---|",
    ])
    for case in stress_cases:
        label = metrics[metrics["case"] == case]["label"].iloc[0]
        pmm = metrics[(metrics["case"] == case) & (metrics["method"] == "PMM")].iloc[0]
        iitm = metrics[(metrics["case"] == case) & (metrics["method"] == "IITM coarse")].iloc[0]
        sel = selector[selector["case"] == case].iloc[0]
        g = gate[gate["case"] == case].iloc[0]
        lines.append(
            f"| {label} | {fmt(pmm['mueller_relerr_vs_ebcm'])} | "
            f"{fmt(iitm['mueller_relerr_vs_ebcm'])} | {sel['selected_method']} | "
            f"{g['reasons']} |"
        )

    lines.extend(["", "## Figures", ""])
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    lines.extend(["", "## Six-Element Mueller Figures", ""])
    for path in six_paths:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    metrics = pd.read_csv(DATA / "validation_metrics.csv")
    mueller = pd.read_csv(DATA / "validation_mueller.csv")
    gate = pd.read_csv(DATA / "validation_gate.csv")
    selector = pd.read_csv(DATA / "validation_selector.csv")
    convergence = pd.read_csv(DATA / "validation_convergence.csv")

    paths = [
        plot_mueller_error_bars(metrics),
        plot_error_decomposition(metrics),
        plot_convergence(convergence),
        plot_gate_scorecard(gate),
        plot_stress_gate_bars(metrics, gate),
        plot_residual_vs_error(metrics),
        plot_runtime(metrics),
    ]
    six_paths = []
    for case in ordered_cases(mueller):
        six_paths.append(plot_six_elements(mueller, case))

    write_report(metrics, gate, selector, paths, six_paths)
    for path in paths + six_paths:
        print(f"Wrote {path}")
    print(f"Updated {REPORT}")


if __name__ == "__main__":
    main()
