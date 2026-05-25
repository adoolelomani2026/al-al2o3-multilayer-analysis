from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data" / "raw" / "pl"
OUT_DIR = REPO_ROOT / "outputs" / "pl"
OMIT_WINDOW = (500.0, 620.0)


def read_scan(path: Path) -> tuple[dict[str, str], np.ndarray, np.ndarray]:
    metadata: dict[str, str] = {}
    rows: list[tuple[float, float]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            parts = [part.strip() for part in line.split(",")]
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    if parts[0] and parts[1]:
                        metadata[parts[0]] = parts[1]

    data = np.array(rows, dtype=float)
    return metadata, data[:, 0], data[:, 1]


def fit_baseline(
    x: np.ndarray,
    y: np.ndarray,
    exclude: tuple[float, float],
    degree: int = 3,
) -> np.ndarray:
    mask = (x < exclude[0]) | (x > exclude[1])
    fit_mask = mask.copy()

    for _ in range(6):
        coeff = np.polyfit(x[fit_mask], y[fit_mask], degree)
        baseline = np.polyval(coeff, x)
        residual = y - baseline
        fit_residual = residual[mask]
        mad = np.median(np.abs(fit_residual - np.median(fit_residual)))
        sigma = 1.4826 * mad if mad > 0 else max(np.std(fit_residual), 1.0)
        fit_mask = mask & (residual < 2.5 * sigma)

        if fit_mask.sum() < degree + 5:
            fit_mask = mask
            break

    coeff = np.polyfit(x[fit_mask], y[fit_mask], degree)
    return np.polyval(coeff, x)


def centered_moving_average(y: np.ndarray, window: int = 5) -> np.ndarray:
    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(y, (window // 2, window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def hampel_filter(y: np.ndarray, window: int = 7, n_sigma: float = 4.0) -> np.ndarray:
    out = y.copy()
    half = window // 2

    for i in range(len(y)):
        lo = max(0, i - half)
        hi = min(len(y), i + half + 1)
        segment = y[lo:hi]
        median = np.median(segment)
        mad = np.median(np.abs(segment - median))
        sigma = 1.4826 * mad
        if sigma > 0 and abs(y[i] - median) > n_sigma * sigma:
            out[i] = median

    return out


def process_emission_scan(
    path: Path,
) -> dict[str, np.ndarray | dict[str, str] | tuple[float, float]]:
    metadata, x, y = read_scan(path)
    exclude = OMIT_WINDOW
    display_ready = y.copy()

    display_mask = (x < exclude[0]) | (x > exclude[1])
    display_ready[display_mask] = hampel_filter(
        display_ready[display_mask],
        window=7,
        n_sigma=4.0,
    )

    baseline = fit_baseline(x, display_ready, exclude=exclude, degree=3)
    corrected = display_ready - baseline
    corrected = np.where(corrected > 0, corrected, 0.0)
    smooth = centered_moving_average(corrected, window=9)

    mask = (x < exclude[0]) | (x > exclude[1])
    norm = smooth.copy()
    norm_max = np.max(norm[mask]) if np.any(mask) else np.max(norm)
    if norm_max > 0:
        norm = norm / norm_max

    return {
        "metadata": metadata,
        "wavelength_nm": x,
        "raw_counts": y,
        "despiked_counts": display_ready,
        "baseline_counts": baseline,
        "corrected_counts": corrected,
        "smoothed_corrected_counts": smooth,
        "normalized_intensity": norm,
        "scatter_window_nm": exclude,
    }


def save_cleaned_csv(
    data: dict[str, dict[str, np.ndarray | tuple[float, float] | dict[str, str]]],
) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "pl_emission_cleaned.csv"
    x = data["scan1"]["wavelength_nm"]

    with out.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "wavelength_nm",
                "scan1_raw_counts",
                "scan1_despiked_counts",
                "scan1_baseline_counts",
                "scan1_corrected_counts",
                "scan1_smoothed_corrected_counts",
                "scan1_normalized_intensity",
                "scan2_raw_counts",
                "scan2_despiked_counts",
                "scan2_baseline_counts",
                "scan2_corrected_counts",
                "scan2_smoothed_corrected_counts",
                "scan2_normalized_intensity",
            ]
        )
        for i, wavelength in enumerate(x):
            writer.writerow(
                [
                    f"{wavelength:.1f}",
                    f"{data['scan1']['raw_counts'][i]:.6f}",
                    f"{data['scan1']['despiked_counts'][i]:.6f}",
                    f"{data['scan1']['baseline_counts'][i]:.6f}",
                    f"{data['scan1']['corrected_counts'][i]:.6f}",
                    f"{data['scan1']['smoothed_corrected_counts'][i]:.6f}",
                    f"{data['scan1']['normalized_intensity'][i]:.8f}",
                    f"{data['scan2']['raw_counts'][i]:.6f}",
                    f"{data['scan2']['despiked_counts'][i]:.6f}",
                    f"{data['scan2']['baseline_counts'][i]:.6f}",
                    f"{data['scan2']['corrected_counts'][i]:.6f}",
                    f"{data['scan2']['smoothed_corrected_counts'][i]:.6f}",
                    f"{data['scan2']['normalized_intensity'][i]:.8f}",
                ]
            )
    return out


def save_notes(
    data: dict[str, dict[str, np.ndarray | tuple[float, float] | dict[str, str]]],
) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "pl_processing_notes.txt"
    lines = [
        "PL emission processing workflow",
        "1. Parsed instrument-exported text/CSV scans directly",
        "2. Applied a Hampel despiking filter outside the excitation-dominated window",
        "3. Estimated a smooth polynomial baseline outside the excitation-dominated window",
        "4. Saved baseline-subtracted data in the cleaned CSV for reference",
        (
            "5. Final PL figure is diagnostic-style, not UV-Vis-style, because the scan is "
            "dominated by excitation/scatter artifacts"
        ),
        (
            "6. Top panel shows the raw primary emission scan; bottom panel shows off-band "
            "baseline-subtracted intensity only"
        ),
        "",
    ]

    for label, key in (("Scan 1", "scan1"), ("Scan 2", "scan2")):
        lo, hi = data[key]["scatter_window_nm"]
        lines.append(f"{label} omitted excitation-dominated window: {lo:.1f}-{hi:.1f} nm")

    out.write_text("\n".join(lines), encoding="ascii")
    return out


def make_plot(
    data: dict[str, dict[str, np.ndarray | tuple[float, float] | dict[str, str]]],
) -> list[Path]:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
        }
    )

    fig = plt.figure(figsize=(7.6, 6.2))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 0.95],
        width_ratios=[1.1, 1.0],
        hspace=0.28,
        wspace=0.08,
    )
    ax_top = fig.add_subplot(gs[0, :])
    ax_bot_left = fig.add_subplot(gs[1, 0])
    ax_bot_right = fig.add_subplot(gs[1, 1], sharey=ax_bot_left)
    fig.patch.set_facecolor("white")
    ax_top.set_facecolor("white")
    ax_bot_left.set_facecolor("white")
    ax_bot_right.set_facecolor("white")

    label = "Al-Al$_2$O$_3$ emission"
    x = data["scan1"]["wavelength_nm"]
    y_despiked = data["scan1"]["despiked_counts"].copy()
    y_clean = data["scan1"]["smoothed_corrected_counts"].copy()
    lo, hi = data["scan1"]["scatter_window_nm"]
    y_plot = np.where(y_despiked > 0.0, y_despiked, 0.0)
    y_clean = np.where(y_clean > 0.0, y_clean, 0.0)

    left_mask = x <= lo
    right_mask = x >= hi

    ax_top.plot(x, y_plot, linewidth=2.1, color="#0b3954", label=label)
    ax_top.axvspan(lo, hi, color="#d9d9d9", alpha=0.4, linewidth=0)
    ax_top.set_xlim(300, 800)
    ax_top.set_ylabel("Raw intensity (counts)")
    ax_top.legend(frameon=False, loc="upper right", handlelength=2.4)

    ax_bot_left.plot(x[left_mask], y_clean[left_mask], linewidth=2.1, color="#0b3954")
    ax_bot_right.plot(x[right_mask], y_clean[right_mask], linewidth=2.1, color="#0b3954")

    ymax = max(np.nanmax(y_clean[left_mask]), np.nanmax(y_clean[right_mask]))
    ax_bot_left.set_xlim(300, lo)
    ax_bot_right.set_xlim(hi, 800)
    ax_bot_left.set_ylim(0, ymax * 1.1)

    for ax in (ax_top, ax_bot_left, ax_bot_right):
        ax.grid(True, color="#d9d9d9", linewidth=0.8, alpha=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(direction="out", length=4.5, width=1)

    ax_bot_right.spines["left"].set_visible(False)
    ax_bot_right.tick_params(left=False, labelleft=False)
    ax_bot_left.set_ylabel("Off-band PL\n(counts)")
    fig.supxlabel("Emission wavelength (nm)", y=0.06)

    ax_top.text(
        (lo + hi) / 2,
        ax_top.get_ylim()[1] * 0.93,
        "Excitation-dominated region\n(500-620 nm) excluded",
        ha="center",
        va="top",
        fontsize=9,
        color="#555555",
    )

    d = 0.015
    kwargs_left = dict(
        transform=ax_bot_left.transAxes,
        color="#555555",
        clip_on=False,
        linewidth=1.0,
    )
    kwargs_right = dict(
        transform=ax_bot_right.transAxes,
        color="#555555",
        clip_on=False,
        linewidth=1.0,
    )
    ax_bot_left.plot((1 - d, 1 + d), (-d, +d), **kwargs_left)
    ax_bot_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs_left)
    ax_bot_right.plot((-d, +d), (-d, +d), **kwargs_right)
    ax_bot_right.plot((-d, +d), (1 - d, 1 + d), **kwargs_right)

    ax_top.annotate(
        "Raw scan shown as measured\n(light despiking only)",
        xy=(690, np.interp(690, x, y_plot)),
        xytext=(645, ax_top.get_ylim()[1] * 0.6),
        arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "#555555"},
        fontsize=8.5,
        color="#555555",
    )
    ax_bot_left.text(
        0.02,
        0.93,
        "Off-band PL (baseline-subtracted)",
        transform=ax_bot_left.transAxes,
        fontsize=9,
        color="#555555",
        va="top",
    )
    fig.text(
        0.125,
        0.01,
        (
            "Diagnostic PL figure: top = raw primary scan, bottom = off-band cleaned view; "
            "excitation-dominated region (500-620 nm) excluded"
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#4d4d4d",
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.14)

    target = OUT_DIR / "pl_emission_publication.png"
    fig.savefig(target, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [target]


def main() -> None:
    scan1 = process_emission_scan(DATA_DIR / "Emission ScanAl-Al2O3. Ex 565. BW 5.CSV")
    scan2 = process_emission_scan(DATA_DIR / "Emission ScanAl-Al2O3. Ex 565. BW 5_second.CSV")
    data = {"scan1": scan1, "scan2": scan2}

    outputs = make_plot(data)

    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
