from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
SCAN_DIR = (
    REPO_ROOT
    / "data"
    / "raw"
    / "uvvis"
    / "Scan - Lambda 1050 Friday, January 23, 2026 10_47 AM Arabian Standard Time"
)
OUT_DIR = REPO_ROOT / "outputs" / "uvvis"
BOUNDARIES = (379.5, 728.5)
FIT_WINDOWS = {
    379.5: ((360.0, 378.0), (381.0, 399.0)),
    728.5: ((710.0, 728.0), (729.0, 747.0)),
}


def read_percent_reflectance(filename: str) -> tuple[np.ndarray, np.ndarray]:
    path = SCAN_DIR / filename
    with path.open(newline="") as handle:
        reader = csv.reader(handle, skipinitialspace=True)
        next(reader)
        rows = [(float(row[0]), float(row[1])) for row in reader if len(row) >= 2]

    data = np.array(sorted(rows), dtype=float)
    return data[:, 0], data[:, 1]


def fit_at_boundary(
    wavelength: np.ndarray,
    values: np.ndarray,
    boundary: float,
) -> tuple[float, float]:
    left_window, right_window = FIT_WINDOWS[boundary]
    left_mask = (wavelength >= left_window[0]) & (wavelength <= left_window[1])
    right_mask = (wavelength >= right_window[0]) & (wavelength <= right_window[1])

    left_fit = np.polyfit(wavelength[left_mask], values[left_mask], 1)
    right_fit = np.polyfit(wavelength[right_mask], values[right_mask], 1)
    return np.polyval(left_fit, boundary), np.polyval(right_fit, boundary)


def normalize_and_stitch(
    wavelength: np.ndarray,
    raw_percent_r: np.ndarray,
    baseline_percent_r: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    normalized = raw_percent_r / baseline_percent_r * 100.0
    cleaned = normalized.copy()

    uv_pred, visible_pred = fit_at_boundary(wavelength, normalized, 379.5)
    visible_pred_left, nir_pred = fit_at_boundary(wavelength, normalized, 728.5)

    uv_scale = visible_pred / uv_pred
    nir_scale = visible_pred_left / nir_pred

    cleaned[wavelength <= BOUNDARIES[0]] *= uv_scale
    cleaned[wavelength >= BOUNDARIES[1]] *= nir_scale

    meta = {
        "uv_scale": uv_scale,
        "nir_scale": nir_scale,
        "uv_boundary_nm": BOUNDARIES[0],
        "nir_boundary_nm": BOUNDARIES[1],
    }
    return normalized, cleaned, meta


def save_cleaned_csv(
    wavelength: np.ndarray,
    baseline: np.ndarray,
    raw_curves: dict[str, np.ndarray],
    normalized_curves: dict[str, np.ndarray],
    cleaned_curves: dict[str, np.ndarray],
) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUT_DIR / "uvvis_cleaned_reflectance.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "wavelength_nm",
                "baseline_percent_reflectance",
                "S1_raw_percent_reflectance",
                "S2_raw_percent_reflectance",
                "S1_baseline_normalized_percent_reflectance",
                "S2_baseline_normalized_percent_reflectance",
                "S1_cleaned_percent_reflectance",
                "S2_cleaned_percent_reflectance",
            ]
        )
        for i, nm in enumerate(wavelength):
            writer.writerow(
                [
                    f"{nm:.1f}",
                    f"{baseline[i]:.6f}",
                    f"{raw_curves['S1'][i]:.6f}",
                    f"{raw_curves['S2'][i]:.6f}",
                    f"{normalized_curves['S1'][i]:.6f}",
                    f"{normalized_curves['S2'][i]:.6f}",
                    f"{cleaned_curves['S1'][i]:.6f}",
                    f"{cleaned_curves['S2'][i]:.6f}",
                ]
            )
    return output


def save_metadata(scales: dict[str, dict[str, float]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUT_DIR / "uvvis_cleaning_notes.txt"
    lines = [
        "UV-Vis cleaning workflow",
        "1. Baseline-normalized reflectance: sample / baseline * 100",
        "2. Corrected detector-switch discontinuities by piecewise scaling",
        "3. Middle segment (380-728 nm) used as the anchor region",
        "",
    ]

    for label, meta in scales.items():
        lines.extend(
            [
                f"{label}",
                f" UV-segment scale (<= {meta['uv_boundary_nm']:.1f} nm): {meta['uv_scale']:.6f}",
                f" NIR-segment scale (>= {meta['nir_boundary_nm']:.1f} nm): {meta['nir_scale']:.6f}",
                "",
            ]
        )

    output.write_text("\n".join(lines), encoding="ascii")
    return output


def style_axes(ax: plt.Axes, ylabel: str = "Reflectance (%)") -> None:
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(ylabel)
    ax.grid(True, color="#d9d9d9", linewidth=0.8, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4.5, width=1)


def make_plot(
    wavelength: np.ndarray,
    curves: dict[str, np.ndarray],
    stem: str,
    note: str | None = None,
    ylim: tuple[float, float] | None = None,
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
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    colors = {
        "Al-Al$_2$O$_3$-5 nm": "#0b3954",
        "Al-Al$_2$O$_3$-Si": "#c1121f",
    }

    for label, values in curves.items():
        ax.plot(wavelength, values, linewidth=2.3, color=colors[label], label=label)

    ax.set_xlim(200, 800)
    if ylim is None:
        ymin = min(np.min(v) for v in curves.values())
        ymax = max(np.max(v) for v in curves.values())
        ylim = (
            np.floor((ymin - 1.5) / 2) * 2,
            min(102, np.ceil((ymax + 1.0) / 2) * 2),
        )
    ax.set_ylim(*ylim)
    ax.set_yticks(np.arange(np.ceil(ylim[0] / 5) * 5, ylim[1] + 0.1, 5))
    ax.legend(frameon=False, loc="upper right", handlelength=2.4)
    style_axes(ax)

    if note:
        fig.text(0.125, 0.02, note, ha="left", va="bottom", fontsize=9, color="#4d4d4d")

    plt.tight_layout(rect=(0, 0.05 if note else 0, 1, 1))

    target = OUT_DIR / f"{stem}.png"
    fig.savefig(target, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [target]


def main() -> None:
    wavelength, baseline = read_percent_reflectance(
        "100% or 0 Absorbance Baseline.Correction.Raw.csv"
    )
    wavelength_s1, raw_s1 = read_percent_reflectance("S1.Sample.Raw.csv")
    wavelength_s2, raw_s2 = read_percent_reflectance("S2.Sample.Raw.csv")

    if not (
        np.array_equal(wavelength, wavelength_s1)
        and np.array_equal(wavelength, wavelength_s2)
    ):
        raise ValueError("Input scans do not share the same wavelength grid.")

    _norm_s1, cleaned_s1, _meta_s1 = normalize_and_stitch(wavelength, raw_s1, baseline)
    _norm_s2, cleaned_s2, _meta_s2 = normalize_and_stitch(wavelength, raw_s2, baseline)

    corrected_outputs = make_plot(
        wavelength,
        {
            "Al-Al$_2$O$_3$-5 nm": cleaned_s1,
            "Al-Al$_2$O$_3$-Si": cleaned_s2,
        },
        stem="uvvis_reflectance_corrected_publication_v2",
        note="Reflectance normalized to baseline and corrected for detector stitching",
        ylim=(64, 90),
    )
    for path in corrected_outputs:
        print(path)


if __name__ == "__main__":
    main()
