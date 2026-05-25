from __future__ import annotations

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


HC_EV_NM = 1239.841984
REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = REPO_ROOT / "data" / "raw" / "cl" / "CL Spectrum_5kV_13nA_400 nm.msa"
OUTPUT_DIR = REPO_ROOT / "outputs" / "cl"


def style_axes(ax, x_label: str, y_label: str, y_max: float | None = None) -> None:
    ax.set_facecolor("white")
    ax.set_xlabel(x_label, fontsize=11, fontweight="semibold", labelpad=8)
    ax.set_ylabel(y_label, fontsize=11, fontweight="semibold", labelpad=8)
    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=4,
        width=0.9,
        labelsize=10,
        colors="#222222",
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="out",
        length=2.5,
        width=0.7,
        colors="#444444",
    )

    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_linewidth(1.0)
        ax.spines[side].set_color("#222222")

    ax.grid(True, which="major", color="#d9d9d9", linewidth=0.7, alpha=0.9)
    ax.grid(True, which="minor", color="#eeeeee", linewidth=0.5, alpha=0.9)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(x=0.01)

    if y_max is not None:
        ax.set_ylim(0, y_max)


def parse_msa(path: Path) -> tuple[dict[str, str], pd.DataFrame]:
    metadata: dict[str, str] = {}
    rows: list[tuple[float, float]] = []
    in_spectrum = False

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            if ":" in line:
                key, value = line.split(":", 1)
                metadata[key.strip()] = value.strip()
            if line.startswith("#SPECTRUM"):
                in_spectrum = True
            elif line.startswith("#ENDOFDATA"):
                in_spectrum = False
            continue

        if in_spectrum and "," in line:
            x_nm, counts = [float(part.strip()) for part in line.split(",", 1)]
            rows.append((x_nm, counts))

    df = pd.DataFrame(rows, columns=["wavelength_nm", "counts"])
    return metadata, df


def main() -> None:
    metadata, df = parse_msa(INPUT_FILE)
    if df.empty:
        raise ValueError(f"No spectrum data found in {INPUT_FILE}")

    df["energy_eV"] = HC_EV_NM / df["wavelength_nm"]
    df = df.sort_values("energy_eV").reset_index(drop=True)
    df["normalized_counts"] = df["counts"] / df["counts"].max()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    norm_plot_path = OUTPUT_DIR / "cl_spectrum_ev_normalized.png"

    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "semibold",
            "axes.labelcolor": "#222222",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=300)
    fig.patch.set_facecolor("white")
    ax.plot(
        df["energy_eV"],
        df["normalized_counts"],
        color="#0d7a64",
        linewidth=2.1,
        solid_capstyle="round",
    )
    style_axes(ax, "Energy (eV)", "Normalized Intensity (a.u.)", y_max=1.05)
    fig.tight_layout(pad=0.8)
    fig.savefig(norm_plot_path, bbox_inches="tight")
    plt.close(fig)

    peak_row = df.loc[df["counts"].idxmax()]
    print(f"Saved normalized plot: {norm_plot_path}")
    print(
        "Peak intensity: "
        f"{peak_row['counts']:.0f} counts at "
        f"{peak_row['wavelength_nm']:.2f} nm ({peak_row['energy_eV']:.3f} eV)"
    )


if __name__ == "__main__":
    main()
