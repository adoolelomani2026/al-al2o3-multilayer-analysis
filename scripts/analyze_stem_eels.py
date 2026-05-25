from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data" / "raw" / "eels"
OUTPUT_DIR = REPO_ROOT / "outputs" / "eels"


def load_curve(filename: str) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(DATA_DIR / filename, delimiter=",")
    return data[:, 0], data[:, 1]


def main() -> None:
    energy_eps1, eps1 = load_curve("Epsilon 1.csv")
    energy_eps2, eps2 = load_curve("Epsilon 2.csv")
    energy_loss, loss = load_curve("Energy-loss function.csv")

    mask = (energy_eps1 >= 0) & (energy_eps1 <= 30)
    mask_loss = (energy_loss >= 0) & (energy_loss <= 30)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(9, 8.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1.0, 1.0], "hspace": 0.08},
    )

    colors = {
        "eps1": "#0f4c5c",
        "eps2": "#c1121f",
        "loss": "#2b9348",
        "accent": "#d4a373",
    }

    axes[0].plot(energy_eps1[mask], eps1[mask], color=colors["eps1"], linewidth=2.2)
    axes[0].axhline(0, color="#555555", linewidth=1, alpha=0.8)
    axes[0].set_ylabel(r"$\varepsilon_1$")
    axes[0].set_title("Dielectric Reconstruction and Energy-Loss Response (EELS)", pad=12)

    sign = np.sign(eps1[mask])
    crossings = np.where(np.diff(sign) != 0)[0]
    if len(crossings):
        idx = crossings[np.argmin(np.abs(energy_eps1[mask][crossings] - 15))]
        x0 = energy_eps1[mask][idx]
        x1 = energy_eps1[mask][idx + 1]
        y0 = eps1[mask][idx]
        y1 = eps1[mask][idx + 1]
        crossing_energy = x0 - y0 * (x1 - x0) / (y1 - y0)
        axes[0].axvline(
            crossing_energy,
            color=colors["accent"],
            linestyle="--",
            linewidth=1.5,
        )
        axes[0].annotate(
            r"$\varepsilon_1 = 0$",
            xy=(crossing_energy, 0),
            xytext=(crossing_energy + 1.8, np.nanmax(eps1[mask]) * 0.58),
            arrowprops={"arrowstyle": "->", "color": colors["accent"], "lw": 1.1},
            color=colors["accent"],
            fontsize=10,
        )

    axes[1].plot(energy_eps2[mask], eps2[mask], color=colors["eps2"], linewidth=2.2)
    axes[1].set_ylabel(r"$\varepsilon_2$")

    peak_idx = np.argmax(loss[mask_loss])
    peak_energy = energy_loss[mask_loss][peak_idx]
    peak_value = loss[mask_loss][peak_idx]

    axes[2].plot(energy_loss[mask_loss], loss[mask_loss], color=colors["loss"], linewidth=2.2)
    axes[2].scatter([peak_energy], [peak_value], color=colors["accent"], s=35, zorder=3)
    axes[2].annotate(
        "plasmon peak\n~15 eV",
        xy=(peak_energy, peak_value),
        xytext=(peak_energy + 1.2, peak_value * 0.82),
        arrowprops={"arrowstyle": "->", "color": colors["accent"], "lw": 1.1},
        color=colors["accent"],
        fontsize=10,
    )
    axes[2].set_ylabel(r"Im$[-1/\varepsilon]$")
    axes[2].set_xlabel("Energy loss (eV)")

    uv_start = 3.1
    uv_end = 6.2
    for ax in axes:
        ax.axvspan(uv_start, uv_end, color="#f6bd60", alpha=0.12, linewidth=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].text(
        uv_start + 0.15,
        axes[0].get_ylim()[1] * 0.83,
        "Near-UV\n(200-400 nm)",
        color="#8d5524",
        fontsize=9,
    )

    fig.text(
        0.015,
        0.015,
        (
            "Story for slide: measured energy-loss response -> reconstructed dielectric "
            "function -> plasmon behavior"
        ),
        fontsize=9,
        color="#444444",
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / "EELS_summary_combined.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
