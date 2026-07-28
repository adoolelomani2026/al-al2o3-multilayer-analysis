from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "raw" / "eels"
FIG_DIR = ROOT / "outputs" / "eels"


def load_curve(name: str, value_name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / name, header=None, names=["energy_eV", value_name])
    df = df[df["energy_eV"] >= 0].copy()
    return df.reset_index(drop=True)


def crossing_points(df: pd.DataFrame, ycol: str) -> list[float]:
    x = df["energy_eV"].to_numpy()
    y = df[ycol].to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    flips = np.where(np.signbit(y[:-1]) != np.signbit(y[1:]))[0]
    crossings = []
    for i in flips:
        x0, x1 = x[i], x[i + 1]
        y0, y1 = y[i], y[i + 1]
        if y1 == y0:
            continue
        crossings.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    return crossings


def main() -> None:
    eps1 = load_curve("Epsilon 1.csv", "eps1")
    eps2 = load_curve("Epsilon 2.csv", "eps2")
    loss = load_curve("Energy-loss function.csv", "loss")

    xmax = 30
    eps1 = eps1[eps1["energy_eV"] <= xmax].copy()
    eps2 = eps2[eps2["energy_eV"] <= xmax].copy()
    loss = loss[loss["energy_eV"] <= xmax].copy()

    zero_crossings = [x for x in crossing_points(eps1, "eps1") if 10 <= x <= 20]
    first_cross = zero_crossings[0] if zero_crossings else 14.5
    second_cross = zero_crossings[1] if len(zero_crossings) > 1 else 15.25
    plasmon_idx = loss["loss"].idxmax()
    plasmon_energy = float(loss.loc[plasmon_idx, "energy_eV"])
    plasmon_value = float(loss.loc[plasmon_idx, "loss"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 11,
            "axes.titlesize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(6.7, 5.8),
        sharex=True,
        gridspec_kw={"hspace": 0.1},
        constrained_layout=False,
    )

    near_uv = (3.1, 6.2)
    colors = {"eps1": "#0E5A74", "eps2": "#C71F1D", "loss": "#2A9D4B", "annot": "#D08A47"}

    def panel_tag(ax, text: str, x: float, y: float) -> None:
        ax.text(
            x,
            y,
            text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.5,
            fontweight="bold",
            color="#202020",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": (1, 1, 1, 0.86), "edgecolor": "none"},
            zorder=10,
            clip_on=True,
        )

    for ax in axes:
        ax.axvspan(*near_uv, color="#EFD9B5", alpha=0.20, lw=0, zorder=0)
        ax.set_xlim(0, xmax)
        ax.grid(True, color="#D6D6D6", lw=0.8, alpha=0.8)
        ax.grid(True, which="minor", color="#EAEAEA", lw=0.6, alpha=0.95)
        ax.minorticks_on()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].plot(eps1["energy_eV"], eps1["eps1"], color=colors["eps1"], lw=2.2)
    panel_tag(axes[0], "(a)", 0.09, 0.97)
    axes[0].axhline(0, color="#6E6E6E", lw=1.0)
    axes[0].axvline(first_cross, color=colors["annot"], lw=1.6, ls="--", alpha=0.9)
    axes[0].set_ylabel(r"$\varepsilon_1$")
    axes[0].set_ylim(-0.55, 3.3)
    axes[0].annotate(
        r"$\varepsilon_1 = 0$",
        xy=(first_cross, 0),
        xytext=(17.1, 1.55),
        color=colors["annot"],
        fontsize=11.5,
        arrowprops={"arrowstyle": "->", "color": colors["annot"], "lw": 1.25},
    )
    axes[0].text(
        5.6,
        3.02,
        "Near-UV\n(200-400 nm)",
        ha="center",
        va="top",
        fontsize=8.0,
        color="#8B6238",
    )

    axes[1].plot(eps2["energy_eV"], eps2["eps2"], color=colors["eps2"], lw=2.0)
    panel_tag(axes[1], "(b)", 0.02, 0.97)
    axes[1].set_ylabel(r"$\varepsilon_2$")
    axes[1].set_ylim(0, 2.65)

    axes[2].plot(loss["energy_eV"], loss["loss"], color=colors["loss"], lw=2.0)
    panel_tag(axes[2], "(c)", 0.02, 0.97)
    axes[2].scatter([plasmon_energy], [plasmon_value], color=colors["annot"], s=26, zorder=5)
    axes[2].annotate(
        "Bulk-plasmon peak\n(~15 eV)",
        xy=(plasmon_energy, plasmon_value),
        xytext=(16.5, 1.12),
        color=colors["annot"],
        fontsize=8.8,
        ha="left",
        arrowprops={"arrowstyle": "->", "color": colors["annot"], "lw": 1.2},
    )
    axes[2].set_ylabel(r"$\mathrm{Im}[-1/\varepsilon]$")
    axes[2].set_xlabel("Energy loss (eV)")
    axes[2].set_ylim(0, 1.42)

    for ax in axes:
        ax.set_xticks(np.arange(0, 31, 5))

    fig.subplots_adjust(left=0.11, right=0.995, top=0.975, bottom=0.105, hspace=0.12)
    FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(FIG_DIR / "eels_uv_assessment.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_DIR / "eels_uv_assessment.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {FIG_DIR / 'eels_uv_assessment.png'}")
    print(f"Saved {FIG_DIR / 'eels_uv_assessment.pdf'}")
    print(f"Zero crossings near plasmon: {first_cross:.2f} eV, {second_cross:.2f} eV")
    print(f"Loss peak: {plasmon_energy:.2f} eV")


if __name__ == "__main__":
    main()
