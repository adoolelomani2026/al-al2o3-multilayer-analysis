from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "data" / "raw" / "ellipsometry" / "Aneesa.xlsx"
PROCESSED_DIR = ROOT / "data" / "processed" / "emt"
FIGURE_DIR = ROOT / "outputs" / "emt"

NOMINAL_F_AL = 5.0 / (5.0 + 7.0)
SENSITIVITY_AL_THICKNESSES_NM = np.array([6.471, 7.056, 7.005, 7.004, 7.003])
SENSITIVITY_AL2O3_THICKNESSES_NM = np.array([6.851, 7.058, 7.013, 6.998, 6.999])
SENSITIVITY_F_AL = SENSITIVITY_AL_THICKNESSES_NM.sum() / (
    SENSITIVITY_AL_THICKNESSES_NM.sum() + SENSITIVITY_AL2O3_THICKNESSES_NM.sum()
)


def load_optical_function(sheet_name: str) -> pd.DataFrame:
    data = pd.read_excel(WORKBOOK, sheet_name=sheet_name, header=1)
    data.columns = [str(column).strip() for column in data.columns]
    return data[["nm", "n", "k"]].apply(pd.to_numeric, errors="coerce").dropna()


def calculate_emt(
    eps_al: np.ndarray, eps_al2o3: np.ndarray, filling_fraction: float
) -> tuple[np.ndarray, np.ndarray]:
    eps_parallel = filling_fraction * eps_al + (1.0 - filling_fraction) * eps_al2o3
    eps_perpendicular = 1.0 / (
        filling_fraction / eps_al + (1.0 - filling_fraction) / eps_al2o3
    )
    return eps_parallel, eps_perpendicular


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    al = load_optical_function("Al")
    al2o3 = load_optical_function("Al2O3")
    if not np.allclose(al["nm"], al2o3["nm"]):
        raise ValueError("Al and Al2O3 workbook wavelengths do not match.")

    wavelength = al["nm"].to_numpy()
    eps_al = (al["n"].to_numpy() + 1j * al["k"].to_numpy()) ** 2
    eps_al2o3 = (al2o3["n"].to_numpy() + 1j * al2o3["k"].to_numpy()) ** 2

    cases = {
        "Nominal design": (NOMINAL_F_AL, *calculate_emt(eps_al, eps_al2o3, NOMINAL_F_AL)),
        "Secondary sensitivity": (
            SENSITIVITY_F_AL,
            *calculate_emt(eps_al, eps_al2o3, SENSITIVITY_F_AL),
        ),
    }

    output = pd.DataFrame({"wavelength_nm": wavelength})
    summary_rows = []
    for case_name, (filling_fraction, eps_parallel, eps_perpendicular) in cases.items():
        key = case_name.lower().replace(" ", "_")
        product = eps_parallel.real * eps_perpendicular.real
        output[f"Re_eps_parallel_{key}"] = eps_parallel.real
        output[f"Im_eps_parallel_{key}"] = eps_parallel.imag
        output[f"Re_eps_perpendicular_{key}"] = eps_perpendicular.real
        output[f"Im_eps_perpendicular_{key}"] = eps_perpendicular.imag
        output[f"real_part_product_{key}"] = product
        summary_rows.append(
            {
                "geometry": case_name,
                "f_Al": filling_fraction,
                "Al_total_nm": (
                    25.0 if case_name == "Nominal design" else SENSITIVITY_AL_THICKNESSES_NM.sum()
                ),
                "Al2O3_total_nm": (
                    35.0
                    if case_name == "Nominal design"
                    else SENSITIVITY_AL2O3_THICKNESSES_NM.sum()
                ),
                "Re_eps_parallel_min": eps_parallel.real.min(),
                "Re_eps_parallel_max": eps_parallel.real.max(),
                "Re_eps_perpendicular_min": eps_perpendicular.real.min(),
                "Re_eps_perpendicular_max": eps_perpendicular.real.max(),
                "real_part_product_min": product.min(),
                "real_part_product_max": product.max(),
                "type_II_at_all_wavelengths": bool(
                    np.all((eps_parallel.real < 0) & (eps_perpendicular.real > 0))
                ),
            }
        )

    output.to_csv(PROCESSED_DIR / "emt_tensor_dual_filling_fraction.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(
        PROCESSED_DIR / "emt_tensor_dual_filling_fraction_summary.csv", index=False
    )

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "figure.dpi": 150,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.0), constrained_layout=True)
    colors = {"parallel": "#285f9b", "perpendicular": "#b54332"}
    styles = {"Nominal design": "-", "Secondary sensitivity": "--"}

    for case_name, (filling_fraction, eps_parallel, eps_perpendicular) in cases.items():
        suffix = (
            rf"$f_{{\rm Al}}={filling_fraction:.3f}$"
        )
        axes[0].plot(
            wavelength,
            eps_parallel.real,
            color=colors["parallel"],
            linestyle=styles[case_name],
            linewidth=2.1,
            label=rf"In-plane $\mathrm{{Re}}\,\varepsilon_{{\parallel}}$; {suffix}",
        )
        axes[0].plot(
            wavelength,
            eps_perpendicular.real,
            color=colors["perpendicular"],
            linestyle=styles[case_name],
            linewidth=2.1,
            label=rf"Out-of-plane $\mathrm{{Re}}\,\varepsilon_{{\perp}}$; {suffix}",
        )
        axes[1].plot(
            wavelength,
            eps_parallel.imag,
            color=colors["parallel"],
            linestyle=styles[case_name],
            linewidth=2.1,
            label=rf"In-plane $\mathrm{{Im}}\,\varepsilon_{{\parallel}}$; {suffix}",
        )
        axes[1].plot(
            wavelength,
            eps_perpendicular.imag,
            color=colors["perpendicular"],
            linestyle=styles[case_name],
            linewidth=2.1,
            label=rf"Out-of-plane $\mathrm{{Im}}\,\varepsilon_{{\perp}}$; {suffix}",
        )
        axes[2].plot(
            wavelength,
            eps_parallel.real * eps_perpendicular.real,
            color="#394867",
            linestyle=styles[case_name],
            linewidth=2.2,
            label=rf"{case_name}; {suffix}",
        )

    axes[0].set_title("(a) Real tensor components")
    axes[0].set_ylabel("Real permittivity")
    axes[1].set_title("(b) Loss components")
    axes[1].set_ylabel("Imaginary permittivity")
    axes[2].set_title("(c) Type-II sign test")
    axes[2].set_ylabel(
        r"$\mathrm{Re}\,\varepsilon_{\parallel}\,"
        r"\mathrm{Re}\,\varepsilon_{\perp}$"
    )
    axes[2].axhspan(
        min(
            (
                values[1].real * values[2].real
            ).min()
            for values in cases.values()
        )
        - 5,
        0,
        color="#dce3f0",
        alpha=0.7,
        zorder=0,
    )
    axes[2].text(
        0.94,
        0.93,
        "Type II",
        transform=axes[2].transAxes,
        ha="right",
        va="top",
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "none"},
    )

    for axis in axes:
        axis.axhline(0, color="black", linestyle=":", linewidth=0.9)
        axis.set_xlim(wavelength.min(), wavelength.max())
        axis.set_xlabel("Wavelength (nm)")
        axis.minorticks_on()
        axis.grid(True, which="major", color="#d9d9d9", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        axis.legend(frameon=False, loc="best")

    for suffix in ("png", "pdf"):
        fig.savefig(
            FIGURE_DIR / f"emt_tensor_estimate.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
