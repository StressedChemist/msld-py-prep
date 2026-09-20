#!/usr/bin/env python3
"""Plot lambda and signed-volume traces from either solvated test runner."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", default="lambda_and_signed_volume.png")
    args = parser.parse_args()

    with (args.run / "lambda_stereo.csv").open() as handle:
        rows = list(csv.DictReader(handle))

    time = [float(row["time_ps"]) for row in rows]
    lambdas = [float(row["lambda_meso"]) for row in rows]
    meso = [float(row["meso_signed_volume_nm3"]) for row in rows]
    threo = [float(row["threo_signed_volume_nm3"]) for row in rows]

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(10, 6),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1.3]},
        constrained_layout=True,
    )
    axes[0].plot(time, lambdas, color="#3b6fb6", linewidth=0.65)
    axes[0].set_ylabel(r"$\lambda_{meso}$")
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].set_yticks([0.0, 0.5, 1.0])
    axes[0].set_title("Solvated meso-erythritol/threitol lambda test")

    axes[1].plot(time, meso, label="meso C3", color="#c44e52", linewidth=0.75)
    axes[1].plot(time, threo, label="threo C3", color="#55a868", linewidth=0.75)
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("Simulation time (ps)")
    axes[1].set_ylabel(r"C3 signed volume (nm$^3$)")
    axes[1].legend(frameon=True, framealpha=0.9, edgecolor="none", ncol=2)

    figure.savefig(args.run / args.output, dpi=180)


if __name__ == "__main__":
    main()
