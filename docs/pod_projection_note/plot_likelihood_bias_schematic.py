#!/usr/bin/env python3
"""Draw a schematic of a POD likelihood bias versus its chi2 offset."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    z = np.linspace(-2.0, 2.0, 401)
    direct = z**2
    delta = -1.32 + 0.08 * (z - 1.0) ** 2
    pod = direct + delta
    z_pod = z[np.argmin(pod)]
    delta_zero = np.interp(0.0, z, delta)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    ax = axes[0]
    ax.plot(z, direct, lw=2.5, color="#1f77b4", label=r"direct $\chi^2(z)$")
    ax.plot(z, pod, lw=2.5, color="#7a3e9d", label=r"POD $\chi^2(z)$")
    ax.axvline(0.0, color="0.45", lw=1, ls=":")
    ax.axvline(z_pod, color="#7a3e9d", lw=1.2, ls="--")
    ax.plot(0.0, 0.0, "o", color="#1f77b4")
    ax.plot(z_pod, pod.min(), "o", color="#7a3e9d")
    ax.annotate(
        fr"direct minimum: $z=0$",
        xy=(0.0, 0.0),
        xytext=(-1.85, 2.8),
        arrowprops={"arrowstyle": "->", "color": "#1f77b4"},
        color="#1f77b4",
    )
    ax.annotate(
        fr"POD minimum: $z\simeq {z_pod:.2f}$",
        xy=(z_pod, pod.min()),
        xytext=(0.45, -0.5),
        arrowprops={"arrowstyle": "->", "color": "#7a3e9d"},
        color="#7a3e9d",
    )
    ax.set(
        xlabel=r"scan coordinate $z$ [direct-fit $\sigma$]",
        ylabel=r"raw $\chi^2$",
        title="Likelihoods that are actually minimised",
        xlim=(-2, 2),
    )
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="upper right")

    ax = axes[1]
    ax.plot(
        z,
        delta,
        lw=2.5,
        color="#c65d28",
        label=r"$\Delta\chi^2=\chi^2_{\rm POD}-\chi^2_{\rm direct}$",
    )
    ax.axvline(0.0, color="0.45", lw=1, ls=":")
    ax.axvline(1.0, color="#c65d28", lw=1.2, ls="--")
    ax.plot(1.0, delta.min(), "o", color="#c65d28")
    ax.annotate(
        r"minimum of $\Delta\chi^2$ is away from zero",
        xy=(1.0, delta.min()),
        xytext=(-1.9, -0.95),
        arrowprops={"arrowstyle": "->", "color": "#c65d28"},
        color="#c65d28",
    )
    ax.annotate(
        fr"$\Delta\chi^2(0)={delta_zero:.2f}$ is mostly an offset",
        xy=(0.0, delta_zero),
        xytext=(-1.9, -0.55),
        arrowprops={"arrowstyle": "->", "color": "#c65d28"},
        color="#c65d28",
    )
    ax.text(
        0.04,
        0.95,
        r"The small local slope is divided by"
        "\n"
        r"the direct-fit curvature when the"
        "\n"
        r"full POD curve is minimised.",
        transform=ax.transAxes,
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.7"},
    )
    ax.set(
        xlabel=r"scan coordinate $z$ [direct-fit $\sigma$]",
        ylabel=r"$\Delta\chi^2$",
        title="Difference curve shown in the lower pad",
        xlim=(-2, 2),
    )
    ax.grid(alpha=0.25)

    output = Path(__file__).with_name("likelihood_bias_schematic.pdf")
    fig.savefig(output)
    print(output)


if __name__ == "__main__":
    main()
