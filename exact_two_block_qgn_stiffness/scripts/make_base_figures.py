#!/usr/bin/env python3
"""Regenerate the two analytic base-model figures used by the draft."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

q = np.linspace(-np.pi, np.pi, 800)
params = [
    (0.5, 2.0, r"Block 1: $\xi_1=1/2$, $K_{11}=2$"),
    (1.0, 5.0, r"Block 2: $\xi_2=1$, $K_{22}=5$"),
]
fig, ax = plt.subplots(figsize=(6.4, 3.7))
for xi, kaa, label in params:
    form = j0(2 * xi * np.sin(q / 2))
    energy = kaa * (1 - np.abs(form)) / 4
    ax.plot(q, energy, label=label)
ax.set_xlabel(r"$Q_x$ at $Q_y=0$")
ax.set_ylabel(r"$E_{\mathrm{pair},a}(Q_x)$ for $U=1$")
ax.set_xlim(-np.pi, np.pi)
ax.grid(True, alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "pair_branches.pdf", bbox_inches="tight")
plt.close(fig)

p = np.linspace(0, 2, 1201)
d = np.abs(p - 1)
stiffness = d * (1 - d) / 16
fig, ax = plt.subplots(figsize=(6.4, 3.7))
ax.plot(p, stiffness, label=r"$D_s^{\mathrm{gs}}(p)$")
ax.axvline(1, linewidth=0.9, linestyle="--")
ax.set_xlabel(r"total pair density $p=n/V$")
ax.set_ylabel(r"$D_s^{\mathrm{gs}}/U$")
ax.set_xlim(0, 2)
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "unresolved_response_gs.pdf", bbox_inches="tight")
plt.close(fig)
