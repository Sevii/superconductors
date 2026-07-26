#!/usr/bin/env python3
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_bridge_covariance import projected_bridge_matrix, middle_jacobi, asymptotic_middle_stiffness

# Fixed-kernel versus covariant-kernel projected bridge at a representative k.
kx, ky = 0.73, 1.11
avalues = np.linspace(-0.08, 0.08, 401)
fixed = []
covariant = []
for ay in avalues:
    fixed.append(abs(projected_bridge_matrix(kx, ky, ay, "X", covariant_kernel=False)[0,1]))
    covariant.append(abs(projected_bridge_matrix(kx, ky, ay, "X", covariant_kernel=True)[0,1]))
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(avalues, fixed, label="fixed bridge kernel (superseded)")
ax.plot(avalues, covariant, label="Peierls-covariant kernel")
ax.set_xlabel(r"uniform twist $A_y$")
ax.set_ylabel(r"$|[P(A)B(A)P(A)]_{12}|$")
ax.set_title("Bridge covariance audit")
ax.grid(True, alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "bridge_covariance_audit.pdf", bbox_inches="tight")
fig.savefig(OUT / "bridge_covariance_audit.png", dpi=220, bbox_inches="tight")
plt.close(fig)

# Middle-filling finite-size convergence for the verifier parameters.
U, gx, gy = 1.0, 0.8, 0.3
gx2, gy2 = gx*gx, gy*gy
volumes = np.arange(4, 121, 2)
stiffness = []
for V in volumes:
    q = middle_jacobi(int(V), U, gx2, gy2)
    stiffness.append(np.linalg.eigvalsh(q)[0] / (4.0 * V))
limit = asymptotic_middle_stiffness(U, gx2, gy2)
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(volumes, stiffness, label=r"$\lambda_{\min}(Q_y^{(V)})/(4V)$")
ax.axhline(limit, linestyle="--", label=r"endpoint limit $J/8$")
ax.set_xlabel("block capacity V")
ax.set_ylabel("middle-filling stiffness")
ax.set_title(r"Endpoint localization for $U=1$, $g_X=0.8$, $g_Y=0.3$")
ax.grid(True, alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "middle_stiffness_convergence.pdf", bbox_inches="tight")
fig.savefig(OUT / "middle_stiffness_convergence.png", dpi=220, bbox_inches="tight")
plt.close(fig)
