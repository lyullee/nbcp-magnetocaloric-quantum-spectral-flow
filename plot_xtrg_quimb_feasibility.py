"""Plot the N=9 purified-MPO feasibility boundary."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
data = json.loads(
    (ROOT / "results" / "xtrg_quimb_n9_purified_stabilized_convergence.json")
    .read_text(encoding="utf-8")
)
runs = data["runs"]
d = np.array([r["max_bond"] for r in runs])
ds = np.array([r["final"]["errors"]["entropy_kB_per_spin"] for r in runs])
de = np.array([r["final"]["errors"]["energy_K"] for r in runs])
dz = np.array([r["final"]["errors"]["logZ"] for r in runs])
wall = np.array([r["wall_s"] for r in runs])
mineig = np.array([r["diagnostic"]["final_min_density_eigenvalue"] for r in runs])

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))

axes[0].axhline(0, color="0.35", lw=0.8)
axes[0].axhspan(-0.001, 0.001, color="#54a24b", alpha=0.13, label="±0.001 target")
axes[0].plot(d, ds, "o-", lw=1.6, color="#4c78a8")
for x, y in zip(d, ds):
    axes[0].annotate(f"{y:+.4f}", (x, y), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)
axes[0].set(xlabel="Maximum bond dimension D", ylabel="Entropy error (kB/spin)", title="90 mK entropy error")
axes[0].legend(frameon=False, fontsize=8)

axes[1].plot(d, np.abs(de), "o-", label="|energy error| (K)", color="#f58518")
axes[1].plot(d, np.abs(dz), "s--", label="|logZ error|", color="#e45756")
axes[1].plot(d, np.abs(ds), "^:", label="|entropy error|", color="#4c78a8")
axes[1].set_yscale("log")
axes[1].set(xlabel="Maximum bond dimension D", ylabel="Absolute error", title="Errors are non-monotonic")
axes[1].legend(frameon=False, fontsize=8)

axes[2].plot(d, wall, "o-", color="#72b7b2")
for x, y in zip(d, wall):
    axes[2].annotate(f"{y:.0f} s", (x, y), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)
axes[2].set(xlabel="Maximum bond dimension D", ylabel="Wall time (s)", title="Local N=9 cost")
axes[2].text(0.03, 0.96, f"min rho eigenvalue: {mineig.min():.1e}", transform=axes[2].transAxes, va="top", fontsize=8)

fig.suptitle("Purified MPO feasibility test: M0, N=9, B=1.50 T", y=1.03, fontsize=12)
fig.tight_layout()
out = ROOT / "figures" / "xtrg_quimb_n9_feasibility.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print(out)
