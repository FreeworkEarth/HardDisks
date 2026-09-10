import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import chi2 as C2

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a1a", "#6b6b6b", "#e2e2e2"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})

# ---------- Fig 1: chi2 distributions ----------
fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
ax = axes[0]
xx = np.linspace(0.02, 12, 600)
styles = {1: "-", 2: "--", 3: "-.", 5: ":"}
for k, ls in styles.items():
    ax.plot(xx, C2.pdf(xx, k), ls=ls, color=INK, lw=1.6, label=f"dof = {k}   (mean = {k})")
ax.set_ylim(0, 0.6); ax.set_xlim(0, 12)
ax.set_xlabel("χ² value"); ax.set_ylabel("probability density")
ax.set_title("What χ² looks like when the hypothesis is TRUE", fontsize=10.5, loc="left")
ax.legend(frameon=False, fontsize=8.5); ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)

ax = axes[1]
k = 1
ax.plot(xx, C2.pdf(xx, k), "-", color=INK, lw=1.6, label="dof = 1")
m = xx >= 7.62
ax.fill_between(xx[m], 0, C2.pdf(xx[m], k), color=ORANGE, alpha=0.8, lw=0, label="p = 0.006  (η = 0.67, χ² = 7.62)")
m2 = xx >= 3.84
ax.fill_between(xx[m2], 0, C2.pdf(xx[m2], k), color=ORANGE, alpha=0.25, lw=0, label="p = 0.05 tail starts at χ² = 3.84")
ax.axvline(0.07, color=BLUE, lw=2, label="η = 0.65: χ² = 0.07 (p = 0.79)")
ax.axvline(7.62, color=ORANGE, lw=1.2, ls="--")
ax.set_ylim(0, 0.15); ax.set_xlim(0, 12)
ax.set_xlabel("χ² value"); ax.set_title("Our two cases, dof = 1 (y-axis zoomed to the tail)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5); ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.text(0.99, 0.01, "p-value = area of the tail to the right of the observed χ²", ha="right", fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("fig_chi2_dist.pdf")

# ---------- Fig 2: the two extrapolations ----------
cases = {
    0.65: (np.array([400, 900, 1600.]), np.array([8.4710, 8.4373, 8.4221]), np.array([0.0043, 0.0031, 0.0032]), 8.4080),
    0.67: (np.array([400, 900, 1600.]), np.array([9.2867, 9.2246, 9.2340]), np.array([0.0081, 0.0087, 0.0053]), 9.3625),
}
fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]})
for j, (eta, (N, Z, e, kr)) in enumerate(cases.items()):
    x = 1 / np.sqrt(N); w = 1 / e**2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ np.diag(w) @ A); p = cov @ A.T @ np.diag(w) @ Z
    r = (Z - A @ p) / e; chi = (r**2).sum(); pv = 1 - C2.cdf(chi, 1)
    ax = axes[0, j]
    xl = np.linspace(0, 0.06, 50)
    ax.plot(xl, p[0] + p[1] * xl, "-", color=INK, lw=1.4, label=f"weighted fit  Z = Z∞ + a/√N")
    ax.errorbar(x, Z, yerr=e, fmt="o", color=BLUE, ms=7, capsize=4, lw=1.4, zorder=3, label="measured Z (N = 400, 900, 1600)")
    ax.errorbar([0], [p[0]], yerr=[np.sqrt(cov[0, 0])], fmt="s", color=ORANGE, ms=8, capsize=4, zorder=4, label=f"Z∞ = {p[0]:.4f} ± {np.sqrt(cov[0,0]):.4f}")
    ax.axhline(kr, color=MUTED, ls="--", lw=1.2, label=f"Kolafa–Rottner  {kr:.4f}")
    for xi, Ni, off in zip(x, N, [(6, -14), (6, 8), (-30, -14)]): ax.annotate(f"N = {int(Ni)}", (xi, Z[list(N).index(Ni)]), textcoords="offset points", xytext=off, fontsize=8, color=MUTED)
    ax.set_title(f"η = {eta}:  χ² = {chi:.2f}, dof = 1, p = {pv:.3f}\n→ 1/√N line {'ACCEPTED' if pv > 0.046 else 'REJECTED'}", fontsize=10, loc="left")
    ax.set_ylabel("Z = PA / (N k_BT)"); ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    axr = axes[1, j]
    axr.axhspan(-1, 1, color=GRID, alpha=0.9, lw=0); axr.axhspan(-2, 2, color=GRID, alpha=0.45, lw=0)
    axr.axhline(0, color=INK, lw=1)
    axr.bar(x, r, width=0.003, color=[ORANGE if abs(v) > 2 else BLUE for v in r])
    for xi, v in zip(x, r): axr.annotate(f"{v:+.2f}", (xi, v), textcoords="offset points", xytext=(0, 5 if v > 0 else -12), ha="center", fontsize=8.5)
    axr.set_ylim(-3.2, 3.2); axr.set_ylabel("(Z − fit) / error"); axr.set_xlabel("1/√N       (0 = infinite box)")
    axr.set_xlim(-0.004, 0.06); axr.grid(True, color=GRID, lw=0.7, axis="x"); axr.set_axisbelow(True)
    for s in ("top", "right"): axr.spines[s].set_visible(False)
    axr.text(0.058, 2.5, f"χ² = Σ(bars)² = {chi:.2f}", ha="right", fontsize=8.5, color=INK)
fig.text(0.99, 0.005, "data: 260908_pressure_final_analysis.md table B; error = max(block SEM, seed SEM); bands = ±1 and ±2 error bars", ha="right", fontsize=7.5, color=MUTED)
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig("fig_extrapolation.pdf")

# ---------- Fig 3: wall layer sketch ----------
fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
eta, r = 0.65, 0.5
for ax, n in zip(axes, [100, 400, 1600]):
    A = n * np.pi * r * r / eta; L = np.sqrt(A)
    ax.add_patch(Rectangle((0, 0), L, L, fc="white", ec=INK, lw=1.4))
    ax.add_patch(Rectangle((0, 0), L, L, fc=ORANGE, alpha=0.35, lw=0))
    ax.add_patch(Rectangle((1, 1), L - 2, L - 2, fc="white", lw=0))
    ax.set_xlim(-2, 46); ax.set_ylim(-2, 46); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"N = {n},  L = {L:.0f} σ\nwall layer (1 σ) = {100*4*L/A:.0f} % of area,  1/√N = {1/np.sqrt(n):.3f}", fontsize=9.5)
fig.suptitle("Same η = 0.65, bigger boxes: the wall layer (orange) shrinks like perimeter/area ∝ 1/√N", fontsize=10.5)
fig.tight_layout(); fig.savefig("fig_wall_layer.pdf")
print("figs written")
