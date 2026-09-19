#!/usr/bin/env python3
"""##CHRIS 2026-09-17: the two Level 2 panels. Analysis only; numbers from paper2_level2_20260917.py.

Left  : excess work over the quasi-static value against piston speed, with the fitted A u^2 and, for
        contrast, the line the zero-lag (Enskog) friction would have predicted.
Right : the running Green-Kubo integral, showing it cancel at the sound traversal time L/c_s.
"""
import glob, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper2_level2_20260917 import work_cells, zeta_running, WQS, DWQS, DX, NEW
import tests_20260913 as T

OUT = os.path.join(T.PLOTS, "260917_level2_work_and_friction")
L_RIGHT, ETA_RIGHT, CS_RIGHT = 35.320, 0.1112, 1.788   # compressed compartment during the hold

cells = work_cells()
u = np.array(list(cells)); y = np.array([c[0] for c in cells.values()]) - WQS
e = np.array([c[1] for c in cells.values()]); w = 1.0 / e ** 2
X = np.vstack([np.ones_like(u), u ** 2]).T
C = np.linalg.inv(X.T @ (X * w[:, None])); a, A = C @ (X.T @ (w * y)); da, dA = np.sqrt(np.diag(C))

lags, Cf, run, self_term = zeta_running(sorted(glob.glob(f"{NEW}/u0.03/ev_*.csv"))
                                        + sorted(glob.glob(f"{NEW}/u0.05/ev_*.csv")))

fig, (ax, axz) = plt.subplots(1, 2, figsize=(11.0, 4.4))

ug = np.linspace(0, 0.21, 200)
ax.axhline(0.0, color="0.75", lw=0.8)
ax.plot(ug, self_term * DX * ug, color="tab:red", lw=1.4, ls="--",
        label=r"$\zeta_{t\to0}\,\Delta x\,u$ (zero-lag friction, %.1f$u$)" % (self_term * DX))
ax.plot(ug, a + A * ug ** 2, color="tab:blue", lw=1.6,
        label=r"fit $W(0)+Au^2$, $A = %.1f \pm %.1f$" % (A, dA))
ax.plot(ug, 25.0 * ug ** 2, color="tab:green", lw=1.2, ls=":", label=r"$N_s m u^2/2$ (no free parameter)")
ax.errorbar(u, y, yerr=e, fmt="o", ms=4.0, color="k", elinewidth=1.0, capsize=2.5, zorder=5,
            label="measured (60–100 seeds each)")
ax.set_xlabel(r"piston speed $u$  [$\sigma/\tau$]")
ax.set_ylabel(r"$\langle W\rangle - W_{\rm qs}^{\rm finite}$  [$k_BT$]")
ax.set_xlim(0, 0.21); ax.set_ylim(-0.25, 1.3)
ax.set_title("Excess work is second order in the piston speed", fontsize=10)
ax.legend(fontsize=8, loc="upper left", framealpha=0.95)

axz.axhline(0.0, color="0.75", lw=0.8)
axz.plot(lags, run, color="tab:purple", lw=1.6)
axz.axvline(L_RIGHT / CS_RIGHT, color="tab:orange", lw=1.2, ls="--")
axz.annotate(r"$L/c_s = %.1f\,\sigma$" % (L_RIGHT / CS_RIGHT), xy=(L_RIGHT / CS_RIGHT, 1.6),
             xytext=(L_RIGHT / CS_RIGHT + 2.5, 2.0), fontsize=9, color="tab:orange")
axz.annotate(r"zero-lag (Enskog) value %.2f" % self_term, xy=(0.6, self_term),
             xytext=(6.0, 2.35), fontsize=9, color="0.3",
             arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
axz.set_xlabel(r"upper cut-off $t_{\rm cut}$  [$\sigma$]")
axz.set_ylabel(r"$\beta\int_0^{t_{\rm cut}}\langle\delta F(0)\delta F(t)\rangle\,dt$")
axz.set_xlim(0, 60); axz.set_ylim(-1.1, 2.9)
axz.set_title(r"The friction integral cancels at the sound traversal time", fontsize=10)
k = (lags >= 20) & (lags <= 60)
axz.hlines(run[k].mean(), 20, 60, color="0.35", lw=1.2, ls="-.")
axz.text(41, run[k].mean() + 0.10, r"mean over one recurrence $2L/c_s$: %+.2f" % run[k].mean(),
         fontsize=8, color="0.35", ha="right")
axz.text(0.97, 0.06, "190 holds of 180 $\\sigma$, $\\eta = %.3f$" % ETA_RIGHT, transform=axz.transAxes,
         ha="right", fontsize=8, color="0.35")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print(f"written {OUT}.png / .pdf   A = {A:.2f} ± {dA:.2f}, W(0) - W_qs = {a:+.4f} ± {da:.4f}")
