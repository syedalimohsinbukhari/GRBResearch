"""Data - fit residual for GRB080916C's original 6-pulse model (norris1-6, before norris2-1/norris2-2 were
split off norris2's overshadowing tail -- see fitter.py's P0 comment for that history). Diagnostic only,
not a paper figure: shows why norris2-1 (~1.5s) and norris2-2 (~2.4s) were added, by plotting exactly the
residual bumps that motivated them.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d

from grb_research import update_style
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import NorrisFitter, norris_pulse

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_080916C = LC_DIR / "GRB080916009"
ENERGY_LOW, ENERGY_HIGH = 10, 400

START1 = -1
END1 = 70

update_style()

dat = [f for f in __import__("os").listdir(f"{GRB_080916C}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = [i.split(".")[0] for i in dat if "n" in i]

nai_data = [lightcurve_data(f"{GRB_080916C}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]
t1, r1, b1 = nai_data[0]

mask_ = np.logical_and(t1 > START1, t1 < END1)
y = (r1 - b1)[mask_]
y_norm = y
Y_MAX_CTS_PER_S = np.max(y)
# y_norm = y / Y_MAX_CTS_PER_S
t_w = t1[mask_]

# The ORIGINAL 6-pulse model (norris1-6), before norris2-1/norris2-2 were added -- this is what
# fitter.py's P0 comment calls "the already-converged 6-pulse fit" used to seed the 8-pulse refit.
P0_6 = [
    (1511, -0.7, 2.34, 0.471),  # norris1: precursor
    (640, 0.3, 0.88, 6),  # norris2: broad envelope/tail
    (400, 5.3, 0.8, 0.3),  # norris3: ~5.7-5.9s peak
    # (0.3 * Y_MAX_CTS_PER_S, 1.3, 43, 14),  # norris4: broad TR3 pedestal
    # (0.1, 20, 9, 1),
    (100, 52, 1, 0.7),  # norris5
    (100, 61, 0.3, 3),  # norris6
    # (0.3, 66, 0.3, 1),
]

nf = NorrisFitter(t_w, y_norm)
nf.fit(p0=P0_6)

pred = sum(norris_pulse(t_w, nf.params[i * 4: (i + 1) * 4]) for i in range(len(P0_6)))
resid = y_norm - pred
# 5-bin moving average -- the same smoothing used when this residual was originally inspected to seed
# norris2-1/norris2-2's p0 (see fitter.py's P0 comment); damps the per-bin Poisson noise enough to show
# the two coherent bumps without smoothing away real structure at this timescale (5 bins ~ 0.32s).
resid_smooth = uniform_filter1d(resid, size=5)

print("=== converged 6-pulse parameters ===")
for i in range(len(P0_6)):
    A, ts, tau1, tau2 = nf.params[i * 4: (i + 1) * 4]
    print(f"  norris{i + 1}: A={A:.4f} ts={ts:9.4f} tau1={tau1:10.4f} tau2={tau2:8.4f}")

# --- Plot: data+fit on top, residual below, zoomed to the early phase where norris2-1/norris2-2 live.
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(12, 8), sharex=True, height_ratios=[2, 1], gridspec_kw={"hspace": 0.08}
)

nf.plot_fit(show_individuals=True, axis=ax_top)
# ax_top.plot(t_w, y_norm, color="tab:blue", lw=1.0, alpha=0.6, label="10-400 keV NaI\nBackground Subtracted")
# ax_top.plot(t_w, pred, color="black", lw=1.8, label="6-pulse total fit (norris1-6)")
# ax_top.set_ylabel("Normalized count rate")
# ax_top.set_title("GRB080916C: data - fit residual, original 6-pulse model")
# ax_top.legend(loc="upper right")

ax_bot.axhline(0, color="gray", lw=0.8, ls="--")
ax_bot.plot(t_w, resid, color="tab:red", lw=0.8, alpha=0.4, label="raw residual (per-bin)")
ax_bot.plot(t_w, resid_smooth, color="tab:red", lw=1.8, label="smoothed residual (5-bin moving average)")
ax_bot.axvspan(1.41, 1.54, color="tab:green", alpha=0.15, label="norris2-1 region (~1.5s)")
ax_bot.axvspan(2.43, 2.62, color="tab:purple", alpha=0.15, label="norris2-2 region (~2.4s)")
ax_bot.set_xlabel("Time since trigger [s]")
ax_bot.set_ylabel("Residual\n(data - fit)")
ax_bot.legend(loc="upper right", fontsize=8)
ax_bot.set_xlim(-1, 70)

fig_path = Path(__file__).parent / "fitter_EXPERIMENT_residual_GRB080916009"
save_fig(fig, fig_path)
print(f"saved {fig_path}.png/.pdf")
