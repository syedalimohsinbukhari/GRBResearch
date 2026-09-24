"""GRB231129C, UNNORMALIZED method, fit over the light curve's own full x.min()/x.max() range.

Same raw-counts/explicit-x_scale scipy.optimize.least_squares path as
experiments/normalized_vs_unnormalized_fit/fit_comparison.py's Method 2 and
GRB080916C/fitter_unnormalized.py -- see that file's docstring for the full rationale (why x_scale='jac'
was rejected, and how the covariance is derived from the least_squares Jacobian the same way
scipy.optimize.curve_fit does internally). Same 5-pulse P0, same dominant-flux LAT-photon assignment and
episode-matching logic as fitter_normalized.py in this folder -- only the fitting method differs.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from _common import (  # noqa: E402 -- sets up sys.path, must import before norris_fit
    T_FULL, Y_RAW_FULL, Y_MAX_CTS_PER_S, FIT_WINDOW, DAT_NAI, GRB_231129, GRB_PAPER_NAME,
    N_PULSES, P0, MAX_NFEV, X_SCALE_TS_TAU1, X_SCALE_TAU2, T05, T95, PLOT_PAD_S, FitResult,
    assign_pulses, photon_summary_for, build_results_df, add_lat_photon_overlay, save_fig,
)
from norris_fit import norris_pulse  # noqa: E402
from grb_research import LINE_WIDTH  # noqa: E402


def total_model(t, params, n_pulses):
    return sum(norris_pulse(t, params[i * 4:(i + 1) * 4]) for i in range(n_pulses))


p0_raw = np.array([(a * Y_MAX_CTS_PER_S, ts, tau1, tau2) for (a, ts, tau1, tau2) in P0]).flatten()
lb = np.tile([0.0, T_FULL.min(), 1e-4, 1e-4], N_PULSES)
ub = np.tile([np.inf, T_FULL.max(), np.inf, np.inf], N_PULSES)
x_scale = np.tile([Y_MAX_CTS_PER_S, X_SCALE_TS_TAU1, X_SCALE_TS_TAU1, X_SCALE_TAU2], N_PULSES)


def residuals(params):
    return total_model(T_FULL, params, N_PULSES) - Y_RAW_FULL


res = least_squares(residuals, p0_raw, bounds=(lb, ub), x_scale=x_scale, max_nfev=MAX_NFEV)
if not res.success:
    raise RuntimeError(f"unnormalized fit did not converge: {res.message}")

# Covariance from the Jacobian at the solution -- see GRB080916C/fitter_unnormalized.py's docstring for
# the formula (matches scipy.optimize.curve_fit's own internal derivation, absolute_sigma=False).
dof = len(res.fun) - len(res.x)
s_sq = (2 * res.cost / dof) if dof > 0 else np.inf
try:
    covariance = np.linalg.inv(res.jac.T @ res.jac) * s_sq
except np.linalg.LinAlgError:
    covariance = np.full((len(res.x), len(res.x)), np.inf)

fit_result = FitResult(params=res.x, covariance=covariance)

photon_pulse = assign_pulses(fit_result.params, N_PULSES)
photon_summary = photon_summary_for(photon_pulse, N_PULSES)

results_df = build_results_df(
    method_label="unnormalized_xscale", fit_result=fit_result, n_pulses=N_PULSES,
    y_max_cts_per_s=Y_MAX_CTS_PER_S, amplitude_is_physical=True, photon_summary=photon_summary,
)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_231129.name}_unnormalized.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")
print(f"fit window: {FIT_WINDOW} (full x.min()/x.max())")

# --- Plot: data + total fit + individual pulses, all in physical counts/s -- built by hand (no
# NorrisFitter/BaseFitter.plot_fit() here), same convention as GRB080916C/fitter_unnormalized.py.
fig, ax = plt.subplots(figsize=(13, 6.5))
ax.plot(T_FULL, Y_RAW_FULL, color="0.6", lw=LINE_WIDTH * 0.6,
        label=f"10-400 keV NaI ({'+'.join(DAT_NAI)}, summed)\nBackground Subtracted")
model_total = total_model(T_FULL, fit_result.params, N_PULSES)
ax.plot(T_FULL, model_total, color="tab:red", lw=LINE_WIDTH, label="Unnormalized fit (explicit x_scale), total")

shades = plt.cm.Oranges(np.linspace(0.4, 0.85, N_PULSES))
for i in range(N_PULSES):
    par = fit_result.params[i * 4:(i + 1) * 4]
    A, ts, tau1, tau2 = par
    ax.plot(T_FULL, norris_pulse(T_FULL, par), ls="--", lw=LINE_WIDTH * 0.7, color=shades[i],
            label=f"Pulse {i + 1} (A={A:.1f}, t_s={ts:.2f}, tau1={tau1:.2f}, tau2={tau2:.2f})")

ax.set_xlabel("Time since trigger [s]")
ax.set_ylabel("Count rate [counts/s]")
ax.set_title(f"{GRB_PAPER_NAME}: unnormalized (x_scale) fit, full range {FIT_WINDOW[0]:.1f}-{FIT_WINDOW[1]:.1f}s")
add_lat_photon_overlay(ax, y_top_data=Y_RAW_FULL.max())
legend = ax.legend(fontsize="x-small", loc="upper left", bbox_to_anchor=(1.08, 1.0))

# Plot-axis bounding only (T05-{PLOT_PAD_S}s .. T95+{PLOT_PAD_S}s) -- does not change what was fitted.
ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_231129.name}_unnormalized"
save_fig(fig, fig_path, bbox_extra_artists=(legend,))
print(f"wrote {fig_path}.png / .pdf")
