"""Created on Sep 20 12:48:47 2025"""

from pathlib import Path

from matplotlib import pyplot as plt

from .grb_atomic import Parameter, ParameterSet, CovarianceMatrix
from .grb_calculations import FluxFluenceCalculator, ModelResampler, mc_e_iso_sampler
from .grb_constants import (
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
    TICK_FONT_SIZE,
    TITLE_FONT_SIZE,
    ANNOTATION_FONT_SIZE,
    FIGURE_DPI,
    SAVE_DPI,
    DEFAULT_FIGURE_SIZE,
    LINE_WIDTH,
    MARKER_SIZE,
    AXES_LINE_WIDTH,
    TICK_MAJOR_SIZE,
    TICK_MAJOR_WIDTH,
    TICK_MINOR_SIZE,
    TICK_MINOR_WIDTH,
    TICK_DIRECTION,
    GRID_ALPHA,
    GRID_LINESTYLE,
    MODEL_ORDER,
    long_to_short,
    short_to_long,
)
from .grb_core import GRB, GRBCatalog, prepare_grbs
from .grb_cross_validation import DEFAULT_BB_PARAM_NAMES, SIMPLE_TO_BB
from .grb_model import Model, ModelSet
from .grb_sed import SpectralModels
from .grb_stability import ModelComparison
from .grb_styles import GRBPlotStyle
from .grb_time import EpisodeTypes, TimeInterval, TimeIntervalSet
from .grb_utils import (
    EpisodeMarkerResolver,
    break_e_to_e_peak,
    plot_per_episode,
    save_value_error_as_parquet,
)
from .safe_good_best import pick_best_single_model


def find_project_root(marker="results.json"):
    """Find the root of the project."""
    path = Path(__file__).resolve()
    for parent in path.parents:
        if (parent / marker).exists():
            return parent
    raise RuntimeError("Project root not found")


def update_style():
    """Update style for publication-ready figures.

    All values are driven by constants in grb_constants — a single source of truth.
    Changing a constant there automatically updates every figure.
    """
    plt.rcParams.update(
        {
            # Font
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "font.size": LABEL_FONT_SIZE,
            "axes.labelsize": LABEL_FONT_SIZE,
            "axes.titlesize": TITLE_FONT_SIZE,
            "xtick.labelsize": TICK_FONT_SIZE,
            "ytick.labelsize": TICK_FONT_SIZE,
            "legend.fontsize": LEGEND_FONT_SIZE,
            "legend.title_fontsize": LEGEND_TITLE_FONT_SIZE,
            # Figure
            "figure.dpi": FIGURE_DPI,
            "figure.figsize": DEFAULT_FIGURE_SIZE,
            "savefig.dpi": SAVE_DPI,
            "savefig.bbox": "tight",
            # Lines / markers
            "lines.linewidth": LINE_WIDTH,
            "lines.markersize": MARKER_SIZE,
            "axes.linewidth": AXES_LINE_WIDTH,
            # 'scatter.size': MARKER_SIZE,
            # Ticks
            "xtick.major.size": TICK_MAJOR_SIZE,
            "xtick.major.width": TICK_MAJOR_WIDTH,
            "xtick.minor.size": TICK_MINOR_SIZE,
            "xtick.minor.width": TICK_MINOR_WIDTH,
            "xtick.minor.visible": True,
            "xtick.direction": TICK_DIRECTION,
            "ytick.major.size": TICK_MAJOR_SIZE,
            "ytick.major.width": TICK_MAJOR_WIDTH,
            "ytick.minor.size": TICK_MINOR_SIZE,
            "ytick.minor.width": TICK_MINOR_WIDTH,
            "ytick.minor.visible": True,
            "ytick.direction": TICK_DIRECTION,
            # Grid
            "axes.grid": True,
            "grid.alpha": GRID_ALPHA,
            "grid.linestyle": GRID_LINESTYLE,
            "axes.axisbelow": True,
        }
    )
