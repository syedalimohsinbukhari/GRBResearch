"""Created on Mar 26 10:09:59 2026

grb_plot_style.py
=================
Consistent, publication-ready plot styling for GRB spectral analysis.

Style hierarchy
---------------
  GRB → color (Okabe-Ito colorblind-safe palette)
  Temporal episode → marker shape
  Base model → line style and hatch (for filled regions / bars)
  Model extension → marker fill style and edge weight

Usage
-----
  from grb_plot_style import S # S is a pre-built singleton

  # scatter / line plot
  ax.plot(x, y, **S.line_kwargs("GRB110721A", "CPL"))

  # error bars
  ax.errorbar(x, y, xerr=xe, yerr=ye,
              **S.errorbar_kwargs("GRB080916C", episode=2, extension="BASE+BB"))

  # confidence band
  ax.fill_between(x, lo, hi, **S.band_kwargs("GRB150210A"))

  # legend sections
  ax.legend(handles=S.grb_handles() + S.episode_handles() + S.model_handles())

  # quick-look reference table
  S.summary()
"""

from __future__ import annotations

import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Core style class
# ---------------------------------------------------------------------------


class GRBPlotStyle:
    """
    Central style registry for a GRB spectral-analysis figure set.

    Dimensions
    ----------
    GRBs       : GRB080916C | GRB110721A | GRB110731A | GRB150210A
    Episodes: 1 (I) | 2 (II) | 3 (III)
    Base models: PL | CPL | Band | SBPL
    Extensions : BASE | BASE+BB | BASE+PL+BB
    """

    # ------------------------------------------------------------------
    # GRB → color (Okabe–Ito colorblind-safe, high contrast on white)
    # ------------------------------------------------------------------
    # Okabe-Ito colourblind-safe palette. The paper's four-GRB sample
    # (080916C, 131014A, 140206B, 231129C) is assigned the four most mutually
    # distinct hues; the remaining entries belong to earlier, disjoint samples
    # and deliberately reuse those hues, since no figure mixes the two sets.
    GRB_COLORS: dict[str, str] = {
        # current paper sample
        "GRB080916C": "#0072B2",  # sky blue
        "GRB131014A": "#E69F00",  # orange
        "GRB140206B": "#009E73",  # bluish green
        "GRB231129C": "#CC79A7",  # reddish purple
        # earlier samples, not part of this paper
        "GRB110721A": "#D55E00",  # vermillion
        "GRB110731A": "#009E73",  # bluish green
        "GRB150210A": "#CC79A7",  # reddish purple
        "GRB190114C": "#808000",
    }

    GRB_COLORS_ITERABLE = list(GRB_COLORS.values())

    GRB_SHORT: dict[str, str] = {
        "GRB080916C": "080916C",
        "GRB131014A": "131014A",
        "GRB140206B": "140206B",
        "GRB231129C": "231129C",
        "GRB110721A": "110721A",
        "GRB110731A": "110731A",
        "GRB150210A": "150210A",
        "GRB190114C": "190114C",
    }

    # ------------------------------------------------------------------
    # Episode → marker shape
    # ------------------------------------------------------------------
    EPISODE_MARKERS: dict[int, str] = {
        1: "o",  # circle      — Episode I
        2: "s",  # square      — Episode II
        3: "^",  # triangle up — Episode III
    }

    EPISODE_LABELS: dict[int, str] = {1: "Episode I", 2: "Episode II", 3: "Episode III"}

    # ------------------------------------------------------------------
    # Base model → line style (both simple and precise tupled forms)
    # ------------------------------------------------------------------
    MODEL_LINESTYLES: dict[str, tuple] = {
        "PL": (0, ()),  # solid
        "CPL": (0, (6, 2)),  # dashed
        "Band": (0, (5, 1.5, 1, 1.5)),  # dash-dot
        "SBPL": (0, (1.5, 1.5)),  # dotted
    }

    # Simpler string aliases — useful for colorbars, legends, quick tests
    MODEL_LINESTYLES_SIMPLE: dict[str, str] = {"PL": "-", "CPL": "--", "Band": "-.", "SBPL": ":"}

    # Hatch patterns for bar charts, residual boxes, posterior bands
    MODEL_HATCHES: dict[str, str | None] = {"PL": None, "CPL": "//", "Band": "xx", "SBPL": ".."}

    # ------------------------------------------------------------------
    # Extension → marker fill style + edge weight
    #   "full" — solid filled (BASE: default model, the highest confidence)
    #   "none" — open / hollow (BASE+BB: one extra component)
    #   "bottom" — half-filled (BASE+PL+BB: two extra components)
    # ------------------------------------------------------------------
    EXTENSION_FILLSTYLE: dict[str, str] = {"BASE": "full", "BASE+BB": "none", "BASE+PL+BB": "bottom"}

    EXTENSION_EDGE_WIDTH: dict[str, float] = {"BASE": 1.2, "BASE+BB": 1.8, "BASE+PL+BB": 2.4}

    EXTENSION_LABELS: dict[str, str] = {"BASE": "Base only", "BASE+BB": "Base + BB", "BASE+PL+BB": "Base + PL + BB"}

    # ------------------------------------------------------------------
    # Global defaults (tune once here, propagate everywhere)
    # ------------------------------------------------------------------
    MARKERSIZE: float = 7.0
    LINEWIDTH: float = 1.8
    ERRORBAR_LW: float = 1.0
    CAPSIZE: float = 2.5
    CAPTHICK: float = 1.0
    BAND_ALPHA: float = 0.18  # fill_between confidence bands
    BAND_LW: float = 0.0  # no edge on bands

    # z-order layers for consistent stacking
    ZORDER: dict[str, int] = {
        "band": 1,  # shaded confidence regions — bottom
        "line": 2,  # model curves
        "data": 3,  # data points — top
    }

    # ------------------------------------------------------------------
    # Core kwargs builders
    # ------------------------------------------------------------------

    @classmethod
    def color(cls, grb: str) -> str:
        """Return the hex color for a GRB."""
        return cls.GRB_COLORS[grb]

    @classmethod
    def marker(cls, episode: int) -> str:
        """Return the matplotlib marker code for an episode."""
        return cls.EPISODE_MARKERS[episode]

    @classmethod
    def linestyle(cls, model: str, simple: bool = False):
        """Return the linestyle for a base model.

        Parameters
        ----------
        model  : one of PL, CPL, Band, SBPL
        simple : if True return string (e.g., '--'), else tupled offset form
        """
        if simple:
            return cls.MODEL_LINESTYLES_SIMPLE[model]
        return cls.MODEL_LINESTYLES[model]

    @classmethod
    def marker_kwargs(cls, grb: str, episode: int, extension: str = "BASE", **overrides) -> dict:
        """
        Kwargs for ax.plot() / ax.scatter() marker appearance.

        Returns keys: marker, color, fillstyle, markeredgecolor,
                      markeredgewidth, markersize.
        """
        color = cls.GRB_COLORS[grb]
        fillstyle = cls.EXTENSION_FILLSTYLE[extension]
        mew = cls.EXTENSION_EDGE_WIDTH[extension]

        kw = dict(
            marker=cls.EPISODE_MARKERS[episode],
            color=color,
            fillstyle=fillstyle,
            markeredgecolor=color,
            markeredgewidth=mew,
            markersize=cls.MARKERSIZE,
            zorder=cls.ZORDER["data"],
        )
        kw.update(overrides)
        return kw

    @classmethod
    def line_kwargs(cls, grb: str, model: str, simple_ls: bool = False, **overrides) -> dict:
        """
        Kwargs for ax.plot() line appearance (model curves).

        Returns keys: color, linestyle, linewidth, zorder.
        """
        kw = dict(
            color=cls.GRB_COLORS[grb],
            linestyle=cls.linestyle(model, simple=simple_ls),
            linewidth=cls.LINEWIDTH,
            zorder=cls.ZORDER["line"],
        )
        kw.update(overrides)
        return kw

    @classmethod
    def errorbar_kwargs(cls, grb: str, episode: int, extension: str = "BASE", **overrides) -> dict:
        """
        Full kwargs dict for ax.errorbar().

        Returns keys: fmt, color, fillstyle, markeredgecolor,
                      markeredgewidth, markersize, ecolor, elinewidth,
                      capsize, capthick, zorder.
        """
        color = cls.GRB_COLORS[grb]
        fillstyle = cls.EXTENSION_FILLSTYLE[extension]
        mew = cls.EXTENSION_EDGE_WIDTH[extension]

        kw = dict(
            fmt=cls.EPISODE_MARKERS[episode],
            color=color,
            fillstyle=fillstyle,
            markeredgecolor=color,
            markeredgewidth=mew,
            markersize=cls.MARKERSIZE,
            ecolor=color,
            elinewidth=cls.ERRORBAR_LW,
            capsize=cls.CAPSIZE,
            capthick=cls.CAPTHICK,
            zorder=cls.ZORDER["data"],
        )
        kw.update(overrides)
        return kw

    @classmethod
    def band_kwargs(cls, grb: str, **overrides) -> dict:
        """
        Kwargs for ax.fill_between() confidence/credible bands.

        Returns keys: color, alpha, linewidth, zorder.
        """
        kw = dict(color=cls.GRB_COLORS[grb], alpha=cls.BAND_ALPHA, linewidth=cls.BAND_LW, zorder=cls.ZORDER["band"])
        kw.update(overrides)
        return kw

    @classmethod
    def patch_kwargs(cls, model: str, grb: str | None = None, **overrides) -> dict:
        """
        Kwargs for bar / Rectangle / Patch (e.g., residual histograms).

        Returns keys: hatch, edgecolor, facecolor, linewidth.
        Facecolor defaults to a semi-transparent version of the GRB color
        if grb is provided, else white.
        """
        fc = cls._alpha_hex(cls.GRB_COLORS[grb], 0.25) if grb else "white"
        kw = dict(
            hatch=cls.MODEL_HATCHES[model],
            edgecolor=cls.GRB_COLORS[grb] if grb else "#333333",
            facecolor=fc,
            linewidth=1.0,
        )
        kw.update(overrides)
        return kw

    # ------------------------------------------------------------------
    # Compound helper: all kwargs for a single plotted entity
    # ------------------------------------------------------------------

    @classmethod
    def full_style(cls, grb: str, episode: int, model: str, extension: str = "BASE") -> dict:
        """
        Merge marker and line style for ax.plot() calls where the same
        artist carries both marker (data) and line (model) information.
        """
        kw = cls.marker_kwargs(grb, episode, extension)
        kw.update(linestyle=cls.linestyle(model), linewidth=cls.LINEWIDTH)
        return kw

    # ------------------------------------------------------------------
    # Legend handle factories
    # ------------------------------------------------------------------

    @classmethod
    def grb_handles(cls, grbs: list[str] | None = None) -> list:
        """One colored circle per GRB."""
        grbs = grbs or list(cls.GRB_COLORS)
        return [
            mlines.Line2D(
                [],
                [],
                linestyle="None",
                marker="o",
                color=cls.GRB_COLORS[g],
                markerfacecolor=cls.GRB_COLORS[g],
                markeredgecolor=cls.GRB_COLORS[g],
                markersize=cls.MARKERSIZE,
                label=g,
            )
            for g in grbs
        ]

    @classmethod
    def episode_handles(cls, episodes: list[int] | None = None) -> list:
        """One neutral-colored marker shape per episode."""
        episodes = episodes or [1, 2, 3]
        neutral = "#444444"
        return [
            mlines.Line2D(
                [],
                [],
                linestyle="None",
                marker=cls.EPISODE_MARKERS[e],
                color=neutral,
                markerfacecolor=neutral,
                markeredgecolor=neutral,
                markersize=cls.MARKERSIZE,
                label=cls.EPISODE_LABELS[e],
            )
            for e in episodes
        ]

    @classmethod
    def model_handles(cls, models: list[str] | None = None, color: str = "#222222") -> list:
        """One line-style sample per base model."""
        models = models or ["PL", "CPL", "Band", "SBPL"]
        return [
            mlines.Line2D(
                [], [], linestyle=cls.MODEL_LINESTYLES_SIMPLE[m], linewidth=cls.LINEWIDTH, color=color, label=m
            )
            for m in models
        ]

    @classmethod
    def extension_handles(cls, grb: str, episode: int = 1, extensions: list[str] | None = None) -> list:
        """One marker fill-style sample per extension, colored by GRB."""
        extensions = extensions or ["BASE", "BASE+BB", "BASE+PL+BB"]
        handles = []
        for ext in extensions:
            mkw = cls.marker_kwargs(grb, episode, ext)
            h = mlines.Line2D(
                [],
                [],
                linestyle="None",
                marker=mkw["marker"],
                color=mkw["color"],
                fillstyle=mkw["fillstyle"],
                markeredgecolor=mkw["markeredgecolor"],
                markeredgewidth=mkw["markeredgewidth"],
                markersize=mkw["markersize"],
                label=cls.EXTENSION_LABELS[ext],
            )
            handles.append(h)
        return handles

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _alpha_hex(hex_color: str, alpha: float) -> str:
        """Return hex_color blended toward white by (1-alpha)."""
        h = hex_color.lstrip("#")
        r, g, b = (int(h[i : i + 2], 16) / 255 for i in (0, 2, 4))
        r2 = r * alpha + (1 - alpha)
        g2 = g * alpha + (1 - alpha)
        b2 = b * alpha + (1 - alpha)
        return "#{:02x}{:02x}{:02x}".format(int(r2 * 255), int(g2 * 255), int(b2 * 255))

    @classmethod
    def summary(cls) -> None:
        """Print a compact reference table for all style assignments."""
        W = 55
        sep = "─" * W

        print(f"\n{'GRBPlotStyle — Quick Reference':^{W}}")
        print(sep)

        print("\n  GRB  →  Color (Okabe–Ito)")
        for grb, color in cls.GRB_COLORS.items():
            print(f"    {grb:<14s}  {color}")

        print(f"\n  Episode  →  Marker")
        for ep, mk in cls.EPISODE_MARKERS.items():
            print(f"    {cls.EPISODE_LABELS[ep]:<12s}  '{mk}'")

        print(f"\n  Base Model  →  Line Style  |  Hatch")
        for m in ["PL", "CPL", "Band", "SBPL"]:
            ls = cls.MODEL_LINESTYLES_SIMPLE[m]
            h = cls.MODEL_HATCHES[m] or "—"
            print(f"    {m:<6s}  '{ls}'   hatch='{h}'")

        print(f"\n  Extension  →  Fill Style  |  Edge Width")
        for ext in ["BASE", "BASE+BB", "BASE+PL+BB"]:
            fs = cls.EXTENSION_FILLSTYLE[ext]
            mew = cls.EXTENSION_EDGE_WIDTH[ext]
            print(f"    {ext:<14s}  fillstyle='{fs}'   mew={mew}")

        print(f"\n  Global defaults")
        print(f"    markersize={cls.MARKERSIZE}  linewidth={cls.LINEWIDTH}  " f"band_alpha={cls.BAND_ALPHA}")
        print(sep + "\n")
