from __future__ import annotations

# ===================================================================
# ============================== COLOURS ==============================
# ===================================================================

# ==============================
"""
- ILUPAL colour table of IACOB-GBAT (gbat_ilupal.pro, S. Simon-Diaz 2012),
  indices 0-11 (the table repeats every 12 colours):
  black, white, red, blue, green, orange, pink, purple, light blue and three
  greys (light, medium, dark). It differs from the ILUPAL of IACOB-BROAD.

- In the PostScript device, colour 1 (white) is only used through symbfill,
  where it means the default (black) colour.
"""
_PALETTE = [(0, 0, 0), (255, 255, 255), (255, 0, 0), (0, 0, 255), (0, 155, 0), (255, 105, 0),
            (255, 0, 155), (155, 0, 255), (0, 155, 255), (205, 205, 205), (155, 155, 155), (105, 105, 105)]
_COLORS = [tuple(c / 255. for c in rgb) for rgb in _PALETTE]


# ==============================
"""
- Page of the figures: IDL PostScript landscape device (29 x 20 cm) converted
  to an A4 landscape PDF page.
"""
_PAGE = (29.7 / 2.54, 21.0 / 2.54)




# ==============================
def _color(c : int) -> tuple:
    """
    Colour of index c of the ILUPAL table (c modulo 12).
    """

    # #################### OUTPUT ####################
    return _COLORS[int(c) % len(_COLORS)]




# ======================================================================
# ============================== LAYOUT ==============================
# ======================================================================

# ==============================
def _crea_subplot(
                  m    : list[int],
                  xmid : float,
                  ymid : float,
                  xmin : float = 0.04,
                  xmax : float = 0.98,
                  ymin : float = 0.15,
                  ymax : float = 0.96
    ) -> list[tuple[float, float, float, float]]:
    """
    Positions of panels arranged in rows (translation of
    gbat_crea_subplot.pro).

    ----------
    Parameters
    ----------
        m : list[int]
            Number of panels in each row, from top to bottom.

        xmid, ymid : float
            Horizontal and vertical gaps between panels (figure fraction).

        xmin, xmax, ymin, ymax : float, optional
            Limits of the plotting area (figure fraction), with the IDL
            defaults.

    -------
    Returns
    -------
        positions : list[tuple[float, float, float, float]]
            IDL normal positions (x1, y1, x2, y2) of each panel, row by row
            from the top-left one.

    -----
    Notes
    -----
        - Same geometry as _panel_positions of IACOB-BROAD, but with the IDL
          (x1, y1, x2, y2) convention, needed to combine panels as in
          gbat_plot_dist.pro. Rows without panels keep their space.
    """

    # #################### POSITIONS ####################
    ny = len(m)
    ysh = (ymax - ymin - ymid * (ny - 1)) / ny
    positions = []
    for i, ncol in enumerate(m):
        y1, y2 = ymax - (i + 1) * ysh - i * ymid, ymax - i * ysh - i * ymid
        if ncol <= 0:
            continue
        xsh = (xmax - xmin - xmid * (ncol - 1)) / ncol
        positions += [(xmin + j * xsh + j * xmid, y1, xmin + (j + 1) * xsh + j * xmid, y2) for j in range(ncol)]


    # #################### OUTPUT ####################
    return positions




# ==============================
def _axes(
          fig : plt.Figure,
          pos : tuple[float, float, float, float]
    ) -> plt.Axes:
    """
    Axes at an IDL normal position (x1, y1, x2, y2).
    """

    # #################### OUTPUT ####################
    return fig.add_axes([pos[0], pos[1], pos[2] - pos[0], pos[3] - pos[1]])




# ======================================================================
# ============================== SYMBOLS ==============================
# ======================================================================

# ==============================
_MARKERS = {1: '+', 2: (6, 2, 0), 3: '.', 4: 'D', 5: '^', 6: 's', 7: '*', 8: 'o', 9: (4, 1, 0),
            10: 'v', 11: '^', 12: '>', 13: '<'}




# ==============================
def _symbfill(
              ax   : plt.Axes,
              x    : np.ndarray | float,
              y    : np.ndarray | float,
              isim : int,
              tsim : float,
              csim : int,
              fill : bool = False,
              **kwargs
    ) -> None:
    """
    Coloured (filled) symbols (translation of gbat_symbfill.pro).

    ----------
    Parameters
    ----------
        ax : plt.Axes
            Axes of the panel.

        x, y : np.ndarray | float
            Coordinates.

        isim : int
            Symbol: 1 cross, 2 asterisk, 3 dot, 4 romb, 5 triangle, 6 square,
            7 star1, 8 circle, 9 star2, 10-13 upper/lower/right/left limits.

        tsim : float
            Symbol size (IDL units).

        csim : int
            ILUPAL colour (1 -> default black colour).

        fill : bool, optional
            Filled symbol. Defaults to False (outline only).

    -----
    Notes
    -----
        - The IDL user symbols are scaled by tsim twice (coordinates and
          SYMSIZE), so the Matplotlib marker size grows as tsim^2.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### SYMBOL ####################
    color = 'black' if csim == 1 else _color(csim)
    marker = _MARKERS.get(isim, 'o')
    ms = 5.5 * tsim**2
    style = dict(ls='none', marker=marker, ms=ms)
    if fill or isim in (1, 2, 3):
        style.update(color=color, mec=color, mfc=color)
    else:
        style.update(color=color, mec=color, mfc='none', mew=0.8)
    style.update(kwargs)


    # #################### PLOT ####################
    ax.plot(np.atleast_1d(x), np.atleast_1d(y), **style)




# ===========================================================================
# ============================== LINE PROFILES ==============================
# ===========================================================================

# ==============================
def _plot_prof(
               ax         : plt.Axes,
               labels     : dict,
               sigma      : float,
               sigma_real : float,
               w0         : np.ndarray | None,
               f0         : np.ndarray | None,
               wo0        : np.ndarray[float],
               fo0        : np.ndarray[float],
               pe0        : np.ndarray[float]
    ) -> None:
    """
    Panel of one line in the profiles figure: observed spectrum and
    convolved FASTWIND profile of the best-fitting model (plotting part of
    gbat_plot_prof.pro).

    ----------
    Parameters
    ----------
        ax : plt.Axes
            Axes of the panel.

        labels : dict
            Label and limits of the line (see iacob_gbat_input._lineshhe).

        sigma : float
            Noise of the line from the SNR of the continuum.

        sigma_real : float
            Noise of the line from the best-fitting model.

        w0, f0 : np.ndarray | None
            Convolved synthetic profile (None if the model is not available).

        wo0, fo0, pe0 : np.ndarray[float]
            Observed spectrum and mask of the line points.

    -----
    Notes
    -----
        - All the observed points are drawn in light blue and the points used
          in the fit in red; the synthetic profile is the black line.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _fi, _idl2mpl


    # #################### LIMITS ####################
    la, lb = labels['la'], labels['lb']
    nwo0 = (wo0 >= la) & (wo0 <= lb) & (pe0 > 0)
    values = [fo0[nwo0]] if nwo0.any() else []
    if w0 is not None:
        values.append(f0[(w0 >= la) & (w0 <= lb)])
    values = np.concatenate(values) if values else np.array([0., 1.])
    ymin, ymax = values.min(), values.max()
    xr = [int(la), int(lb)]


    # #################### PLOT ####################
    inside = (wo0 >= xr[0] - 1) & (wo0 <= xr[1] + 1)
    _symbfill(ax, wo0[inside], fo0[inside], 8, 0.6, 8, fill=True)
    sel = inside & (pe0 > 0)
    _symbfill(ax, wo0[sel], fo0[sel], 8, 0.6, 2, fill=True)
    if w0 is not None:
        ax.plot(w0, f0, color='black', lw=1.)


    # #################### AXES ####################
    ax.set(xlim=xr, ylim=[ymin - 0.5 * (ymax - ymin), ymax + 0.3 * (ymax - ymin)])
    ax.set_xticks([xr[0], (xr[0] + xr[1]) / 2., xr[1]])
    ax.xaxis.set_major_formatter(lambda v, _: f'{v:.0f}')
    ax.set_title(_idl2mpl(labels['li']), fontsize=9, pad=2)
    ax.set_xlabel(r'$\lambda$ (A)', fontsize=7, labelpad=1)
    ax.tick_params(labelsize=6.5, direction='in', top=True, right=True, length=3)


    # #################### LABELS ####################
    dx = (xr[1] - xr[0]) * 0.03
    ax.text(xr[0] + dx, ymax + 0.12 * (ymax - ymin), 'SNR$_c$=' + _fi(1. / sigma, 3).strip(), fontsize=6)
    ax.text(xr[0] + dx, ymin - 0.35 * (ymax - ymin), r'1/$\sigma_{BFM}$=' + _fi(1. / sigma_real, 3).strip(), fontsize=6)




# ===============================================================================
# ============================== CHI2 DISTRIBUTIONS ==============================
# ===============================================================================

# ==============================
def _plot_dist_panel(
                     ax      : plt.Axes,
                     var     : np.ndarray[float],
                     chi     : np.ndarray[float],
                     varx_sp : np.ndarray[float],
                     vary_sp : np.ndarray[float],
                     var_x   : np.ndarray[float],
                     var_y   : np.ndarray[float],
                     title   : str,
                     fw2     : np.ndarray[bool] | None = None,
                     fw3     : np.ndarray[bool] | None = None
    ) -> None:
    """
    Panel with the chi2 distribution of one parameter (plotting part of
    gbat_plot_dist_one.pro).

    -----
    Notes
    -----
        - Grey stars: all the models; black line and circles: lower envelope
          and its spline fit; dashed red lines: min(chi2) + 0, 1 and 4; red
          circle: best-fitting model. fw2 and fw3 (masks of two extra models)
          are marked with pink and purple stars.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from matplotlib.ticker import MaxNLocator
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _idl2mpl


    # #################### LIMITS ####################
    var_min, var_max = var.min(), var.max()
    cmin = chi.min()
    xsize = var_max - var_min
    if xsize == 0:
        xsize = 10. * var_min


    # #################### PLOT ####################
    ax.set(xlim=[var_min - 0.07 * xsize, var_max + 0.07 * xsize], ylim=[0.8 * cmin, cmin + 5])
    _symbfill(ax, var, chi, 8, 0.7, 9, fill=True)
    ax.plot(varx_sp, vary_sp, color='black', lw=1.3)
    _symbfill(ax, var_x, var_y, 8, 0.8, 1, fill=True)
    for mask, c in ((fw2, 6), (fw3, 7)):
        if mask is not None and np.any(mask):
            _symbfill(ax, var[mask], chi[mask], 7, 1., c, fill=True)
    for dy in (0., 1., 4.):
        ax.axhline(cmin + dy, ls='--', color=_color(2), lw=1.3)
    _symbfill(ax, var[0], chi[0], 8, 1., 2, fill=True, mec='black', mew=0.8)


    # #################### AXES ####################
    ax.set_xlabel(_idl2mpl(title), fontsize=8, labelpad=1)
    ax.tick_params(labelsize=7, direction='in', top=True, right=True, length=3)
    ax.xaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))




# ================================================================================
# ============================== EVOLUTIONARY TRACKS ==============================
# ================================================================================

# ==============================
def _plot_gva(
              ax   : plt.Axes,
              diag : str = 'gtd',
              mli  : tuple[float, float] | None = None,
              xli  : tuple[float, float] | None = None,
              yli  : tuple[float, float] | None = None,
              te   : bool = False,
              line : int = 3
    ) -> None:
    """
    Non-rotating Geneva evolutionary tracks (Z = 0.014) in the Kiel (gtd),
    HR (hrd) or spectroscopic HR (shrd) diagram (translation of
    gbat_plot_gva.pro).

    ----------
    Parameters
    ----------
        ax : plt.Axes
            Axes of the panel.

        diag : str, optional
            'gtd' (log g), 'hrd' (log L) or 'shrd' (log L_spec).

        mli : tuple[float, float] | None, optional
            Range of initial masses [Msun].

        xli, yli : tuple[float, float] | None, optional
            Limits of the axes (IDL defaults if None).

        te : bool, optional
            Teff [kK] on the x-axis instead of log Teff.

        line : int, optional
            IDL line style + 1 of the tracks (3 -> dashed).

    -----
    Notes
    -----
        - The first 210 points of each track (up to the beginning of the He
          burning) are drawn, together with the ZAMS (thick black line) and
          the initial masses.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import numpy as np

    # ########## PANORAMA ##########
    from tracks import trackgene12
    from gbat.io import _readcol


    # #################### TRACKS ####################
    m0 = ['003', '004', '005', '007', '009', '012', '015', '020', '025', '032', '040', '060', '085']
    m1 = ['3', '4', '5', '7', '9', '12', '15', '20', '25', '32', '40', '60', '85']
    if mli is not None:
        keep = [i for i, m in enumerate(m0) if mli[0] <= float(m) <= mli[1]]
        m0, m1 = [m0[i] for i in keep], [m1[i] for i in keep]


    # #################### AXES ####################
    ytit = {'gtd': 'logg', 'hrd': r'logL/L$_{sol}$', 'shrd': r'log $\mathcal{L}/\mathcal{L}_{\odot}$'}[diag]
    if yli is None:
        yli = {'gtd': (4.5, 2.0), 'hrd': (1.7, 6.5), 'shrd': (2.4, 4.7)}[diag]
    xtit = r'T$_{eff}$ [kK]' if te else r'log(T$_{eff}$ [K])'
    if xli is None:
        xli = (57., 9.) if te else (4.8, 3.9)
    ax.set(xlim=xli, ylim=yli)
    ax.set_xlabel(xtit, fontsize=8, labelpad=1)
    ax.set_ylabel(ytit, fontsize=8, labelpad=1)
    ax.tick_params(labelsize=7, direction='in', top=True, right=True, length=3)


    # #################### PLOT ####################
    x0, y0 = [], []
    xsh, ysh = 0.02 * (max(xli) - min(xli)), 0.02 * (max(yli) - min(yli))
    style = {1: '-', 2: ':', 3: '--', 4: '-.'}.get(line, '--')
    for m in m0:
        track = trackgene12(float(m), vini_vcrit=0.0)
        y = {'gtd': track['log_g'], 'hrd': track['log_L'], 'shrd': track['log_Lspec']}[diag]
        x = 10.**track['log_Teff'] * 1e-3 if te else track['log_Teff']
        x0.append(x[0])
        y0.append(y[0])
        ax.plot(x[:210], y[:210], color=_color(10), lw=1.5, ls=style)
    if x0:
        ax.plot(x0, y0, color='black', lw=2.)
        for x, y, lab in zip(x0, y0, m1):
            ax.text(x + xsh, y - ysh, lab, ha='right', color=_color(11), fontsize=8, clip_on=True)
    if diag == 'shrd':
        ax.text(4.75, 3.15, r'$\mathcal{L}$:=T$_{eff}^{4}$/g', fontsize=8)




# ==========================================================================
# ============================== CHI2 MAPS ==============================
# ==========================================================================

# ==============================
def _plot_map(
              ax     : plt.Axes,
              chi2   : np.ndarray[float],
              chi2_x : np.ndarray[float],
              chi2_y : np.ndarray[float],
              xtit   : str = ' ',
              ytit   : str = ' ',
              xticks : bool = True,
              yticks : bool = True
    ) -> None:
    """
    Filled contours of chi2 - min(chi2) at 1, 2 and 3 sigma in the plane of
    two parameters (translation of gbat_plot_map.pro).

    -----
    Notes
    -----
        - Irregular grid (IDL CONTOUR /IRREGULAR) -> Matplotlib triangulation.
          If the points cannot be triangulated (e.g. a parameter with a single
          value), only the points within 3 sigma are drawn.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### RANGES ####################
    chi0 = chi2 - chi2.min()
    rx = [chi2_x.min(), chi2_x.max()]
    ry = [chi2_y.min(), chi2_y.max()]
    ax.set(xlim=[rx[0] - 0.1 * (rx[1] - rx[0]) - 1e-9, rx[1] + 0.1 * (rx[1] - rx[0]) + 1e-9],
           ylim=[ry[0] - 0.1 * (ry[1] - ry[0]) - 1e-9, ry[1] + 0.1 * (ry[1] - ry[0]) + 1e-9])


    # #################### CONTOURS ####################
    levels = chi0.min() + np.array([0., 1., 4., 9.])
    try:
        ax.tricontourf(chi2_x, chi2_y, chi0, levels=[*levels, np.inf],
                       colors=[_color(9), _color(10), _color(11), _color(1)])
        cs = ax.tricontour(chi2_x, chi2_y, chi0, levels=levels, colors='black', linewidths=0.8)
        ax.clabel(cs, levels=levels[1:], fmt=dict(zip(levels[1:], [r'1$\sigma$', r'2$\sigma$', r'3$\sigma$'])),
                  fontsize=6)
    except (ValueError, RuntimeError):
        near = chi0 <= 9
        ax.plot(chi2_x[near], chi2_y[near], '.', color=_color(10))


    # #################### AXES ####################
    ax.set_xlabel(xtit, fontsize=7)
    ax.set_ylabel(ytit, fontsize=7)
    ax.tick_params(labelsize=6, direction='in', top=True, right=True, labelbottom=xticks, labelleft=yticks)




# ==========================================================================
# ============================== INTERACTIVE ==============================
# ==========================================================================
"""
- _xhair copied from iacob_broad_interactive of IACOB-BROAD (gbat_xhair.pro).
"""

# ==============================
def _xhair(
           ax : plt.Axes,
           n  : int = 1
    ) -> list[tuple[float, float]]:
    """
    Collect data coordinates from interactive mouse clicks.

    This function is used for interactively selecting n points 
    on a Matplotlib axis. The returned coordinates are expressed in 
    the data coordinate system of the axis.

    ----------
    Parameters
    ----------
        ax : matplotlib.axes.Axes
            Matplotlib axis on which the points are selected.

        n : int, optional
            Number of points to select. Defaults to 1.

    -------
    Returns
    -------
        coordinates : list(tuple[float, float])
            List containing n (x, y) pairs in data coordinates,
            in the order in which they were selected.

    -----
    Notes
    -----
        - Left clicks are collected with a custom 'button_press_event' handler
          (not matplotlib.pyplot.ginput), so the function waits until the
          requested number of points has been selected.

        - While the toolbar's zoom or pan tool is active, clicks are left to
          that tool (panning/zooming the axes) and are not read as a selected
          point; points are only read while zoom/pan is inactive.

        - The returned coordinates correspond to the axis data coordinates,
          regardless of the figure size, zoom level, or current axis limits.

        - The mouse position is marked with a red crosshair (a vertical and
          a horizontal line spanning the full axes), always visible (not only
          while the mouse is moving) once it has entered ax.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import matplotlib.pyplot as plt


    # #################### FUNCTIONS ####################
    def _on_move(event):
        if event.inaxes is not ax or event.xdata is None:
            return
        crosshair_v.set_xdata([event.xdata, event.xdata])
        crosshair_h.set_ydata([event.ydata, event.ydata])
        crosshair_v.set_visible(True)
        crosshair_h.set_visible(True)
        fig.canvas.draw_idle()

    def _on_click(event):
        toolbar = fig.canvas.toolbar
        if toolbar is not None and toolbar.mode: # Zoom/pan tool active: let it handle the click
            return
        if event.inaxes is not ax or event.button != 1 or event.xdata is None: # Only left clicks inside ax
            return
        points.append((event.xdata, event.ydata))
        markers.append(ax.plot(event.xdata, event.ydata, '+', color='red', ms=10, mew=1.5)[0])
        fig.canvas.draw_idle()


    # #################### INTERACTIVE POINT SELECTION ####################
    fig = ax.figure

    crosshair_v = ax.axvline(ax.get_xlim()[0], color='red', lw=1, visible=False)
    crosshair_h = ax.axhline(ax.get_ylim()[0], color='red', lw=1, visible=False)
    points = []
    markers = []

    cid_move = fig.canvas.mpl_connect('motion_notify_event', _on_move)
    cid_click = fig.canvas.mpl_connect('button_press_event', _on_click)
    try:
        while len(points) < n:
            if not plt.fignum_exists(fig.number):
                raise RuntimeError('Figure closed before selecting the requested points.')
            plt.pause(0.05)
    finally:
        fig.canvas.mpl_disconnect(cid_move)
        fig.canvas.mpl_disconnect(cid_click)
        crosshair_v.remove()
        crosshair_h.remove()
        fig.canvas.draw_idle()
    coordinates = [(float(x), float(y)) for x, y in points]


    # #################### OUTPUT ####################
    return coordinates
