from __future__ import annotations

# ======================================================================
# ============================== LINELIST ==============================
# ======================================================================

# ==============================
def _read_linelist() -> dict[str, float]:
    """
    Reading the line list from file 'pyiacob_broad_lines.dat'.

    ----------
    Parameters
    ----------
        None

    -------
    Returns
    -------
        lines_dict : dict[str, float]
            Dictionary with line names and their central wavelengths [A].
            Lines starting with '#' are comments, and are ignored.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .variables import LINE_LIST


    # #################### READING LINE LIST ####################
    lines_table = np.genfromtxt(LINE_LIST, dtype=None, comments='#', encoding='utf-8', usecols=(0, 1))
    lines_dict = {str(name): float(lam) for name, lam in lines_table}


    # #################### OUTPUT ####################
    return lines_dict




# =================================================================================
# ============================== INTERACTIVE PLOTS 1 ==============================
# =================================================================================

# ==============================
"""
- ILUPAL colour table (I. Monteverde 1998, S. Simon-Diaz 2012), indices 0-13:
  black, white, red, blue, green, orange, pink, purple, light blue, greys, violet
"""
_PALETTE = [(0, 0, 0), (255, 255, 255), (255, 0, 0), (0, 0, 255), (0, 155, 0), (255, 165, 0),
           (255, 153, 255), (155, 0, 255), (0, 155, 255), (240, 240, 240), (205, 205, 205),
           (155, 155, 155), (105, 105, 105), (220, 220, 255)]
_COLORS = [tuple(c / 255. for c in rgb) for rgb in _PALETTE]
_DASH_DOT3 = (0, (3, 1, 1, 1, 1, 1))




# ==============================
def _window() -> plt.Figure:
    """
    Create and initialise the interactive preprocessing window.

    -------
    Returns
    -------
        matplotlib.figure.Figure
            Initialised Matplotlib figure used for the interactive
            preprocessing steps.

    -----
    Notes
    -----
        - The figure size is specified in inches assuming the Matplotlib
          default resolution of 100 dpi:
            1120 × 250 pixels = 11.2 × 2.5 inches.

        - Calling this function repeatedly reuses figure number 2 and clears
          its previous contents.
    """

    # #################### MODULES & PACKAGES ####################
    import matplotlib.pyplot as plt


    # #################### FIGURE INITIALISATION ####################
    fig = plt.figure(num=2) # Create or reuse the preprocessing figure
    fig.set_size_inches(11.2, 2.5) # Set the figure size to approximately 1120 × 250 pixels at 100 dpi.
    fig.clf() # Clear any axes and artists from a previous preprocessing step.


    # #################### OUTPUT ####################
    return fig




# ==============================
def _show(
          fig    : plt.Figure,
          title  : str,
          curves : Sequence[Tuple],
          lines  : np.ndarray[float] | None = None,
          xlim   : np.ndarray[float] | None = None,
          ylim   : np.ndarray[float] | None = None,
          ax     : plt.Axes | None = None
    ) -> plt.Axes:
    """
    Display spectra and spectral-line markers in an interactive window.

    This function prepares a Matplotlib axis for the interactive
    preprocessing steps. Spectral curves are plotted, optional
    reference-line wavelengths are marked with vertical dashed lines,
    and the axis limits and title are updated. The figure canvas is
    explicitly refreshed so that the plot is immediately available for
    subsequent interactive input.

    ----------
    Parameters
    ----------
        fig : matplotlib.figure.Figure
            Figure used for the interactive display.

        title : str
            Title or instruction displayed above the plot.

        curves : sequence of tuple[numpy.ndarray[float], numpy.ndarray[float], dict]
            Sequence of (x, y, kwargs) tuples defining the curves to
            plot. kwargs contains optional Matplotlib plotting
            arguments for each curve.

        lines : np.ndarray[float] | None, optional
            Wavelengths of reference spectral lines. Each wavelength is
            displayed as a vertical dashed line.

        xlim : np.ndarray[float] | None, optional
            Lower and upper limits of the x-axis. If None, Matplotlib
            determines the limits automatically.

        ylim : np.ndarray[float] | None, optional
            Lower and upper limits of the y-axis. If None, Matplotlib
            determines the limits automatically.

        ax : matplotlib.axes.Axes, | None, optional
            Existing axis on which to draw the curves. If None, the
            figure is cleared and a new single-panel axis is created.

    -------
    Returns
    -------
        matplotlib.axes.Axes
            Axis containing the displayed spectra and reference lines.
    """

    # #################### MODULES & PACKAGES ####################
    import matplotlib.pyplot as plt


    # #################### PLOTTING FIGURE ####################
    # ########## AXIS INITIALISATION ##########
    if ax is None: # Create a fresh single-panel axis when no existing axis is supplied.
        fig.clf()
        ax = fig.add_subplot(111)


    # ########## SPECTRUM PLOTTING ##########
    for x, y, kwargs in curves: # Plot each curve using the standard preprocessing line style.
        plot_kwargs = {'color': 'black', 'lw': 0.6, **kwargs,}
        ax.plot(x, y, **plot_kwargs)


    # ########## REFERENCE LINE MARKERS ##########
    if lines is not None:  # Mark the wavelengths from the line list with vertical dashed lines.
        for wavelength in lines:
            ax.axvline(wavelength, ls='--', color='black', lw=0.6)


    # ########## AXIS CONFIGURATION ##########
    ax.set(title=title, xlim=xlim, ylim=ylim) # Apply the requested title and axis limits.


    # #################### FIGURE UPDATE ####################
    fig.canvas.draw_idle() # Force the interactive backend to render the updated figure before returning control to the calling interactive routine.
    plt.pause(0.05)


    # #################### OUTPUT ####################
    return ax




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




# =================================================================================
# ============================== INTERACTIVE PLOTS 2 ==============================
# =================================================================================

# ==============================
def _limits(
            values : list[float] | None,
            ax     : matplotlib.axes.Axes,
            n      : int = 2
    ) -> np.ndarray[float]:
    """
    Return up to n sorted limit values.

    If values is not provided, the limits are obtained interactively
    from n cursor clicks on the given Matplotlib axis using
    _xhair. The selected values are then converted to a NumPy
    array and sorted in ascending order, reproducing the ordering
    expected from the IDL MINMAX convention.

    ----------
    Parameters
    ----------
        values : list[float] | None
            User-supplied limit values. If None, the values are obtained
            interactively from ax using _xhair(ax, n).

        ax : matplotlib.axes.Axes
            Matplotlib axis used for interactive selection when values
            is None.

        n : int, optional
            Maximum number of limit values to use. Defaults to 2.

    -------
    Returns
    -------
        limits : numpy.ndarray[float]
            One-dimensional array containing the first n limit values,
            sorted in ascending order.

    -----
    Notes
    -----
        - When values is None, _xhair(ax, n) returns the (x, y)
          coordinates of each cursor selection. Only the selected
          x-coordinates are used as limits.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### LIMITS ####################
    if values is None:
        values = [x for x, _ in _xhair(ax, n)]


    # #################### OUTPUT ####################
    return np.sort(np.asarray(values)[:n])




# ==============================
def _inside(
            w   : np.ndarray[float],
            lim : np.ndarray[float]
    ) -> np.ndarray[bool]:
    """
    Return a boolean mask selecting values within a given interval.

    The interval is closed, so values equal to either boundary are
    included. The returned mask has the same shape as w and can be
    used directly for boolean indexing.

    ----------
    Parameters
    ----------
        w : numpy.ndarray[float]
            Array of values to test against the specified limits.

        lim : array-like
            Two-element sequence defining the lower and upper limits of the
            interval. The interval is defined as [lim[0], lim[1]].

    -------
    Returns
    -------
        bool_array : numpy.ndarray[bool]
            Boolean array with the same shape as w. An element is
            True when the corresponding value of w satisfies
                lim[0] <= w <= lim[1]
            and False otherwise.
    """

    # #################### COMPUTATION ####################
    bool_array = (w >= lim[0]) & (w <= lim[1]) # Select values within the closed interval [lim[0], lim[1]].


    # #################### OUTPUT ####################
    return bool_array




# ==============================
def _clip_curves(
                 w0    : np.ndarray[float],
                 f0    : np.ndarray[float],
                 wx0   : np.ndarray[float],
                 fx0   : np.ndarray[float],
                 cut   : np.ndarray[bool] | None,
                 shift : float = 0.
    ) -> list[tuple]:
    """
    Curves of a (possibly clipped) line for _show.

    ----------
    Parameters
    ----------
        w0 : numpy.ndarray[float]
            Wavelength array of the original spectrum.

        f0 : numpy.ndarray[float]
            Flux array of the original spectrum.

        wx0 : numpy.ndarray[float]
            Wavelength array of the line (FT and GOF).

        fx0 : numpy.ndarray[float]
            Flux array of the line, with the clipped regions already replaced
            by their interpolation.

        cut : numpy.ndarray[bool] | None
            Mask of the clipped points of wx0. None -> no clipping.

        shift : float, optional
            Constant added to every flux (e.g. to undo the continuum
            normalisation in the plot). Defaults to 0.

    -------
    Returns
    -------
        curves : list[tuple]
            (x, y, kwargs) tuples: original spectrum (black line, without the
            clipped regions), unclipped points of the line (red stars) and
            interpolated points (blue stars).

    -----
    Notes
    -----
        - A point of w0 is hidden if it lies between two clipped points of wx0
          (linear interpolation of the mask > 0.5).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### CURVES ####################
    if cut is None or not np.any(cut):
        cut = np.zeros(wx0.size, dtype=bool)
    hidden = np.interp(w0, wx0, cut.astype(float), left=0., right=0.) > 0.5 # Original data inside the clipped regions
    f0_plot = np.where(hidden, np.nan, f0 + shift)
    curves = [
              (w0, f0_plot, {}),
              (wx0[~cut], fx0[~cut] + shift, dict(color=_COLORS[2], ls='none', marker='*')),
              (wx0[cut], fx0[cut] + shift, dict(color=_COLORS[3], ls='none', marker='*')),
             ]


    # #################### OUTPUT ####################
    return curves




# ==============================
def _clip_line(
               fig   : matplotlib.figure.Figure,
               w0    : np.ndarray[float],
               f0    : np.ndarray[float],
               wx0   : np.ndarray[float],
               fx0   : np.ndarray[float],
               lines : list[tuple[str, float]],
               xr0   : list[float],
               xc    : float,
               ycont : float
    ) -> tuple[np.ndarray[float], np.ndarray[bool] | None]:
    """
    Interactively clip selected regions of a spectral line.

    The user selects a wavelength interval to be clipped. The flux
    values inside the selected region are replaced by a straight line
    joining the nearest unmasked points immediately outside the region.
    The corresponding points are also identified by a boolean mask so
    that they can be excluded later from the goodness-of-fit analysis.

    ----------
    Parameters
    ----------
        fig : matplotlib.figure.Figure
            Figure used for displaying the spectrum during the interactive
            clipping procedure.

        w0 : numpy.ndarray[float]
            Wavelength array of the original spectrum.

        f0 : numpy.ndarray[float]
            Flux array of the original spectrum.

        wx0 : numpy.ndarray[float]
            Wavelength array of the spectrum used for the Fourier-transform
            analysis and clipping.

        fx0 : numpy.ndarray[float]
            Flux array corresponding to wx0. This array is modified in
            the selected clipping regions.

        lines : list[float]
            Wavelengths of the reference lines passed to _show.

        xr0 : list[float]
            Initial x-axis limits passed to _show, only used if fig has no
            axes. Otherwise the wavelength range of the figure on screen (that
            of the previous step, zoom included) is kept in every redraw.

        xc : float
            Wavelength of the line core [A], marked with a blue dashed
            vertical line on every redraw of the plot.

        ycont : float
            Continuum level, marked with a blue dashed horizontal line on
            every redraw of the plot.

    -------
    Returns
    -------
        fx0 : numpy.ndarray[float]
            Flux array with the selected regions replaced by linear
            interpolations between the nearest points outside each region.

        cut : numpy.ndarray[bool] or None
            Boolean mask identifying the points clipped in all the selected
            regions. Returns None if no region was clipped.

    -----
    Notes
    -----
        - The clipping limits are selected interactively using _limits.
          Only points satisfying
            lim[0] < wx0 < lim[1]
          are modified.

        - The replacement flux is obtained from the straight line connecting
          the closest (unclipped) points outside the clipped region:
            wx0 <= lim[0] and wx0 >= lim[1].

        - The user is asked again after every clip, until the answer is 0.

        - In every redraw, the interpolated points are drawn in blue and the
          real data of the clipped regions are not shown (see _clip_curves).
    """

    # #################### MODULES & PACKAGES ####################
    import matplotlib.pyplot as plt
    import numpy as np


    # #################### INITIALISATION ####################
    cut = np.zeros(wx0.size, dtype=bool) # Accumulated mask of clipped points (every region)
    xlim = lambda: fig.axes[-1].get_xlim() if fig.axes else xr0 # Wavelength range of the figure on screen (kept in every redraw, zoom included)

    def mark_core_cont(ax): # Re-mark the core & continuum on every redraw
        ax.axvline(xc, color=_COLORS[3], ls='--')
        ax.axhline(ycont, color=_COLORS[3], ls='--')
        ax.figure.canvas.draw_idle()
        plt.pause(0.05)

    def ask(): # Clip (1) or stop (0); asked again until a valid answer is given
        while True:
            answer = input('Do you want to clip part of the line? (0/1)\n: ').strip()
            if answer in ('0', '1'):
                return answer == '1'
            print(f"Invalid answer '{answer}'. Type 0 (no) or 1 (yes).")

    def redraw(title): # Spectrum with the clipped regions replaced by their interpolation (blue)
        ax = _show(fig, title, _clip_curves(w0, f0, wx0, fx0, cut), lines, xlim())
        mark_core_cont(ax)
        return ax


    # #################### CLIPPING LOOP ####################
    """
    - The question is repeated after every clip until the answer is 0.
    - Every clipped region is accumulated in 'cut' (all of them are removed for the GOF).
    """
    redraw('Do you want to clip part of the line? (see terminal)')
    while ask():
        # ########## SELECT CLIPPING REGION ##########
        ax = redraw('Mark the region of the line you want to clip')
        lim = _limits(None, ax)

        # ########## CLIPPING LIMITS ##########
        inside = (wx0 > lim[0]) & (wx0 < lim[1]) # Points strictly inside the selected interval
        left = np.nonzero((wx0 <= lim[0]) & ~cut)[0] # Closest unclipped points outside the region
        right = np.nonzero((wx0 >= lim[1]) & ~cut)[0]
        if not inside.any() or left.size == 0 or right.size == 0:
            print('Invalid region: it must contain points of the line and have unclipped points on both sides. Try again.')
            redraw('Do you want to clip part of the line? (see terminal)')
            continue
        i1, i2 = left[-1], right[0]

        # ########## REPLACE CLIPPED REGION ##########
        fx0 = fx0.copy() # Straight line between the closest points outside the region
        fx0[inside] = fx0[i1] + (fx0[i2] - fx0[i1]) * (wx0[inside] - wx0[i1]) / (wx0[i2] - wx0[i1])
        cut |= inside

        # ########## DISPLAY MODIFIED SPECTRUM ##########
        redraw('Do you want to clip another part of the line? (see terminal)')


    # #################### OUTPUT ####################
    return fx0, (cut if cut.any() else None)