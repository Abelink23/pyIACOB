from __future__ import annotations

# ===========================================================================
# ============================== PRE-PROCESSING ==============================
# ===========================================================================

# ==============================
class _CleanPlot:
    """
    Interactive window of the pre-processing of a line (IDL window 0 of
    gbat_clean_line.pro). It is only opened when an answer or a mark of the
    user is needed, so that runs with pre-processed spectra ('nocut') do not
    need a display.
    """

    # ==============================
    def __init__(self, title : str, xl : np.ndarray, ref : tuple | None):
        self.title, self.xl, self.ref = title, xl, ref
        self.fig = self.ax = None
        self.box = None # (wmin, wmax, cont_flux, sigma) of the SNR region


    # ==============================
    def draw(
             self,
             wo0       : np.ndarray,
             fo0       : np.ndarray,
             pe0       : np.ndarray | None = None,
             title     : str | None = None
        ) -> plt.Axes:
        """
        (Re)draw the window: spectrum, continuum (dotted red), reference model
        (blue), points of the line (red) and SNR region (blue box).
        """

        # #################### MODULES & PACKAGES ####################
        # ########## PYTHON ##########
        import matplotlib.pyplot as plt
        import numpy as np


        # ########## PANORAMA ##########
        from gbat.io import _fi, _idl2mpl
        from gbat.plots import _color


        # #################### FIGURE ####################
        if self.fig is None or not plt.fignum_exists(self.fig.number):
            self.fig = plt.figure(num='IACOB-GBAT', figsize=(6., 4.))
        xlim = self.ax.get_xlim() if self.ax is not None and self.ax.figure is self.fig else self.xl
        self.fig.clf()
        ax = self.ax = self.fig.add_subplot(111)


        # #################### SPECTRUM ####################
        ss = (wo0 >= self.xl[0]) & (wo0 <= self.xl[1])
        ylim = [0.9 * fo0[ss].min(), min(1.1 * fo0[ss].max(), 2.)] if ss.any() else None
        ax.plot(wo0, fo0, color='black', lw=0.8)
        ax.axhline(1., ls=':', color=_color(2))
        if self.ref is not None:
            ax.plot(*self.ref, color=_color(3))
        if pe0 is not None:
            ax.plot(wo0, np.where(pe0 > 0, fo0, np.nan), color=_color(2), lw=1.)


        # #################### SNR REGION ####################
        if self.box is not None:
            w1, w2, cont, sigma = self.box
            ax.plot([w1, w1, w2, w2, w1], cont + sigma * np.array([-1, 1, 1, -1, -1]), ls='--', color=_color(3))
            ax.text((w1 + w2) / 2., cont + 3. * sigma, 'SNR=' + _fi(cont / sigma, 4).strip(), ha='center')


        # #################### AXES ####################
        ax.set(xlim=xlim, ylim=ylim, title=_idl2mpl(title or self.title), xlabel=r'$\lambda$ (A)',
               ylabel='Normalized flux')
        self.fig.canvas.draw_idle()
        plt.pause(0.05)


        # #################### OUTPUT ####################
        return ax


    # ==============================
    def marks(self, n : int = 2) -> list[tuple[float, float]]:
        """
        n points marked with the mouse (gbat_xhair.pro -> _xhair of
        IACOB-BROAD).
        """

        # #################### PANORAMA ####################
        from gbat.plots import _xhair


        # #################### OUTPUT ####################
        return _xhair(self.ax, n)




# ==============================
def _clean_line(
                wo0      : np.ndarray[float],
                fo0      : np.ndarray[float],
                pe0      : np.ndarray[float],
                sig0     : np.ndarray[float],
                xl       : np.ndarray[float],
                line     : str,
                vsini    : float,
                cfg      : Config,
                zeta_t   : float = 0.,
                vgauss   : float = 0.,
                snr      : float | None = None,
                nocut    : bool = False,
                fw_ref   : str | None = None,
                interactive : bool = False
    ) -> float:
    """
    Pre-processing of a line (translation of gbat_clean_line.pro):
    renormalisation of a region, SNR of the adjacent continuum, limits of the
    line and cleaning/restoring of regions of the line.

    ----------
    Parameters
    ----------
        wo0, fo0, pe0, sig0 : np.ndarray[float]
            Observed spectrum: wavelength [A], normalised flux, mask of the
            points of the line used in the fit (0/1) and mask of the points
            used to compute the SNR (0/1). fo0, pe0 and sig0 are modified in
            place.

        xl : np.ndarray[float]
            Limits of the window of the line [A].

        line : str
            FASTWIND label of the line.

        vsini : float
            Projected rotational velocity [km/s] (only for the reference
            model).

        cfg : Config
            Configuration of the run.

        zeta_t, vgauss : float, optional
            Macroturbulent velocities [km/s] (only for the reference model).

        snr : float | None, optional
            Global SNR. If given, sigma = 1/snr and the whole window is marked
            as SNR region.

        nocut : bool, optional
            The spectrum is pre-processed: do not ask for re-processing.

        fw_ref : str | None, optional
            FASTWIND model of the grid used as reference (micro '005').

        interactive : bool, optional
            Ask, for every line, whether to renormalise a region of the line,
            whether to (re)compute the SNR from a region of the continuum and
            whether to mark again the limits of the line (all marked with the
            mouse), even with 'nocut' or a global 'snr'. Default False (IDL
            behaviour: renormalisation and SNR only asked without 'nocut', and
            the limits only when none are stored).

    -------
    Returns
    -------
        sigma : float
            Noise of the line (1/snr or standard deviation of the SNR region).

    -----
    Notes
    -----
        - As in IDL, even with 'nocut' the limits of the SNR region and of the
          line are asked for when they are not stored in the spectrum. With
          'interactive', stored limits of the line are shown and the user is
          asked whether to mark them again (the old ones are then removed
          inside the window of the line).
        - Renormalisation: straight line through two marked continuum points,
          F' = F - line + 1 between them.
        - The SNR of a stored SNR region is not used (and asked again) when
          it covers 40% or more of the window, i.e. when it comes from a run
          with a global SNR.
        - The standard deviation is the sample one (IDL STDDEV, N - 1).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from broadening.convol import _convol
    from gbat.grid import _find_model, _load_model
    from gbat.input import _ask, _question
    from gbat.io import _fa, _fg, _fi
    from gbat.variables import MICRO_REF


    # #################### REFERENCE MODEL ####################
    title = line
    ref = None
    if fw_ref:
        file = _find_model(cfg, fw_ref)
        if file is None:
            print(f'WARNING ... reference model {fw_ref} not found in {cfg.casa_red}')
        else:
            model = _load_model(file)
            if (line, MICRO_REF) in model['lines']:
                x, y = model['lines'][(line, MICRO_REF)]
                ref = _convol(x, y, cfg.res0, resol=cfg.resol, vsini=vsini, zeta_t=zeta_t, vgauss=vgauss)
            title = line + ' (Ref. model: ' + fw_ref + ')'
    plot = _CleanPlot(title, xl, ref)


    # #################### RENORMALISATION ####################
    ss = (wo0 >= xl[0]) & (wo0 <= xl[1])
    if not nocut or interactive:
        plot.draw(wo0, fo0)
        a = _ask('DO YOU WANT TO RENORMALIZE A REGION OF THE LINE? (0/1)')
        if a == 1:
            _question('MARK TWO SEPARATED POINTS IN THE CONTINUUM TO NORMALIZE (LINEAR CASE) ')
            (x1, y1), (x2, y2) = sorted(plot.marks(2))
            tt1 = (wo0 >= x1) & (wo0 <= x2)
            coef = np.polyfit([x1, x2], [y1, y2], 1)
            fo0[tt1] = fo0[tt1] - np.polyval(coef, wo0[tt1]) + 1.
            plot.draw(wo0, fo0)


    # #################### SNR ALREADY ESTABLISHED? ####################
    def snr_box(mm, sigma, cont):
        plot.box = (wo0[mm].min(), wo0[mm].max(), cont, sigma)

    nsigma = max(np.count_nonzero(sig0[ss] > 0), 1) # IDL n_elements(where(...))
    sigma = cont_flux = None
    if snr:
        sigma, cont_flux = 1. / snr, 1.
        mm = ss.copy()
        sig0[mm] = 1.
        snr_box(mm, sigma, cont_flux)
        print(_fa(line, 11) + _fa(': SNR = ', 8) + _fi(cont_flux / sigma, 4))
        if not nocut or interactive:
            plot.draw(wo0, fo0)
            if _ask('SNR ALREADY COMPUTED. DO YOU WANT TO COMPUTE IT AGAIN? (0/1) '):
                sig0[ss] = 0.
    elif sig0[ss].sum() > 0. and nsigma < 0.4 * np.count_nonzero(ss):
        mm = ss & (sig0 > 0.)
        sigma, cont_flux = np.std(fo0[mm], ddof=1), np.mean(fo0[mm])
        snr_box(mm, sigma, cont_flux)
        print(line + ': SNR = ' + _fg(cont_flux / sigma))
        if not nocut or interactive:
            plot.draw(wo0, fo0)
            if _ask('SNR ALREADY COMPUTED. DO YOU WANT TO COMPUTE IT AGAIN? (0/1) '):
                sig0[ss] = 0.
    else:
        sig0[ss] = 0.


    # #################### SNR OF THE CONTINUUM ####################
    if sig0[ss].sum() <= 0.:
        plot.box = None
        plot.draw(wo0, fo0, title=line)
        _question('MARK THE LIMITS OF THE REGION TO COMPUTE THE SNR  ')
        (x1, _), (x2, _) = sorted(plot.marks(2))
        sig0[(wo0 >= x1) & (wo0 <= x2)] = 1.
        mm = ss & (sig0 > 0.)
        if np.count_nonzero(mm) < 2:
            raise ValueError(f'{line}: the SNR region must contain at least 2 points inside the window of the line.')
        sigma, cont_flux = np.std(fo0[mm], ddof=1), np.mean(fo0[mm])
        snr_box(mm, sigma, cont_flux)
        print(line + ': SNR = ' + _fg(cont_flux / sigma))


    # #################### LIMITS OF THE LINE ####################
    if interactive and pe0[ss].sum() > 0.: # Stored limits: shown in red, asked whether to mark them again
        plot.draw(wo0, fo0, pe0, title=line)
        if _ask('LIMITS OF THE LINE ALREADY MARKED. DO YOU WANT TO MARK THEM AGAIN? (0/1) '):
            pe0[ss] = 0.
    if pe0[ss].sum() == 0.:
        plot.draw(wo0, fo0, title=line)
        _question('MARK THE LIMITS OF THE LINE            ')
        (x1, _), (x2, _) = sorted(plot.marks(2))
        pe0[(wo0 >= x1) & (wo0 <= x2)] = 1.
        if interactive:
            plot.draw(wo0, fo0, pe0, title=line) # Show the new limits


    # #################### CLEAN / RESTORE ####################
    a = 0 if nocut else 1
    while a in (1, 2):
        plot.draw(wo0, fo0, pe0, title=line)
        a = _ask('DO YOU WANT TO CLEAN (1) OR RESTORE (2) A REGION OF THE LINE? 0-no      ', (0, 1, 2))
        if a in (1, 2):
            _question('MARK THE LIMITS OF THE REGION YOU WANT TO CLEAN/RESTORE     ')
            (x1, _), (x2, _) = sorted(plot.marks(2))
            pe0[(wo0 >= x1) & (wo0 <= x2)] = 0. if a == 1 else 1.


    # #################### OUTPUT ####################
    return float(sigma)




# ====================================================================
# ============================== CHI2 ==============================
# ====================================================================

# ==============================
def _chi_one(
             dat_line : np.void,
             count    : int,
             wo0      : np.ndarray[float],
             fo0      : np.ndarray[float],
             pe0      : np.ndarray[float],
             sig0     : np.ndarray[float],
             model    : dict,
             micro    : str,
             vsini    : float,
             cfg      : Config,
             dat      : dict,
             zeta_t   : float = 0.,
             vgauss   : float = 0.,
             fw_ref   : str | None = None,
             snr      : float | None = None,
             nocut    : bool = False,
             interactive : bool = False
    ) -> tuple[float, int]:
    """
    chi2 of one line of one model (translation of gbat_chi_one.pro).

    ----------
    Parameters
    ----------
        dat_line : np.void
            Record of the line in dat_line (line, sigma, np), updated in
            place.

        count : int
            Index of the model in sol. The line is pre-processed only for
            the first model (count = 0).

        wo0, fo0, pe0, sig0 : np.ndarray[float]
            Observed spectrum and masks (see _clean_line).

        model : dict
            FASTWIND model (see iacob_gbat_grid._load_model).

        micro : str
            Microturbulence of the profile, e.g. '005'.

        vsini : float
            Projected rotational velocity [km/s].

        cfg : Config
            Configuration of the run (resol and res0 may be changed by
            _extra_info_obs).

        dat : dict
            Lines of gbat_linesHHe.dat (see _read_lines_hhe).

        zeta_t, vgauss, fw_ref, snr, nocut, interactive
            See _clean_line.

    -------
    Returns
    -------
        chi : float
            sum((F_obs - F_model)^2) over the points of the line (pe0 > 0)
            inside the window lam -/+ 20 A (35 A for the Balmer lines), or
            CHI_ERROR (999) if the line is not in the model.

        err_mes : int
            1 if the line is not in the model, 0 otherwise.

    -----
    Notes
    -----
        - The FASTWIND profile is convolved with the instrumental profile,
          vsini and the macroturbulence, and interpolated (with linear
          extrapolation, IDL INTERPOL) to the observed wavelengths.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from broadening.numerics import _interpol
    from broadening.convol import _convol
    from gbat.input import _extra_info_obs
    from gbat.variables import BALMER_LINES, CHI_ERROR, X1_BALMER, X1_DEFAULT


    # #################### WINDOW OF THE LINE ####################
    line = str(dat_line['line'])
    x0 = float(dat['lam'][np.nonzero(dat['fwlab'] == line)[0][0]])
    x1 = X1_BALMER if line in BALMER_LINES else X1_DEFAULT # [SS] Include also the IR H lines
    xl = x0 + x1 * np.array([-1., 1.])
    _extra_info_obs(cfg, x0)


    # #################### PRE-PROCESSING (ONLY ONCE) ####################
    if count == 0:
        dat_line['sigma'] = _clean_line(wo0, fo0, pe0, sig0, xl, line, vsini, cfg, zeta_t=zeta_t,
                                        vgauss=vgauss, snr=snr, nocut=nocut, fw_ref=fw_ref,
                                        interactive=interactive)


    # #################### POINTS OF THE LINE ####################
    tt = (wo0 >= xl[0]) & (wo0 <= xl[1]) & (pe0 > 0)
    wo, fo = wo0[tt], fo0[tt]
    dat_line['np'] = wo.size


    # #################### CHI2 ####################
    profile = model['lines'].get((line, micro))
    if profile is None or wo.size == 0 or not (profile[0][0] > 0 and profile[1][0] > 0):
        return CHI_ERROR, 1
    xm0, ym0 = _convol(profile[0], profile[1], cfg.res0, resol=cfg.resol, vsini=vsini, zeta_t=zeta_t, vgauss=vgauss)
    ym = _interpol(ym0, xm0, wo)
    chi = float(np.sum((fo - ym)**2))


    # #################### OUTPUT ####################
    return chi, 0
