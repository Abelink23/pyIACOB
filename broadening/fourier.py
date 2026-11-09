from __future__ import annotations
from typing import NamedTuple

# ===============================================================================
# ============================== FOURIER TRANSFORM ==============================
# ===============================================================================

# ==============================
class FT(NamedTuple):
    """
    Fourier transform of a line profile, in vsini units.
    """
    vsini: np.ndarray # vsini scale, 1/(sigma*factor) [km/s] (Inf at zero frequency)
    log_amp: np.ndarray # log10 of the amplitude normalised to the zero frequency
    factor: float # vsini = 1/(sigma*factor), sigma the frequency [1/A]
    vsini_nyq: float # minimum measurable vsini (Nyquist frequency) [km/s]
    ew: float # equivalent width from the zero-frequency amplitude [A]




# ==============================
def vfft(
        wx        : np.ndarray[float],
        fx        : np.ndarray[float],
        xc        : float,
        sigma_FT0 : float | None = None,
        nyq       : bool = False
    ) -> FT:
    """
    Compute the Fourier transform of a normalised spectral-line profile.

    The input line profile is first apodised at both wavelength edges
    using a cosine taper and then extended with continuum points on both
    sides. The Fourier transform of the line depth, 1 - F, is then
    computed and normalised to its zero-frequency amplitude.

    The Fourier frequencies are converted into an equivalent projected
    rotational velocity, v sin(i), using the position of the first
    zero of the rotational broadening function.

    ----------
    Parameters
    ----------
        wx : np.ndarray[float]
            Wavelength array of the spectral-line profile, in Angstroms.
            The wavelength grid must be equidistant.

        fx : np.ndarray[float]
            Normalised flux corresponding to wx. The continuum level is
            assumed to be unity.

        xc : float
            Central wavelength of the spectral line, in Angstroms.

        sigma_FT0 : float | None, optional
            Dimensionless position of the first zero of the Fourier
            transform of the rotational broadening profile. If None,
            the module default SIGMA_FT0 is used.

        nyq : bool, optional
            If True, print the minimum projected rotational velocity
            allowed by the Nyquist frequency of the wavelength grid.
            Defaults to False.

    -------
    Returns
    -------
        FT
            Fourier-transform results containing:

            - vsini : np.ndarray[float]
                Projected rotational velocity corresponding to each Fourier
                frequency, in km/s.

            - log_amp : np.ndarray[float]
                Base-10 logarithm of the normalised Fourier-transform
                amplitude.

            - factor : float
                Conversion factor between Fourier frequency and projected
                rotational velocity.

            - vsini_nyq : float
                Minimum projected rotational velocity corresponding to the
                Nyquist frequency, in km/s.

            - ew : float
                Equivalent width of the line, in Angstroms.

    -----
    Notes
    -----
        - The continuum level is fixed to unity.

        - Before computing the Fourier transform, FT_EDGE_PERC percent of
          the points at each edge of the profile are apodised with a
          half-cosine window. The taper is applied to the line depth while
          preserving the continuum level:
            F' = (F - 1) * window + 1.

        - The profile is subsequently extended on both sides with
          FT_EXTENSION * N continuum points. The wavelength spacing of
          the input grid is preserved during this extension.

        - The Fourier transform is calculated from the line depth,
            F - 1,
        using the NumPy FFT convention and normalised by the total number
        of points. The zero-frequency amplitude is used to determine the
        equivalent width and to normalise the Fourier-transform amplitude.

        - The conversion from Fourier frequency sigma to projected
          rotational velocity is
            v sin(i) = 1 / (sigma * factor)
        where
            factor = xc / (c * sigma_FT0),
        with c given by C_KMS in km/s.

        - The first Fourier-transform point corresponds to zero frequency.
          Its associated v sin(i) is therefore formally infinite.

        - The Nyquist frequency is 1 / (2 * delta_lambda) and is converted
          into the corresponding minimum measurable v sin(i).

        - The equivalent width is obtained from the zero-frequency Fourier
          amplitude and is expressed in the same wavelength units as wx.
          Since wx is given in Angstroms, the returned equivalent width is
          in Angstroms.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .variables import C_KMS, FT_EDGE_PERC, FT_EXTENSION, SIGMA_FT0


    # #################### PARAMETERS ####################
    sigma_FT0 = sigma_FT0 or SIGMA_FT0 # Use the user-provided first-zero parameter or the module default.
    wf = np.asarray(wx) # Convert the input arrays to NumPy arrays.
    ff = np.asarray(fx)
    n = wf.size


    # #################### COSINE TAPER ####################
    # ########## NUMBER OF POINTS TO APODISE ##########
    nxw = int(FT_EDGE_PERC * n / 100. + 0.5) # Determine the number of points to apodise at each edge. IDL NINT is reproduced by adding 0.5 before truncation for positive values.
    if nxw < 1:
        raise ValueError(
            "vfft: Array dimensions must be greater than 0 "
            "(IDL error in WIN_COS)."
        )


    # ########## WINDOW FUNCTION ##########
    wi = 0.5 * (1. - np.cos(np.pi * np.arange(nxw) / nxw)) # Construct the rising half-cosine window.

    window = np.ones(n)
    window[:nxw] = wi
    window[n - nxw:] = wi[::-1]


    # ########## APPLY TAPER ##########
    ff = (ff - 1) * window + 1 # Apply the taper to the line depth while keeping the continuum level equal to unity.


    # #################### CONTINUUM EXTENSION ####################
    # ########## NUMBER OF POINTS TO EXTEND ##########
    plus = FT_EXTENSION * n # Extend the profile with continuum points on both sides.


    # ########## STEP SIZE ##########
    step = wf[1] - wf[0] # Preserve the original wavelength sampling.
    k = np.arange(plus)


    # ########## EXTEND PROFILE ##########
    wf = np.concatenate([wf[0] - (plus - k) * step, wf, wf[-1] + (k + 1) * step])
    ff = np.concatenate([np.ones(plus), ff, np.ones(plus)])


    # #################### FOURIER TRANSFORM ####################
    # ########## PARAMETERS ##########
    nx = wf.size # Total number of points after continuum extension.
    t = wf[3] - wf[2] # Wavelength sampling interval.
    ds = 1. / (t * nx) # Fourier-frequency sampling interval.


    # ########## COMPUTE FOURIER TRANSFORM ##########
    amp = np.abs(np.fft.fft(ff - 1)) / nx # Compute the FFT of the line depth. Division by nx reproduces the normalisation.


    # ########## RESULTS ##########
    # ##### EQUIVALENT WIDTH #####
    ew = amp[0] * nx * t # The zero-frequency amplitude corresponds to the integrated line depth, i.e. the equivalent width.

    # ##### NORMALISE AMPLITUDE #####
    amp = amp / amp[0] # Normalise the Fourier-transform amplitude to its zero-frequency value.

    # ##### NUMBER OF USEFUL POINTS #####
    nr = nx - 1 # Number of useful positive-frequency points.


    # #################### FREQUENCY -> VSINI CONVERSION ####################
    # ########## CONVERSION FACTOR ##########
    factor = xc / C_KMS / sigma_FT0 # Conversion factor between Fourier frequency and projected rotational velocity.


    # ########## FREQUENCY ARRAY ##########
    freq = np.concatenate([[0.], np.arange(1, nr) * ds])  # Construct the frequency array, including the zero-frequency component explicitly.


    # ########## VSINI ARRAY ##########
    with np.errstate(divide='ignore'): # Convert Fourier frequencies into projected rotational velocities.  The zero-frequency point corresponds formally to infinite vsini.
        vsini = 1. / (freq * factor)

    vsini_nyq = 1. / (0.5 / t * factor) # Convert the Nyquist frequency into the minimum measurable vsini.

    if nyq:
        print(f"{'max_sigma_obs':>13} = {vsini_nyq:.1f} km/s (minimum vsini, Nyquist)")


    # #################### OUTPUT ####################
    return FT(
              vsini,
              np.log10(amp[:nr]),
              factor,
              vsini_nyq,
              ew,
             )
