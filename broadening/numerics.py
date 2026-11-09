from __future__ import annotations

"""
pyiacob_broad_numerics — Numerical routines of IDL without a direct NumPy/SciPy equivalent.

Only three operations are needed here because NumPy/SciPy do not provide the
same mathematics:

    - interpol: IDL INTERPOL extrapolates linearly outside the tabulated range;
                numpy.interp clamps to the end values.

    - median_upper: IDL MEDIAN returns the upper of the two central values
                    for an even number of points; numpy.median returns their mean.

    - spline_tension: IDL SPLINE is a cubic spline under tension (Cline 1974);
                    SciPy has no spline under tension.

Everything else of the original IDL code is written directly with NumPy.
"""

# =======================================================================
# ============================== FUNCTIONS ==============================
# =======================================================================

# ==============================
def _interpol(
              v : np.ndarray[float],
              x : np.ndarray[float],
              u : np.ndarray[float]
    ) -> np.ndarray[float]:
    """
    Perform linear interpolation with linear extrapolation.

    This function reproduces the behaviour of the IDL INTERPOL
    routine. For each requested abscissa u, the value is obtained
    by linear interpolation between the two neighbouring tabulated
    points in (x, v).

    ----------
    Parameters
    ----------
        v : numpy.ndarray[float]
            One-dimensional array containing the tabulated values
            corresponding to the abscissae in x.

        x : numpy.ndarray[float]
            One-dimensional array containing the tabulated abscissae.
            Values must be monotonically increasing.

        u : numpy.ndarray[float]
            Abscissae at which the interpolated values are evaluated.
            Values inside the range of x are linearly interpolated,
            while values outside the range are linearly extrapolated.

    -------
    Returns
    -------
        interp : numpy.ndarray[float]
            Array containing the interpolated or extrapolated values
            evaluated at u.

    -----
    Notes
    -----
        - For each value u[i], the interpolation is performed over the
          interval [x[s], x[s + 1]] containing that value, using
            v[s] + (u - x[s]) * (v[s + 1] - v[s]) / (x[s + 1] - x[s])
          For u < x[0], the first interval [x[0], x[1]] is used for
          linear extrapolation. For u > x[-1], the last interval
          [x[-2], x[-1]] is used.

        - This differs from numpy.interp, which by default returns the
          endpoint values for abscissae outside the interpolation range
          rather than performing linear extrapolation.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### INTERPOLATION ####################
    # ########## ARRAYS ##########
    v, x, u = np.asarray(v), np.asarray(x), np.asarray(u)


    # ########## CLIP ##########
    s = np.clip(np.searchsorted(x, u, side='right') - 1, 0, v.size - 2) # Interval of each point: x[s] <= u < x[s+1], clipped to the end intervals


    # ########## COMPUTATION ##########
    interp = (u - x[s]) * (v[s + 1] - v[s]) / (x[s + 1] - x[s]) + v[s]


    # #################### OUTPUT ####################
    return interp




# ==============================
def _median_upper(
                  x : np.ndarray[float]
    ) -> float:
    """
    Return the IDL-style median of an array.

    The input array is flattened, sorted in ascending order, and the
    element at index N // 2 is returned, where N is the total
    number of elements. This reproduces the behaviour of the IDL
    MEDIAN function without the /EVEN keyword.

    For an odd number of elements, this is the conventional median.
    For an even number of elements, the upper of the two central values
    is returned rather than their arithmetic mean, as would be returned
    by numpy.median.

    ----------
    Parameters
    ----------
        x : numpy.ndarray[float]
            Input array. The array may have any shape; all elements are
            flattened before computing the median.

    -------
    Returns
    -------
        median : float
            The median value, defined as the element at position N // 2
            in the sorted, flattened input array.

    -----
    Notes
    -----
        - For N elements sorted as
            x[0] <= x[1] <= ... <= x[N-1]
        the returned value is
            x[N // 2]`.
        Therefore, for an even number of elements, the result is the upper
        central value:
            [1, 2, 3, 4] -> 3

        - In Python, numpy.median returns
            numpy.median([1, 2, 3, 4]) -> 2.5.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### MEDIAN ####################
    median = np.sort(np.ravel(x))[np.size(x) // 2]


    # #################### OUTPUT ####################
    return median




# ==============================
def _spline_tension(
                    x     : np.ndarray[float],
                    y     : np.ndarray[float],
                    t     : np.ndarray[float],
                    sigma : float = 1.0
    ) -> np.ndarray[float]:
    """
    Evaluate a cubic spline under tension at the specified abscissae.

    This function implements the cubic spline-under-tension algorithm used
    by the IDL SPLINE routine, following the formulation of Cline
    (1974, CACM, 17, 218) and the corresponding FITPACK CURV1/CURV2
    routines.

    The spline passes through all input knots (x, y) and provides a
    smooth interpolation between them. The tension parameter sigma
    controls the stiffness of the spline: small values produce a curve
    closer to a conventional cubic spline, whereas increasing sigma
    progressively suppresses curvature and makes the interpolation more
    nearly linear.

    ----------
    Parameters
    ----------
        x : numpy.ndarray[float]
            One-dimensional array containing the knot abscissae. Values must
            be strictly increasing.

        y : numpy.ndarray[float]
            One-dimensional array containing the values associated with the
            knots in x. x and y must contain the same number of
            elements.

        t : numpy.ndarray[float]
            Abscissae at which the spline is evaluated. Scalar values and
            array-like inputs are accepted.

        sigma : float, optional
            Tension parameter. Values smaller than 0.001 are replaced by
            0.001, following the behaviour of the IDL routine. The tension
            is normalised by the mean knot spacing before constructing the
            spline. For the IACOB-BROAD implementation, sigma=5 is used.

    -------
    Returns
    -------
        spline_tension : numpy.ndarray[float]
            Spline values evaluated at t. The returned array has the same
            shape as the input t after conversion to at least one
            dimension.

    -----
    Notes
    -----
        - The endpoint derivatives are obtained from quadratic interpolants
          through the first and last three knots. The spline coefficients are
          then determined by solving the resulting tridiagonal system using
          recursive forward elimination and back substitution, following the
          algorithm described by Cline (1974).

        - Evaluation within each knot interval uses the hyperbolic functions
          associated with the tension-spline formulation. Abscissae outside
          the knot range are evaluated using the first or last spline interval,
          respectively, providing linear/tension-spline extrapolation.

    ----------
    References
    ----------
        Cline, A. K. 1974, Communications of the ACM, 17, 218
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### SPLINE UNDER TENSION ####################
    # ########## INPUT ARRAYS & PARAMETERS ##########
    xx = np.asarray(x, dtype=np.float64)
    yy = np.asarray(y, dtype=np.float64)
    tt = np.atleast_1d(np.asarray(t, dtype=np.float64))

    n = min(xx.size, yy.size) # Number of knots and index of the last knot.
    nm1 = n - 1

    sigma = max(sigma, 0.001) # IDL SPLINE imposes a minimum tension of 0.001.


    # ########## ENDPOINT SLOPES ##########
    """
    - Estimate the spline derivatives at the first and last knots.
    - For three or more knots, the endpoint slopes are obtained from
    - quadratic interpolants through the first and last three points.
    """
    if n == 2:
        slpp1 = slppn = 0.0
    else:
        # Slope at the first knot from a quadratic through x[0:3].
        d1, d2, d12 = xx[1] - xx[0], xx[2] - xx[1], xx[2] - xx[0]
        slpp1 = (-(d12 + d1) / d12 / d1 * yy[0]
                 + d12 / d1 / d2 * yy[1]
                 - d1 / d12 / d2 * yy[2])

        # Slope at the last knot from a quadratic through x[-3:].
        dn = xx[nm1] - xx[nm1 - 1]
        dnm1 = xx[nm1 - 1] - xx[nm1 - 2]
        dnn = xx[nm1] - xx[nm1 - 2]
        slppn = (dn / dnn / dnm1 * yy[nm1 - 2]
                 - dnn / dn / dnm1 * yy[nm1 - 1]
                 + (dnn + dn) / dnn / dn * yy[nm1])


    # ########## TRIDIAGONAL SYSTEM ##########
    """
    - Solve the tridiagonal system that determines the spline
      coefficients. The first n elements of yp contain the solution,
      while the remaining n elements store the elimination factors
      required during forward elimination and back substitution.
    """

    # ##### PARAMS #####
    yp = np.zeros(2 * n)
    sigmap = sigma * nm1 / (xx[nm1] - xx[0]) # Normalise the tension parameter by the mean knot spacing.

    # ##### FORWARD ELIMINATION: FIRST KNOT #####
    delx1 = xx[1] - xx[0] # Initialise the tridiagonal system using the first knot interval.
    dx1 = (yy[1] - yy[0]) / delx1

    dels = sigmap * delx1 # Dimensionless tension across the first interval.

    exps = np.exp(dels) # Hyperbolic functions required by the spline-under-tension formulation, evaluated through exponentials.
    sinhs = 0.5 * (exps - 1.0 / exps)
    sinhin = 1.0 / (delx1 * sinhs)

    diag1 = sinhin * (dels * 0.5 * (exps + 1.0 / exps) - sinhs) # Diagonal and inverse diagonal terms of the tridiagonal system.
    diagin = 1.0 / diag1

    yp[0] = diagin * (dx1 - slpp1)  # Right-hand side and elimination factor for the first knot.
    spdiag = sinhin * (sinhs - dels)
    yp[n] = diagin * spdiag

    dx2 = dx1 # Store the slope of the previous knot interval.

    # ##### FORWARD ELIMINATION: INTERIOR KNOTS #####
    """
    - Process each interior knot and recursively eliminate the
      sub-diagonal elements of the tridiagonal system.
    """
    for i in range(1, nm1):
        delx2 = xx[i + 1] - xx[i] # Width and data slope of the current knot interval.
        dx2 = (yy[i + 1] - yy[i]) / delx2

        dels = sigmap * delx2 # Dimensionless tension across the current interval.

        exps = np.exp(dels) # Evaluate the hyperbolic terms for the current interval.
        sinhs = 0.5 * (exps - 1.0 / exps)
        sinhin = 1.0 / (delx2 * sinhs)

        diag2 = sinhin * (dels * 0.5 * (exps + 1.0 / exps) - sinhs) # Diagonal coefficient of the current row.

        diagin = 1.0 / (diag1 + diag2 - spdiag * yp[n + i - 1]) # Eliminate the contribution from the previous row.

        yp[i] = diagin * (dx2 - dx1 - spdiag * yp[i - 1]) # Store the current solution term.

        spdiag = sinhin * (sinhs - dels) # Compute and store the elimination factor for the next row.
        yp[i + n] = diagin * spdiag

        dx1, diag1 = dx2, diag2 # Advance the interval coefficients for the next iteration.

    # ##### FORWARD ELIMINATION: LAST KNOT #####
    """
    - Apply the final boundary condition using the slope at the
      last knot and complete the forward elimination.
      """
    diagin = 1.0 / (diag1 - spdiag * yp[n + nm1 - 1])
    yp[nm1] = diagin * (slppn - dx2 - spdiag * yp[nm1 - 1])

    # ##### BACK SUBSTITUTION #####
    """
    - Recover the spline coefficients from the last knot towards
      the first using the elimination factors stored in yp[n:].
    """
    for i in range(n - 2, -1, -1):
        yp[i] = yp[i] - yp[i + n] * yp[i + 1]


    # ########## SPLINE EVALUATION (CURV2) ##########
    """
     - Determine the knot interval containing each evaluation point.
       k is the index of the right-hand knot:
        xx[k-1] <= tt <= xx[k]

     - Values outside the knot range use the first or last interval.
    """
    k = np.clip(np.searchsorted(xx[:n], tt, side='right'), 1, nm1)

    del1 = tt - xx[k - 1] # Distances from each evaluation point to the left and right knots, together with the total width of the corresponding interval.
    del2 = xx[k] - tt
    dels = xx[k] - xx[k - 1]

    exps1 = np.exp(sigmap * del1) # Hyperbolic sine terms associated with the distances from the evaluation point to the left and right knots.
    sinhd1 = 0.5 * (exps1 - 1.0 / exps1)

    exps = np.exp(sigmap * del2)
    sinhd2 = 0.5 * (exps - 1.0 / exps)

    exps = exps1 * exps # Hyperbolic sine across the complete knot interval.
    sinhs = 0.5 * (exps - 1.0 / exps)

    spline_tension =  (yp[k] * sinhd1 + yp[k - 1] * sinhd2) / sinhs + ((yy[k] - yp[k]) * del1 + (yy[k - 1] - yp[k - 1]) * del2) / dels # Evaluate the spline-under-tension using the coefficients obtained from the tridiagonal system.


    # #################### OUTPUT ####################
    return spline_tension

