
from __future__ import annotations

# ===================================================================
# ============================== LINES ==============================
# ===================================================================

# ==============================
from pathlib import Path
LINE_LIST = Path(__file__).resolve().parent / "lines.dat"




# ===================================================================================
# ============================== VARIABLES & CONSTANTS ==============================
# ===================================================================================


# #################### SPEED OF LIGHT ####################
"""
- Speed of light in vacuum [km/s].

- IACOB-BROAD uses 2.99e5 everywhere except inside convol, which uses the
  exact value; both are kept as in the IDL code.


SERGIO REVISAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
C_KMS = 2.99e5
C_KMS_CONVOL = 299792.458


# #################### GOF COMPUTATION ####################
# ########## MAXIMUM VSINI AND THETA_RT ##########
"""
- The GOF computation for vsini will run from
    >> vmin (one pixel of the line in velocity units, unless given by the user) to,
    >> vmax = FACT_VS*vsini(zero FT)

- The GOF computation for vmac will run from
    >> vmin (the same as for vsini) to,
    >> vmax = FACT_VRT*vsini(zero FT),

- NOTE: If a HeI is analyzed and an intrinsic
  delta-profile is assumed, a large fact_vrt value
  may be needed
"""
FACT_VS = 2.0 # factor to scale the maximum vsini of the GOF computation
FACT_VRT = 2.5 # factor to scale the maximum vmac of the GOF computation
FACT_VRT_HEI = 4.5 # factor to scale the maximum vmac in the GOF computation for HeI lines


# ########## GRID STEP ##########
"""
- The GOF computation for vsini and vmac will be done in a grid with a step of
    >> STEP_V*vft,
  where vft = vsini(zero FT) for vsini and vmac.
"""
STEP_V = 0.04 # Grid step in vsini and vmac for the GOF computation


# ########## EQUIVALENT WIDTH (EW) GRID ##########
"""
- The equivalent width (EW) grid will be defined as:
    >> EW0*(1 + k*EW_STEP),
  where k ranges from -n_down to n_up (N_EW) for
  normal lines and (N_EW_CUT) for clipped lines.
"""
EW_STEP = 0.05 # Step in EW (relative to the EW provided by FT(sigma=0))
N_EW = (0, 0) # Number of steps with (EW < EW0, EW > EW0) for normal lines
N_EW_CUT = (0, 3) # Number of steps with (EW < EW0, EW > EW0) for clipped lines


# ########## INTRINSIC PROFILE ##########
# ##### DELTA LIKE PROFILE #####
R_DELTA = 200000. # Resolution of the intrinsic delta-like profile (R = lambda/delta_lambda)

# ##### USER-DEFINED PROFILE #####
"""
SERGIO REVISAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
RV_RANGE = 10. # radial-velocity search [km/s] shifts from -RV_RANGE to +RV_RANGE
RV_STEP = 0.1 # radial-velocity step [km/s] for the user-defined intrinsic profile


# ########## FOURIER TRANSFORM (vfft) ##########
SIGMA_FT0 = 0.660 # First zero of the rotational FT (lambda*vsini/c units)
FT_EDGE_PERC = 10 # % of the line points tapered at each edge
FT_EXTENSION = 40 # Continuum extension factor (points = factor * n)


# ########## ROTATIONAL PROFILE ##########
BETA = 1.5 # Limb darkening beta = eps/(1-eps) (eps = 0.6)