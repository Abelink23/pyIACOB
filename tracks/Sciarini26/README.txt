The tracks are subdivided into 3 folders, corresponding to the sections 3.1, 3.2, 3.3 of the paper.

In each folder, the tracks are further subdivided according to the subsections of the paper. Finally, tracks may be further divided according to the mass of the computed model.
_____________________________________

Nomenclature

The name of the files gives informations on the model.
the mass stands after the M (M020: 20 Msol)
the metallicity stands after the Z (z014: Z=0.014; z002: Z=0.002)
the rotation rate stands after the V (V0: no rotation; V4: V/Vc=0.40; Vs: star initialized at synchronization)
for binaries: the mass of the companion stands after the B, the period after the P (B20P3: M2 = 20 Msol, P = 3 days)
angular momentum transport: hydro if not specified, "_magn" for magnetic models
Note: in principle, name of single stars models ends just after rotation (or after "magn" for magnetic models), but in case of ambiguity, we add "no" for "not a binary" (e.g., 
Grid I: stars initiated at synchronization, for singles we give the initial rotation they would have had if in a binary of the corresponding period -> period included in name 
of star even if tides not applied). If no "no" -> binary, tides accounted for
_____________________________________

Evolutionary data (single or binaries)

For each model the file "name_star.dat" is provided, usual genec models evolutionary data.

The different output quantities are described below:

col.   1: model number
col.   2: age [yr]
col.   3: mass [Msol]
col.   4: log L [Lsol]
col.   5: log Teff [K]
col.   6: 1H surface [mass fraction]
col.   7: 4He         "
col.   8: 3He         "
col.   9: 12C         "
col.  10: 13C         "
col.  11: 14N         "
col.  12: 16O         "
col.  13: 17O         "
col.  14: 18O         "
col.  15: 20Ne        "
col.  16: 22Ne        "
col.  17: convective core mass Mcc [M/Mtot]
col.  18: log Teff [K] not corrected for the thickness of the wind (WR)
col.  19: log Mdot [Msol yr-1]
col.  20: log ρc [g cm-3]
col.  21: log Tc [K]
col.  22: 1H centre [mass fraction]
col.  23: 4He         "
col.  24: 3He         "
col.  25: 12C         "
col.  26: 13C         "
col.  27: 14N         "
col.  28: 16O         "
col.  29: 17O         "
col.  30: 18O         "
col.  31: 20Ne        "
col.  32: 22Ne        "
col.  33: 7Be         "
col.  34: 8B          "
col.  35: neutrino flux from 7Be(e−, ν)7Li
col.  36: neutrino flux from 8B(e+, ν)8Be
col.  37: neutrino flux from 7Be(e−, ν)7Li in SNU units
col.  38: neutrino flux from 8B(e+, ν)8Be in SNU units
col.  39: Ωsurf/Ωcrit
col.  40: Ωsurf [rad s−1]
col.  41: Ωcenter [rad s−1]
col.  42: Rpol/Req polar to equatorial radius ratio
col.  43: 26Al surface [mass fraction]
col.  44: 26Al centre [mass fraction]
col.  45: mass-loss correction factor due to rotation Mdot(Ω) / Mdot(O)
col.  46: layer with max of CNO iCNO,max
col.  47: mass of the iCNO,max layer [Msol]
col.  48: sum of CNO elements in iCNO,max layer
col.  49: specific angular momentum within 3 Msol [1016 cm2 s−1]
col.  50: specific angular momentum within 5 Msol [1016 cm2 s−1]
col.  51: Vcrit,1 first root (before convergence) [km s−1]
col.  52: Vcrit,2 second root (before convergence) [km s−1]
col.  53: Veq equatorial velocity (before convergence) [km s−1]
col.  54: Ωsurf/Ωcrit ΤΩ limit (before convergence)
col.  55: ΤEdd Eddington factor (before convergence)
col.  56: Vcrit,1 first root (after convergence) [km s−1]
col.  57: Vcrit,2 second root (after convergence) [km s−1]
col.  58: Veq equatorial velocity (after convergence) [km s−1]
col.  59: Ωsurf/Ωcrit ΤΩ limit (after convergence)
col.  60: ΤEdd Eddington factor (after convergence)
col.  61: ΔMneed mass needed to be lost at equator (critical velocity) [Msol]
col.  62: Mdot,need mass-loss rate needed to lose ΔMneed [Msol yr-1]
col.  63: ΔLrad+aniso+mech angular momentum lost per timestep (wind + mechanical) [1053 g cm2 s−1]
col.  64: I total moment of inertia [1057 g cm2]
col.  65: Ltot total angular momentum [1053 g cm2 s−1]
col.  66: Erot rotational kinetic energy [E51]
col.  67: Epot gravitational potential energy [E51]
col.  68: Egaz thermal energy of the gas [E51]
col.  69: Erad radiative energy [E51]
col.  70-109: limits of convective zones (layer boundaries)
col. 110: Ltot total angular momentum (star + envelope) [1053 g cm2 s−1]
col. 111: percentage of main sequence phase
col. 112: True if star is O or B
col. 113: True if star is Red Supergiant
col. 114: True if star is WolF Rayet
col. 115: apsidal motion constant k2
_____________________________________

Orbital evolution (binaries)

For binaries, an other file is provided: "name_star.period_evol.dat", containing the orbital evolution.

The different output quantities are described below:

col.   1: model number
col.   2: age [yr]
col.   3: radius R [Rsol]
col.   4: Ωsurf [rad s−1]
col.   5: Ωorb [rad s−1]
col.   6: orbital period Porb [days]
col.   7: semi-major axis a [cm]
col.   8: eccentricity e      
col.   9: dynamical tides torque Jdot,dyn [g cm2 s−2]
col.  10: equilibrium tides torque Jdot,eq [g cm2 s−2]
col.  11: dynamical tides contribution to eccentricity variation edot,dyn [s−1]
col.  12: equilibrium tides contribution to eccentricity variation edot,eq [s−1]
_____________________________________
