from __future__ import annotations

"""
gbat.variables — User configuration and constants of IACOB-GBAT.

This module replaces the user-editable parts of gbat_chi_input.pro and
gbat_extra_info_obs.pro. As in the IDL version, steps [1] to [5] below must
be checked (and modified if needed) before running IACOB-GBAT.
"""

# =====================================================================================
# ============================== USER CONFIGURATION (IDL: gbat_chi_input.pro) ==============================
# =====================================================================================

# #################### [1] PDF VIEWER ####################
"""
- Program used to open the final figure when 'ps=True' (IDL: ps_viewer).
"""
PS_VIEWER = 'evince'


# #################### [2] LOCATION OF THE FASTWIND GRIDS ####################
"""
- Folder holding the FASTWIND grids (IDL: casa). The grid of a metallicity is
  expected in
    >> CASA/iac_gbatgrid_HHe_<METAL><METAL_EXTRA>/
  with the list of models in
    >> CASA/iac_gbatgrid_HHe_<METAL><METAL_EXTRA>/iac_gbatgrid_HHe_<METAL><METAL_EXTRA>_allHe.lis

- CASA can be a local folder or the link of a Nextcloud public share (see
  iacob_gbat_remote). With a link, the files of the grid are downloaded the
  first time they are needed (only the models of each subgrid) and kept in
  GRID_CACHE. The share may hold the folder iac_gbatgrid_HHe_<METAL>/ or be
  that folder itself.

- Each model is an IDL save file '<sub-folder>/<model>.xdr' (variables m, mag,
  vsyn and l). A pickle file '<model>.idl' with the same variables is also
  accepted (see iacob_gbat_grid._load_model).

- In the share, each set of models (sub-folder of the list of models, e.g.
  iac_gbatgrid_HHe_Z100_He08_xdr) may be a .tgz file with its .xdr files:
  only the sets needed by each subgrid are downloaded and unpacked. Sets can
  also be chosen and downloaded in advance: iacob_gbat_remote.download_grid().

- To check the access to the share: iacob_gbat_remote.check_grid() and
  iacob_gbat_remote.list_grid().
"""

CASA = 'https://cloud.iac.es/index.php/s/4TptbGxFZE9eQGs'
#CASA = '/net/nas/proyectos/hots/AAA_WORKING_2030/grids/fastwind/user_grids/ssimon/' # Local folder (IDL default)

GRID_PASSWORD = ''                     # Password of the share, if it is protected
GRID_CACHE = '~/.iacob_gbat_grid_cache' # Local copy of the downloaded files of the grid
GRID_THREADS = 8                       # Parallel downloads of the models of a subgrid
GRID_TIMEOUT = 120.                    # Timeout of each download [s]
GRID_TGZ_EXT = ('.tgz', '.tar.gz')     # Extensions of the sets of models in the share (<set folder><ext>)
GRID_KEEP_TGZ = False                  # Keep the .tgz files in the cache after unpacking them


# #################### [3] METALLICITY ####################
METAL = 'Z100'      # Other available options 'Z050', 'Z020'
METAL_EXTRA = ''    # In case there is an extra label in the metal folder (e.g. '_v106')


# #################### [4] CHARACTERISTICS OF THE OBSERVATIONS ####################
"""
- resol      : Resolving power of the observations.
- res0       : Resolution in A/pix for the convolved FASTWIND profiles.
               ** Recommended: ~ 0.3-0.5 of the A/pix of the observed spectrum.
- extra_info : (0/1) Use the wavelength-dependent values of EXTRA_INFO_OBS.

- To specify different values of resol & res0 for different spectral ranges,
  configure EXTRA_INFO_OBS below and select extra_info=1.
"""
INSTRUMENTS = {
               #  inst       : (resol,  res0, extra_info)
               'vfts'        : ( 8000., 0.10, 1), # FOR VFTS analyses (see also EXTRA_INFO_OBS)
               'bloems'      : ( 8000., 0.10, 0), # FOR BLOeMS analyses
               'hermes85'    : (85000., 0.02, 0), # FOR HERMES analyses
               'feros48'     : (48000., 0.05, 0), # FOR FEROS analyses
               'fies67'      : (67000., 0.05, 0), # FOR FIES hi-res analyses
               'fies46'      : (46000., 0.05, 0), # FOR FIES med-res analyses
               'fies25'      : (25000., 0.08, 0), # FOR FIES low-res analyses
               'ges'         : (20000., 0.08, 0), # For GES analyses
               'weave-lr'    : ( 5000., 0.10, 0), # For WEAVE (LR) analyses
               'megara'      : (20000., 0.08, 0), # For MEGARA analyses
              }


# ########## DETAILED CHARACTERISATION (IDL: gbat_extra_info_obs.pro) ##########
"""
- Example of detailed characterisation of the observations: VFTS case.
- Each entry is (lambda_min, lambda_max, resol, res0) [A]. The ranges are
  checked in order, and a later matching range overrides an earlier one, as
  in the IDL code. If the line is outside all of them, the previous values
  are kept.
"""
EXTRA_INFO_OBS = [
                  (3900., 4510.,  7000., 0.10),
                  (4510., 6470.,  8000., 0.10),
                  (6470., 7000., 17000., 0.06),
                 ]


# #################### [5] LINES CONSIDERED BY IACOB-GBAT ####################
"""
- File with the lines of the analysis: 'gbat_linesHHe_<obs_run>.dat' in the
  working folder. Create your own by modifying gbat_linesHHe_benchmark.dat.
  As in IDL, it is copied to 'gbat_linesHHe.dat', the file read by the code.
"""
LINES_FILE = 'gbat_linesHHe_{obs_run}.dat'
LINES_FILE_ACTIVE = 'gbat_linesHHe.dat'




# ===================================================================================
# ============================== VARIABLES & CONSTANTS ==============================
# ===================================================================================

# #################### SPEED OF LIGHT ####################
"""
- Speed of light [km/s] used for the radial-velocity correction (2.99e5 in
  IDL). The convolutions use the exact value C_KMS_CONVOL of IACOB-BROAD.
"""
C_KMS = 2.99e5
C_KMS_CONVOL = 299792.458 # Exact value, used in the convolutions (gbat_convol.pro)


# #################### FOLDERS ####################
"""
- Sub-folders of the working folder (IDL: relative paths). Except for
  summary, tables_launcher, trash, all are organised in
  sub-folders of the observing run (obs_run).
"""
GBAT_CUT_DIR = 'spectra_cut'           # Pre-processed spectra (4 columns)
GBAT_BFM_DIR = 'spectra_bfm'           # Synthetic spectra of the best-fitting model
GBAT_FILES_DIR = 'xdr_files'           # Results of step 1 (chi_global)
GBAT_DIST_DIR = 'xdr_dist'             # chi2 distributions (chi_step2)
GBAT_RESULTS_DIR = 'xdr_results'       # Results for each parameter (plot_dist_one)
GBAT_FIGURES_DIR = 'output_figures'    # PDF summary figures
GBAT_TRASH_DIR = 'trash'               # Temporary files
GBAT_SUMMARY_DIR = 'summary'           # Summary tables of a list of stars
GBAT_LAUNCHER_DIR = 'tables_launcher'  # Input tables of the launcher


# #################### SAVE FILES ####################
"""
- Results are saved as pickle files with extension '.idl' (as IACOB-BROAD),
  together with a readable '.json' copy if JSON_COPIES is True. IDL '.xdr'
  files with the same name are read when the '.idl' file does not exist.
"""
SAVE_EXT = '.idl'
IDL_EXT = '.xdr'
JSON_COPIES = True


# #################### DEFAULT RANGES OF PARAMETERS (gbat_chi_global.pro) ####################
TEFF_RANGE = (19000., 61000.)  # [K]
GRAV_RANGE = (2.1, 4.5)        # [dex]
Q_RANGE = (11., 16.)           # -log Q
YHE_RANGE = (0.04, 0.31)       # N(He)/N(H)
BETA_RANGE = (0.7, 2.1)
MICRO_DEFAULT = ('001', '005', '009', '015', '020', '025', '030') # [km/s]


# #################### MICROTURBULENCES OF THE LAUNCHER (gbat_launcher.pro) ####################
MICRO_HIGH = ('001', '005', '009', '015', '020', '025', '030') # mean(micro range) >= 10 km/s
MICRO_LOW = ('001', '003', '005', '007', '009', '011', '015', '020', '025', '030')
LAUNCHER_SNR = 150. # Global SNR used by the launcher


# #################### LINE WINDOWS (gbat_chi_one.pro) ####################
X1_DEFAULT = 20. # Half-width of the window of the line [A]
X1_BALMER = 35.  # Half-width for the optical Balmer lines [A]
BALMER_LINES = ('HALPHA', 'HBETA', 'HGAMMA', 'HDELTA')
X1_PROF = 20.    # Half-width of the panels of the profiles figure [A] (gbat_lineshhe.pro)


# #################### CHI2 (gbat_chi_one.pro) ####################
CHI_ERROR = 999. # chi2 of a line not found in the model


# #################### MODEL USED FOR THE PRE-PROCESSING ####################
MICRO_REF = '005' # Microturbulence of the reference model (fw_ref)


# #################### ENVELOPE OF THE CHI2 DISTRIBUTIONS (gbat_plot_dist_one.pro) ####################
CURV = 80.   # Tension of the spline fit (high values mean linear fit)
NSTEP = 100  # Number of steps of the spline


# #################### STELLAR PARAMETERS (gbat_rlm.pro) ####################
MSUN = 1.989e33   # [g]
RSUN = 6.9598e10  # [cm]
G_CGS = 6.6742e-8 # [cgs]
TSUN = 5780.      # [K]
LSP_CONST = 10.61 # log L/Lsun = 4 log Teff - log g - 10.61
