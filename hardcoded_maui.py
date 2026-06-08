
# MAUI GRIDS FROM FASTWIND MODELS DEFINED AS
# full_name, short_name, color, order, log(Teff) limits, log(Ls) limits
# e.g.: nlte_10.OBSg, OBSgs, g, 0, 4.399 4.544 4.544 4.399 4.399, 3.488 3.488 4.386 4.386 3.488
#
# The limits are cuts a posteriori in MAUI's user_defined_model_set_fastwind_....pro in the PRO_BASE_USER folder and are introduced manually in gen_gridlimits()
# Note: O9BGs_CNOSiMg has Teff > 16000K, logQs <= -12.5 and lgf <= 2.4
# Note: O9BSg_CNOSiMg has the upper limit of the logL (lgf) constraint <= 3.54 (1.85) also the upper limit of logQs <= -12.5
#
dic_maui_grids = {
'all':                                                                             ['Grids coverage', 'dodgerblue', 0, '4.543 4.290 4.290 4.146 4.146 4.543 4.543', '2.391 2.391 3.092 3.092 4.391 4.391 2.391'],
'nlte_10.1.6_SOLAR_expoclump_2019-10-24':                                          ['BSgs_CNOSiMg',            'b', 1, '4.190 4.477 4.477 4.190 4.190', '3.785 3.785 4.391 4.391 3.785'],
'nlte_10.1.6_bdwarfs_SOLAR_2020-01-29':                                            ['BDws_CNOSIMg',       'orange', 2, '4.290 4.543 4.543 4.290 4.290', '2.391 2.391 3.889 3.889 2.391'],
'nlte_10.4.7_OB.Sg_SOLAR_2021-01-23':                                              ['OBSgs_hot_NOSi',          'g', 3, '4.399 4.544 4.544 4.399 4.399', '3.488 3.488 4.386 4.386 3.488'],
'nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06' :                      ['BSgs_cool_NOSi',          'r', 4, '4.146 4.322 4.322 4.146 4.146', '3.092 3.092 4.391 4.391 3.092'],
'astar2013_SOLAR_2_LMC_4_grid_2019-10-24_2019-10-24' :                          ['ASgs_CNOMgSTiFe_Kurucz','purple', 5, '3.900 4.114 4.114 3.900 3.900', '3.142 3.142 4.292 4.292 3.142'],
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1_2021-05-05' :       ['BSg_CNOSiMg',      'DeepPink', 6, '4.146 4.477 4.477 4.146 4.146', '3.392 3.392 4.386 4.386 3.392'],
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19' :   ['O9BSg_CNOSiMg',   'turquoise', 7, '4.146 4.543 4.543 4.146 4.146', '3.540 3.540 4.394 4.394 3.540'],
'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21' : ['O9BGs_CNOSiMg',        'lime', 8, '4.204 4.543 4.543 4.204 4.204', '2.937 2.937 3.791 3.791 2.937']
}

ad_hoc_limits = {
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19': {'Teff_DW': 1.4, 'lgf_UP': 1.85, 'logQs_UP': -12.2},
'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21': {'Teff_DW': 1.6, 'lgf_UP': 2.40, 'logQs_UP': -12.5}
}

dic_maui_param={'Teff': 'Teff',
                'logg': 'logg',
                'lgf': 'lgf',
                'He': 'He',
                'Micro': 'Micro',
                'logQs': 'logQs',
                'beta': 'beta',
                'C': 'C',
                'N': 'N',
                'O': 'O',
                'Mg': 'Mg',
                'Si': 'Si',
                'S': 'S',
                'Fe': 'Fe',
                'Ti': 'Ti',
                'fcl': 'fcl',
                'vcl': 'vcl'
}

# taken from de Burgos et al. 2024, A&A, 687, A228
dic_maui_uncertainties = {
'Teff': 0.05, # in kK
'logg': 0.07,
'logQs': 0.12,
'beta': 0.41,
'He': 0.02,
'micro': 1.5,
'Si': 0.08
}

# Define the regions used with weight=1 that are evaluated in the chi2 calculation.
# They must match the exact use_defined_model_fitness.pro used in the analysis.
# If only a region within the window is used, then the rest of the window has weight 0.
# In this case, I add -0.1 and +0.1 to the window limits defined in '*_lines_for_chi2_*'
# Otherwise I specify the exact region with weight 1 within the window.
mask_maui_SR = [
    (6521.00, 6532.00), # Halpha, mask NII lines in the blue wing [??]
    (6575.75, 6585.29), #         mask CII lines in the red wing
    (4330.50, 4335.00), # Hgamma, mask NIII/SIII? lines in the blue wind
    (4344.42, 4355.83), #         mask OII lines in the red wing
    (4083.20, 4085.50), # Hdelta, mask metal lines
    (4086.29, 4090.53),
    (4092.27, 4093.99),
    (4096.50, 4098.00),
    (3953.09, 3955.05), # Hepsil, mask metal lines
    (3960.70, 3962.24),
    (3963.38, 3965.76),
    (3967.00, 3969.00),
    (3972.53, 3974.10),
    (4465.00, 4468.50), # HeI 4471, mask OII lines in the blue wind
    (4923.50, 4926.00), # HeI 4922, mask an OII line in the red wing
    (5013.20, 5014.50), # HeI 5015, mask a SII line
    (5017.50, 5019.50), #           mask some Ni I/II lines?
    (4544.00, 4546.00), # HeII 4541, mask AlIII lines
    (4478.87, 4480.20), # MgII 4481, mask the AlIII blend
    (4128.71, 4130.10), # SiIII 4130
    (4131.40, 4134.08),
    (4547.60, 4550.20), # SiIII 4552 (remove continuum)
    (4555.00, 4558.96),
    (4563.08, 4565.00), # SiIII 4567 (remove continuum)
    (4571.00, 4571.35),
    (4571.30, 4572.50), # SiIII 4575 (remove continuum)
    (4576.00, 4579.32),
    (4112.50, 4113.80), # SiIV 4116 (first val should match lmin)
    (4117.70, 4124.70), # (last val should match lmax)
]

mask_maui_FR = [
    (6521.00, 6532.00), # Halpha, mask NII line in the blue wing [??]
    (6575.75, 6585.29), #         mask CII lines in the red wing
    (4330.50, 4335.00), # Hgamma, mask NIII/SIII? lines in the blue wind
    (4344.42, 4355.83), #         mask OII lines in the red wing
    (4083.20, 4085.50), # Hdelta, mask metal lines
    (4086.29, 4090.53),
    (4092.27, 4093.99),
    (4096.50, 4098.00),
    (3953.09, 3955.05), # Hepsil, mask metal lines
    (3960.70, 3962.24),
    (3963.38, 3965.76),
    (3967.00, 3969.00),
    (3972.53, 3974.10),
    (4465.00, 4468.50), # HeI 4471, mask OII lines in the blue wind
    (4923.50, 4926.00), # HeI 4922, mask an OII line in the red wing
    (5013.20, 5014.50), # HeI 5015, mask a SII line
    (5017.50, 5019.50), #           mask some Ni I/II lines?
    (4544.00, 4546.00), # HeII 4541, mask AlIII lines
    (4478.87, 4480.20), # MgII 4481, mask the AlIII blend
    (4128.71, 4130.10), # SiIII 4130
    (4131.40, 4134.08),
    (4547.40, 4549.00), # SiIII 4552 (remove continuum)
    (4556.00, 4558.90),
    (4576.00, 4581.40), # SiIII 4567 + 4575 (remove continuum)
    (4112.50, 4113.80), # SiIV 4116 (first val should match lmin)
    (4117.70, 4126.10), # (last val should match lmax)
]