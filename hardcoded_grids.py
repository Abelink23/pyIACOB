
# MAUI GRIDS FROM FASTWIND MODELS DEFINED AS
# full_name, short_name, color, order, log(Teff) limits, log(Ls) limits
# e.g.: nlte_10.OBSg, OBSgs, g, 0, 4.399 4.544 4.544 4.399 4.399, 3.488 3.488 4.386 4.386 3.488
#
# The limits are cuts a posteriori in MAUI's user_defined_model_set_fastwind_....pro in the PRO_BASE_USER folder and are introduced manually in gen_gridlimits()
# Note: O9BGs_CNOSiMg has Teff > 16000K, logQs <= -12.5 and lgf <= 2.4
# Note: O9BSg_CNOSiMg has the upper limit of the logL (lgf) constraint <= 3.54 (1.85) also the upper limit of logQs <= -12.5
#
grids_dic = {
'all':                                                                              ['Grids coverage', 'dodgerblue', 0, '4.543 4.290 4.290 4.146 4.146 4.543 4.543', '2.391 2.391 3.092 3.092 4.391 4.391 2.391'],
'nlte_10.1.6_SOLAR_expoclump_2019-10-24':                                           ['BSgs_CNOSiMg',            'b', 1, '4.190 4.477 4.477 4.190 4.190', '3.785 3.785 4.391 4.391 3.785'],
'nlte_10.1.6_bdwarfs_SOLAR_2020-01-29':                                             ['BDws_CNOSIMg',       'orange', 2, '4.290 4.543 4.543 4.290 4.290', '2.391 2.391 3.889 3.889 2.391'],
'nlte_10.4.7_OB.Sg_SOLAR_2021-01-23':                                               ['OBSgs_hot_NOSi',          'g', 3, '4.399 4.544 4.544 4.399 4.399', '3.488 3.488 4.386 4.386 3.488'],
'nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06' :                       ['BSgs_cool_NOSi',          'r', 4, '4.146 4.322 4.322 4.146 4.146', '3.092 3.092 4.391 4.391 3.092'],
'astar2013_SOLAR_2_LMC_4_grid_2019-10-24_2019-10-24' :                          ['ASgs_CNOMgSTiFe_Kurucz', 'purple', 5, '3.900 4.114 4.114 3.900 3.900', '3.142 3.142 4.292 4.292 3.142'],
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1_2021-05-05' :        ['BSg_CNOSiMg',      'DeepPink', 6, '4.146 4.477 4.477 4.146 4.146', '3.392 3.392 4.386 4.386 3.392'],
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19' :    ['O9BSg_CNOSiMg',   'turquoise', 7, '4.146 4.543 4.543 4.146 4.146', '3.540 3.540 4.394 4.394 3.540'],
'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21' :  ['O9BGs_CNOSiMg',        'lime', 8, '4.204 4.543 4.543 4.204 4.204', '2.937 2.937 3.791 3.791 2.937']
}

ad_hoc_limits = {
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19': {'Teff_DW': 1.4, 'lgf_UP': 1.85, 'logQs_UP': -12.2},
'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21': {'Teff_DW': 1.6, 'lgf_UP': 2.40, 'logQs_UP': -12.5}
}