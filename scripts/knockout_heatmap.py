#!/usr/bin/python
import matplotlib.pyplot as plt
import logging
import numpy as np
from itertools import combinations
from scipy.io import savemat

from src.optknock import OptKnock
from src.models import *
from src.analysis_toolbox import *
from src.html_writer import HtmlWriter

def main():
    analyze('wild_type', None, None)
    analyze('rubisco', 'RBC', 'RBC,PRK')
    analyze('deoxyribose', 'DXS', 'DXS')
    analyze('MOG', 'MCS', 'MCS,MCL')

def analyze(title, target_reaction, knockins):
    carbon_uptake_rate = 50 # mmol C / (gDW*h)
    max_number_of_knockouts = 2

    # the following gene pairs are combined into one "knockout" since they have no intermediate flux options:
    # G6PDH2r,PGL - PP-oxidative shunt
    # GAPD,PGK - 
    # ICL,MALS - glyoxlyate shunt
    # ALCD2x,ACALD - ethanol secretion
    single_ko_list = ['GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PGM','ENO','PYK','PPS','PDH','G6PDH2r,PGL,GND','RPE','RPI','TKT1,TKT2','TALA','ICL,MALS','PFL','LDH_D','ALCD2x,ACALD']
    carbon_sources = ['g6p',  'xu5p_D', 'r5p',
                      'f6p',  'dhap',   '2pg',
                      '6pgc', 'pyr',    'ac',
                      'electrons']
    
    core_model = init_wt_model('core', {}, BM_lower_bound=0.1)
    #knockin_reactions(core_model, 'EDD,EDA', 0, 1000)
    knockin_reactions(core_model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc,EX_pyr', 0, 0)

    wt_model = clone_model(core_model)
    if knockins is not None:
        knockin_reactions(wt_model, knockins, 0, 1000)
    
    sys.stdout.write("There are %d single knockouts\n" % len(single_ko_list))
    sys.stdout.write("There are %d carbon sources: %s\n" % (len(carbon_sources), ', '.join(carbon_sources)))
    yield_data = np.zeros((len(single_ko_list), len(single_ko_list), len(single_ko_list), len(carbon_sources)))
    slope_data = np.zeros((len(single_ko_list), len(single_ko_list), len(single_ko_list), len(carbon_sources)))
    
    for i, ko1 in enumerate(single_ko_list):
        sys.stdout.write('delta-' + ko1 + ' and : ')
        for j, ko2 in enumerate(single_ko_list):
            sys.stdout.write(ko2 + ', ')
            for l, ko3 in enumerate(single_ko_list):
                sys.stdout.write(ko3 + ', ')
                for k, carbon_source in enumerate(carbon_sources):
                    temp_model = clone_model(wt_model)

                    if carbon_source == 'electrons':
                        knockin_reactions(temp_model, 'RED', 0, carbon_uptake_rate*2)
                    else:
                        # find out how many carbon atoms are in the carbon source
                        # and normalize the uptake rate to be in units of mmol carbon-source / (gDW*h) 
                        nC = 0
                        for cs in carbon_source.split(','):
                            met = wt_model.metabolites[wt_model.metabolites.index(cs + '_c')]
                            nC += met.formula.elements['C']
                        uptake_rate = carbon_uptake_rate / float(nC)
            
                        for cs in carbon_source.split(','):
                            set_exchange_bounds(temp_model, cs, lower_bound=-uptake_rate)
                    
                    if ko1 != '':
                        knockout_reactions(temp_model, ko1)
                    if ko2 != '' and ko2 != ko1:
                        knockout_reactions(temp_model, ko2)
                    if ko3 != '' and ko3 != ko1 and ko3 != ko2:
                        knockout_reactions(temp_model, ko3)
                        
                    yield_data[i,j,l,k] = OptKnock(temp_model).solve_FBA() or np.nan
                    
                    if target_reaction is not None and np.isnan(yield_data[i,j,l,k]):
                        slope_data[i,j,l,k] = np.nan
                    else:
                        slope_data[i,j,l,k] = OptKnock(temp_model).get_slope(target_reaction)
        sys.stdout.write('\n')

    savemat('res/data3D_%s.mat' % title,
            {'yields': yield_data, 'slopes':slope_data,
            'cs': carbon_sources, 'ko': single_ko_list}, oned_as='row')

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    main()
