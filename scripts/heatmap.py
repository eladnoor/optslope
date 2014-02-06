#!/usr/bin/python
import matplotlib.pyplot as plt
import logging
import numpy as np
from itertools import combinations, product
from scipy.io import savemat

from src.optknock import OptKnock
from src.models import *
from src.analysis_toolbox import *
from src.html_writer import HtmlWriter

def main():
    #analyze('wild_type', None, None)
    analyze('rubisco', 'RBC', 'RBC,PRK', dimension=2)
    analyze('SBPase', 'SBP', 'SBP,SBA', dimension=2)
    analyze('SBPase_with_rbc', 'SBP', 'SBP,SBA,RBC,PRK', dimension=2)
    #analyze('deoxyribose', 'DXS', 'DXS')
    #analyze('MOG', 'MCS', 'MCS,MCL')

def analyze(title, target_reaction, knockins, dimension=3):
    carbon_uptake_rate = 50 # mmol C / (gDW*h)
    max_number_of_knockouts = 2

    # the following gene pairs are combined into one "knockout" since they have no intermediate flux options:
    # G6PDH2r,PGL,GND - PP-oxidative shunt
    # GAPD,PGK        - GAP dehydrogenase + 3PG kinase
    # ICL,MALS        - glyoxlyate shunt
    # ALCD2x,ACALD    - ethanol secretion
    single_ko_list = ['GLCpts',
                      'PGI',
                      'PFK',
                      'FBP',
                      'FBA',
                      'TPI',
                      'GAPD,PGK', # GAP dehydrogenase + 3PG kinase
                      'PGM',             
                      'ENO',
                      'PYK',
                      'PPS',    
                      'PDH',
                      'G6PDH2r,PGL,GND', # PP-oxidative shunt
                      'RPE',
                      'RPI',
                      'TKT1,TKT2',
                      'TALA',
                      'ICL,MALS', # glyoxlyate shunt
                      'PFL',
                      'LDH_D',
                      'ALCD2x,ACALD' # ethanol secretion
                      ]
 
    carbon_sources = ['g6p',  'xu5p_D', 'r5p',
                      'f6p',  'dhap',   '2pg',
                      '6pgc', 'pyr',    'ac',
                      'electrons']
    
    core_model = init_wt_model('core', {}, BM_lower_bound=0.1)
    #knockin_reactions(core_model, 'EDD,EDA', 0, 1000)
    knockin_reactions(core_model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc', 0, 0)

    wt_model = clone_model(core_model)
    if knockins is not None:
        knockin_reactions(wt_model, knockins, 0, 1000)
    
    sys.stdout.write("There are %d single knockouts\n" % len(single_ko_list))
    sys.stdout.write("There are %d carbon sources: %s\n" % (len(carbon_sources), ', '.join(carbon_sources)))

    n_ko = len(single_ko_list)
    yield_data = np.zeros((n_ko**dimension, len(carbon_sources)))
    slope_data = np.zeros((n_ko**dimension, len(carbon_sources)))
    
    indices = range(n_ko)
    for i, ko_list in enumerate(product(*([single_ko_list]*dimension))):
        sys.stdout.write('KO: ' + ', '.join(ko_list) + '; carbon source: ')
        for j, carbon_source in enumerate(carbon_sources):
            sys.stdout.write(carbon_source + ', ')
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

            for ko in set(ko_list): # need to be careful not to KO the same gene twice
                knockout_reactions(temp_model, ko)

            yield_data[i,j] = OptKnock(temp_model).solve_FBA() or np.nan
            
            if target_reaction is not None and np.isnan(yield_data[i,j]):
                slope_data[i,j] = np.nan
            else:
                slope_data[i,j] = OptKnock(temp_model).get_slope(target_reaction)

        sys.stdout.write('\n')
    savemat('res/data_%dD_%s.mat' % (dimension, title),
            {'yields': yield_data, 'slopes':slope_data,
            'cs': carbon_sources, 'ko': single_ko_list}, oned_as='row')

if __name__ == "__main__":
    main()
