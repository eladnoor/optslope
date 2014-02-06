#!/usr/bin/python
from optknock import OptKnock
from models import *
import matplotlib.pyplot as plt
from analysis_toolbox import *
from html_writer import HtmlWriter
import logging
import numpy as np
from itertools import combinations
from scipy.io import savemat

def main():
    carbon_uptake_rate = 50 # mmol C / (gDW*h)
    max_number_of_knockouts = 2

    carbon_sources = ['electrons', 'g6p', 'f6p', '2pg', 'pyr', 'ac', 'succ', '6pgc', 'r5p', 'xu5p_D', 'dhap']
    #carbon_sources = ['g6p', '6pgc']

    single_ko_list = ['ENO','FBA','FBP','G6PDH2r','GAPD','GLCpts','GND','PDH','PFK','PFL','PGI','PGK','PGL','PGM','PPC','PPCK','PYK','RPE','RPI','TALA','TKT1,TKT2','TPI','EDD','EDA']
    #single_ko_list = ['ENO','FBA','FBP']

    core_model = init_wt_model('core', {}, BM_lower_bound=0.1)
    knockin_reactions(core_model, 'EDD,EDA', 0, 1000)
    knockin_reactions(core_model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc', 0, 0)

    target_reaction = 'RBC'
    extra_reactions = 'RBC,PRK'

    #target_reaction = 'DXS'
    #extra_reactions = 'DXS'

    for knockins in [ '', extra_reactions ]:
        wt_model = clone_model(core_model)
        if knockins != '':
            knockin_reactions(wt_model, knockins, 0, 1000)
        
        all_ko = []
        for r in xrange(0, max_number_of_knockouts+1):
            for ko_list in combinations(single_ko_list, r):
                all_ko.append(','.join(ko_list))
                
        sys.stdout.write("There are %d single knockouts\n" % len(single_ko_list))
        sys.stdout.write("There are %d knockout combinations\n" % len(all_ko))
        sys.stdout.write("There are %d carbon sources: %s\n" % (len(carbon_sources), ', '.join(carbon_sources)))
        yield_data = []    
        slope_data = []    
        for carbon_source in carbon_sources:
            temp_model = clone_model(wt_model)

            sys.stdout.write(carbon_source + ' ')

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
            
            yield_data_row = []
            slope_data_row = []
            
            for knockout in all_ko:
                sys.stdout.write('.')
                ko_model = clone_model(temp_model)
                if knockout != '':
                    knockout_reactions(ko_model, knockout)
                
                max_biomass_ko = OptKnock(ko_model).solve_FBA()
                
                if max_biomass_ko is None:
                    yield_data_row.append(np.nan)
                else:
                    yield_data_row.append(max_biomass_ko)

                if (max_biomass_ko is None or 
                    target_reaction not in ko_model.reactions):
                    slope_data_row.append(np.nan)
                else:
                    slope_data_row.append(OptKnock(ko_model).get_slope(target_reaction))
                    
            yield_data.append(yield_data_row)
            slope_data.append(slope_data_row)
            sys.stdout.write('\n')
        yield_data = np.array(yield_data)
        slope_data = np.array(slope_data)
        
        savemat('res/data_%s.mat' % knockins,
                {'yields': yield_data, 'slopes':slope_data,
                'cs': carbon_sources, 'ko': all_ko}, oned_as='row')

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    main()
