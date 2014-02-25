#!/usr/bin/python
import logging
import numpy as np
import matplotlib.pyplot as plt

from src import models
from src import analysis_toolbox

def main():
    #draw_PPP('rubisco_single',
    #         target_reaction='RBC',
    #         knockins=['RBC,PRK'],
    #         carbon_sources=['g6p', '6pgc', 'xu5p_D', 'dhap'],
    #         knockouts=['RPI', 'GAPD'])

    #draw_PPP('MOG',
    #         target_reaction='MCS',
    #         knockins=['MCS,MCL'],
    #         carbon_sources=['g6p', 'ac', 'succ'],
    #         knockouts=['CS', 'PYK'])

    draw_PPP('sedoheptulose_bypass',
             target_reaction='RBC',
             knockins=['RBC,PRK', 'RBC,PRK,SBP,SBA'],
             carbon_sources=['xu5p_D', 'dhap', 'electrons'],
             knockouts=['G6PDH2r,PFK', 'G6PDH2r,PFK,TALA', 'G6PDH2r,TALA'],
             rbc_upper_bound=100)
    
    draw_PPP('rubisco_double',
             target_reaction='RBC',
             knockins=['RBC,PRK'],
             carbon_sources=[ 'xu5p_D', 'dhap', 'pyr'],
             knockouts=['RPI,G6PDH2r', 'PFK,G6PDH2r', 'PGM,G6PDH2r'],
             rbc_upper_bound=100)

    draw_PPP('rubisco_triple',
             target_reaction='RBC',
             knockins=['RBC,PRK'],
             carbon_sources=['xu5p_D,pyr', 'dhap,pyr', 'pyr'],
             knockouts=['PGM,PFK,G6PDH2r', 'PGM,RPI,G6PDH2r'],
             rbc_upper_bound=15)

    draw_PPP('rubisco_triple_SI',
             target_reaction='RBC',
             knockins=['RBC,PRK'],
             carbon_sources=['xu5p_D,ac', 'dhap,ac', 'ac'],
             knockouts=['ICL,MALS,PFK,G6PDH2r', 'ICL,MALS,RPI,G6PDH2r'],
             rbc_upper_bound=100)
             
def draw_PPP(title, target_reaction, knockins, carbon_sources, knockouts,
             rbc_upper_bound=100):
    for ki in knockins:
        wt_model = models.init_wt_model('ecoli_core', {}, BM_lower_bound=0.1)
        models.knockin_reactions(wt_model, ki)
        models.knockin_reactions(wt_model, 'EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc', 0, 0)

        carbon_uptake_rate = 50
        
        Ny = len(carbon_sources)
        Nx = len(knockouts)
        
        fig, axarr = plt.subplots(Ny, Nx, figsize=(Nx * 5, Ny * 5))
        
        for y, carbon_source in enumerate(carbon_sources):
            for x, knockout in enumerate(knockouts):
                logging.info('Calculating PPP for delta-%s with %s as a carbon source'
                             % (knockout, carbon_source))

                temp_model = models.clone_model(wt_model)
                if carbon_source == 'electrons':
                    models.knockin_reactions(temp_model, 'RED', 0, carbon_uptake_rate*2)
                else:
                    # find out how many carbon atoms are in the carbon source
                    # and normalize the uptake rate to be in units of mmol carbon-source / (gDW*h) 
                    nC = 0
                    for cs in carbon_source.split(','):
                        met = wt_model.metabolites[wt_model.metabolites.index(cs + '_c')]
                        nC += met.formula.elements['C']
                    uptake_rate = carbon_uptake_rate / float(nC)
        
                    for cs in carbon_source.split(','):
                        models.set_exchange_bounds(temp_model, cs, lower_bound=-uptake_rate)

                #if knockout == 'TKT1,TKT2':
                #    set_exchange_bounds(temp_model, 'e4p', lower_bound=-1)

                ax = axarr[y, x]
                ax.set_title(r'$v$(%s) = %.1f [mmol g(DW)$^{-1}$ h$^{-1}$]' %
                             (carbon_source, uptake_rate), fontsize=12)
                ax.set_xlabel(r'Biomass production [h$^{-1}$]', fontsize=12)
                ax.set_ylabel(r'%s flux [mmol g(DW)$^{-1}$ h$^{-1}$]' % target_reaction, fontsize=12)

                wt_PPP, wt_slope = analysis_toolbox.get_PPP(temp_model, target_reaction)
                ax.fill_between(wt_PPP[:,0].flat,
                                [max(wt_PPP[i,1], 0) for i in xrange(wt_PPP.shape[0])],
                                [min(wt_PPP[i,2], rbc_upper_bound) for i in xrange(wt_PPP.shape[0])],
                                facecolor='#E0E0E0', linewidth=0)

                models.knockout_reactions(temp_model, knockout)
                ko_PPP, slope = analysis_toolbox.get_PPP(temp_model, target_reaction)
                if slope is None:
                    ko_text = 'Not feasible at any %s flux' % target_reaction
                    ko_color = '#FF9073'
                else:
                    slope = np.round(slope, 1)
                    ko_text = ('Slope = %g' % slope)
                    if slope > 0:
                        ko_color = '#00B64F'
                    else:
                        ko_color = '#FF7060'

                ax.fill_between(ko_PPP[:,0].flat, ko_PPP[:,1].flat, ko_PPP[:,2].flat,
                                facecolor=ko_color, linewidth=1)
                
                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0, rbc_upper_bound)
                ax.text(0.5, 0.7*rbc_upper_bound, carbon_source, color='black', fontsize=15, ha='center')
                ax.text(0.5, 0.6*rbc_upper_bound, r'$\Delta$%s' % knockout, color='black', fontsize=15, ha='center')
                ax.text(0.5, 0.5*rbc_upper_bound, ko_text, color=ko_color, fontsize=15, ha='center')
        
        fig.tight_layout()
        fig.savefig('res/PPP_%s_%s.svg' % (ki, title))

if __name__ == "__main__":
    main()
