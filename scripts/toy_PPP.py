# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 16:44:51 2014

@author: eladn
"""

#!/usr/bin/python
import matplotlib.pyplot as plt

from src import models
from src import analysis_toolbox

def draw_PPP():
    wt_model = models.init_wt_model('toy2', {}, BM_lower_bound=0)

    target_reaction = 'r6'
    target_upper_bound = 2
    uptake_reactios = ['r1', 'r2']
    knockouts = 'r3,r5'
    
    Nrows = len(uptake_reactios)

    fig, ax = plt.subplots(Nrows, 2, figsize=(8,Nrows*4), sharex=True)

    for i, uptake_reaction in enumerate(uptake_reactios):
        print 'Uptake reaction: %s' % uptake_reaction
        
        temp_model = models.clone_model(wt_model)
        other_uptake_reactions = set(uptake_reactios).difference([uptake_reaction])
        models.knockout_reactions(temp_model, ','.join(other_uptake_reactions))

        wt_PPP, wt_slope = analysis_toolbox.get_PPP(temp_model, target_reaction)
        ax[i,0].fill_between(wt_PPP[:,0].flat, wt_PPP[:,1].flat, wt_PPP[:,2].flat,
                        facecolor='#BBBBBB', linewidth=1)
        ax[i,0].set_xlim(0.0, 1.0)
        ax[i,0].set_ylim(0, target_upper_bound)
        ax[i,0].set_title('wild-type')
        ax[i,0].set_ylabel('Target reaction flux')
        ax[i,0].set_xlabel('Biomass flux')

        models.knockout_reactions(temp_model, knockouts)
        ko_PPP, slope = analysis_toolbox.get_PPP(temp_model, target_reaction)
        ax[i,1].fill_between(ko_PPP[:,0].flat, ko_PPP[:,1].flat, ko_PPP[:,2].flat,
                        facecolor='#BBBBBB', linewidth=1)
        
        ax[i,1].set_xlim(0.0, 1.0)
        ax[i,1].set_ylim(0, target_upper_bound)
        ax[i,1].set_title('knockout of V1 and V3')
        ax[i,1].set_xlabel('Biomass flux')

    fig.savefig('res/toy2.svg')     
if __name__ == "__main__":
    draw_PPP()
