#!/usr/bin/python

from copy import deepcopy
import matplotlib.pyplot as plt

from src.analysis_toolbox import plot_multi_PPP
from src.models import init_wt_model, knockout_reactions, knockin_reactions
from src.optknock import OptKnock
from src.html_writer import HtmlWriter

def main():
    main_html = HtmlWriter('res/fba.html')
    main_html.write('<h1>Flux Balance Analysis</h1>\n')
    
    ko_reactions = ''
    ki_reactions = ''
    PPP_reaction = ''
    
    # full model carbon sources: glc, fru, xyl__D, succ, ac, rib__D, pyr
    
    #########################################
    # Core model for testing glycolysis KOs #
    #########################################
    #model = init_wt_model('ecoli_core', {'ac' : -10}); ko_reactions = 'PGM'; ki_reactions = 'RED';
    #model = init_wt_model('ecoli_core', {'glc' : -10}); ko_reactions = 'PFK'; ki_reactions = 'RED';

    #######################################
    # Full model Rubisco optimization KOs #
    #######################################
    #model = init_wt_model('iJO1366', {'ac' : -10}); ko_reactions = 'PGM,TRSARr,HPYRRx,HPYRRy'; ki_reactions = 'RED';
    #model = init_wt_model('iJO1366', {'rib_D' : -10}); ko_reactions = 'TKT1'; ki_reactions = 'RBC,PRK'; PPP_reaction = 'RBC';
    model = init_wt_model('iJO1366', {'xyl_D' : -10}); ko_reactions = 'RPI'; ki_reactions = 'RBC,PRK'; PPP_reaction = 'RBC';
    #model = init_wt_model('iJO1366', {'xyl_D' : -10}); ko_reactions = 'G6PDH2r,PFK,F6PA,FRUK,PFK_3,DRPA'; ki_reactions = 'RBC,PRK';
    #model = init_wt_model('iJO1366', {'fru' : -10, 'rib_D' : -10}); ko_reactions = 'G6PDH2r,PFK,F6PA,FRUK,PFK_3,DRPA,TKT1,TKT2'; ki_reactions = 'PKT';
    
    ###############################
    # Shikimate generating strain #
    ###############################
    #model = init_wt_model('iJO1366', {'fru' : -10});
    #ko_reactions = 'G6PDH2r,TALA'; ki_reactions = 'EX_3dhsk_c'; PPP_reaction = 'EX_3dhsk_c';
    #ko_reactions = 'G6PDH2r,PFK,F6PA,FRUK,PFK_3,DRPA'; ki_reactions = 'PKT'; PPP_reaction = 'PKT';

    ##########################################################
    # Testing no growth when electrons are provided directly #
    ##########################################################
    #model = init_wt_model('iJO1366', {}); ki_reactions = 'RED';
    #ko_reactions = "POR5,MCITL2";
    #ko_reactions = "POR5,FTHFLi,GART,RPI"; 
    
    models = {'WT' : model}

    if ko_reactions:
        for k in models.keys():
            m = deepcopy(models[k])
            knockout_reactions(m, ko_reactions)
            models[k + ' -%s' % ko_reactions] = m

    if ki_reactions:
        for k in models.keys():
            m = deepcopy(models[k])
            knockin_reactions(m, ki_reactions)
            models[k + ' +%s' % ki_reactions] = m

    # Run the optimization for the objective reaction and medium composition
    # set in the file.
    main_html.write('<table border="1">\n')
    main_html.write('<tr><td><b>Model Name</b></td><td><b>Growth Yield</b></td></tr>\n')
    growths = {}
    for name, model in sorted(models.iteritems()):
        print "solving %50s model" % name,
        ok = OptKnock(model)
        ok.prepare_FBA_primal()
        ok.solve()
        growths[name] = ok.get_objective_value()
    
        if growths[name] is None:
            main_html.write('<tr><td>%s</td><td>infeasible</td></tr>\n' % name)
        else:
            print ': f = %.3g' % growths[name]
            main_html.write('<tr><td>')
            html = main_html.branch(name)
            main_html.write('</td><td>%.3g</td></tr>\n' % growths[name])
            html.write('<h1>Model name: %s</h1>\n' % name)
            html.write('<h2>Growth Yield: %.3g</h2>\n' % growths[name])
            ok.draw_svg(html)
            ok.model_summary(html)
    main_html.write('</table>\n')

    if PPP_reaction:
        print 'Calculating Phenotypic Phase Plane for phosphoketolase ...'
        fig, ax = plt.subplots(1, figsize=(6,6))
        plot_multi_PPP(models, PPP_reaction, ax)
        ax.set_title('Phenotypic Phase Plane')
        main_html.embed_matplotlib_figure(fig, width=400, height=400)

if __name__ == "__main__":
    main()
