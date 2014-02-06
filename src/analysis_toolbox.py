import numpy as np
import matplotlib.pyplot as plt
from optknock import OptKnock

################################################################################
#                          HTML output tools                                   #
################################################################################
def model_summary(model, solution, html):
    reaction2flux_dict = dict([(model.reactions[i], solution.x[i])
                               for i in xrange(len(model.reactions))])

    display_exchange_reactions(model, reaction2flux_dict, html)
    
    for m in model.metabolites:
        display_metabolite_reacitons(model, m, reaction2flux_dict, html)
    
def display_exchange_reactions(model, reaction2flux_dict, html):
    # metabolite
    html.write('<br />\n')
    html.write('<a name="EXCHANGE"></a>\n')
    html.write('Exchange reactions: <br />\n')
    html.write('<br />\n')
    
    titles = ['Sub System', 'Reaction Name', 'Reaction ID',
              'Reaction', 'LB', 'UB', 'Reaction Flux']

    # fluxes
    rowdicts = []
    for r in model.reactions:
        if r.subsystem not in ['', 'Exchange']:
            continue
        if abs(reaction2flux_dict[r]) < 1e-10:
            continue

        direction = np.sign(reaction2flux_dict[r])
        d = {'Sub System': 'Exchange', 'Reaction Name': r.name,
             'Reaction ID': r.id,
             'Reaction': display_reaction(r, None, direction),
             'LB': '%g' % r.lower_bound,
             'UB': '%g' % r.upper_bound,
             'Reaction Flux': '%.2g' % abs(reaction2flux_dict[r]),
             'sortkey': reaction2flux_dict[r]}

        rowdicts.append(d)
        
    # add a zero row (separating forward and backward) and sort the 
    # rows according to the net flux
    rowdicts.append({'sortkey': 0})
    rowdicts.sort(key=lambda x: x['sortkey'])

    # add color to the rows
    max_flux = max([abs(d['sortkey']) for d in rowdicts])
    rowcolors = [color_gradient(d['sortkey']/max_flux) for d in rowdicts]

    html.write_table(rowdicts, titles, rowcolors=rowcolors)

def display_metabolite_reacitons(model, m, reaction2flux_dict, html):
    # metabolite
    html.write('<br />\n')
    html.write('<a name="%s"></a>\n' % m.id)
    html.write('Metabolite name: ' + m.name + '<br />\n')
    html.write('Metabolite ID: ' + m.id + '<br />\n')
    html.write('Compartment: ' + m.compartment + '<br />\n')
    html.write('<br />\n')

    titles = ['Sub System', 'Reaction Name', 'Reaction ID',
              'Reaction', 'LB', 'UB', 'Reaction Flux', 'Net Flux']
    
    # fluxes
    rowdicts = []
    for r in m.get_reaction():
        if abs(reaction2flux_dict[r]) < 1e-10:
            continue
            
        direction = np.sign(reaction2flux_dict[r])
        net_flux = reaction2flux_dict[r] * r.get_coefficient(m)
        d = {'Sub System': r.subsystem, 'Reaction Name': r.name,
             'Reaction ID': r.id,
             'Reaction': display_reaction(r, m, direction),
             'LB': '%g' % r.lower_bound,
             'UB': '%g' % r.upper_bound,
             'Reaction Flux': '%.2g' % abs(reaction2flux_dict[r]),
             'Net Flux': '%.2g' % net_flux,
             'sortkey': -net_flux}

        rowdicts.append(d)
    
    if rowdicts == []:
        return

    # add a zero row (separating forward and backward) and sort the 
    # rows according to the net flux
    rowdicts.append({'sortkey': 0})
    rowdicts.sort(key=lambda x: x['sortkey'])

    # add color to the rows
    max_flux = max([abs(d['sortkey']) for d in rowdicts])
    rowcolors = [color_gradient(d['sortkey']/max_flux) for d in rowdicts]

    html.write_table(rowdicts, titles, rowcolors=rowcolors)

def display_reaction(r, m_bold=None, direction=1):
    """
        Returns a string representation of a reaction and highlights the
        metabolite 'm' using HTML tags.
    """
    s_left = []
    s_right = []
    for m in r.get_reactants() + r.get_products():
        if m == m_bold:
            s_met = "<a href='#%s'><b>%s</b></a>" % (m.id, m.id)
        else:
            s_met = "<a href='#%s'>%s</a>" % (m.id, m.id)
        
        coeff = r.get_coefficient(m)
        if abs(coeff) == 1:
            s_coeff = ""
        else:
            s_coeff = "%g " % abs(coeff)
        
        if coeff < 0:
            s_left += [s_coeff + s_met]
        else:
            s_right += [s_coeff + s_met]
    
    if direction == 1:
        return ' + '.join(s_left) + ' &#8651; ' + ' + '.join(s_right)
    else:
        return ' + '.join(s_right) + ' &#8651; ' + ' + '.join(s_left)

def color_gradient(x):
    """
        Returns a color in Hex-RGB format between white and red if x is positive
        or between white and green if x is negative
    """
    grad = 220 - abs(x)*80
    if x > 0:
        return '%.2x%.2x%.2x' % (255, grad, grad)
    elif x < 0:
        return '%.2x%.2x%.2x' % (grad, 255, grad)
    else:
        return '%.2x%.2x%.2x' % (100, 100, 100)

################################################################################
#                           plotting tools                                     #
################################################################################
def plot_phase(model, pivot_reaction_id):
    ok = OptKnock(model)
    r_pivot = ok.get_reaction_by_id(pivot_reaction_id)
    ok.prepare_FBA_primal()
    
    data = []
    for f in np.arange(0, 100, 1):
        ok.var_v[r_pivot].lowBound = f
        ok.var_v[r_pivot].upBound = f
        ok.solve()
        obj = ok.get_objective_value() or 0
        data.append((f, obj))
        print "RBC = %.3g, BM = %.3g" % (f, obj)
    
    data = np.matrix(data)
    fig = plt.figure()
    plt.plot(data[:,1], data[:,0], '-', figure=fig)
    return fig

def get_PPP(model, reaction_id):
    ok = OptKnock(model)
    PPP_data = ok.get_PPP_data(reaction_id)
    if PPP_data is None:
        print 'model is infeasible'
        PPP_data = np.matrix([[0, 0, 0.05], [0.05, 0, 0]])
        slope = None
    else:
        # calculate the "slope" - i.e. close to 0 biomass production, what is the
        # minimal flux in the target reaction that is required?
        slope = PPP_data[1, 1] / PPP_data[1, 0]

    return PPP_data, slope

def plot_multi_PPP(model_dict, reaction_id, ax):
    """
        Draw a comparative Phenotypic Phase Plane plot.
        Arguments:
            models - a list of cobrapy Model objects
            labels - a list of strings for use as the legend
    """
    colors = ['red', 'blue', 'green', 'magenta', 'orange', 'cyan']

    PPP_data_dict = {}
    for label, model in model_dict.iteritems():
        # verify that this model has the pivot reaction:
        if reaction_id not in model.reactions:
            print 'model "%s" does not contain the reaction "%s"' % (label, reaction_id)
            continue
        PPP_data = OptKnock(model).get_PPP_data(reaction_id, num_steps=20)

        # verify that this model is feasible (i.e. biomass yield is more than minimal threshold):
        if PPP_data is None:
            print 'model "%s" is infeasible' % label
            continue
        PPP_data_dict[label] = PPP_data        

    for i, (label, PPP_data) in enumerate(PPP_data_dict.iteritems()):
        if label == 'wild-type':
            ax.fill_between(PPP_data[:,0].flat, PPP_data[:,1].flat, PPP_data[:,2].flat,
                facecolor='grey', alpha=0.1, linewidth=0)
        else:
            ax.fill_between(PPP_data[:,0].flat, PPP_data[:,1].flat, PPP_data[:,2].flat,
                            facecolor=colors[i], alpha=0.1, linewidth=0)

    ax.set_xlabel(r'Biomass production [h$^{-1}$]')
    ax.set_ylabel(r'%s flux [mmol g(DW)$^{-1}$ h$^{-1}$]' % reaction_id)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.grid()

    for i, (label, PPP_data) in enumerate(PPP_data_dict.iteritems()):
        if label == 'wild-type':
            ax.plot(PPP_data[:,0].flat, PPP_data[:,1].flat, lw=0, label=label)
        else:
            ax.plot(PPP_data[:,0].flat, PPP_data[:,1].flat, color=colors[i], lw=1, label=label)

    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.7)

    for i, (label, PPP_data) in enumerate(PPP_data_dict.iteritems()):
        if label == 'wild-type':
            ax.plot(PPP_data[:,0].flat, PPP_data[:,2].flat, lw=0, label=label)
        else:
            ax.plot(PPP_data[:,0].flat, PPP_data[:,2].flat, color=colors[i], lw=1, label=label)

