import numpy as np
from copy import deepcopy
from pulp import LpProblem, LpMaximize, LpMinimize, LpVariable, LpAffineExpression,\
                 solvers, LpContinuous, LpBinary, LpStatusOptimal, lpSum, LpStatus
from draw_flux import DrawFlux
from cobra.core.Solution import Solution

M = 1000

class OptKnock(object):

    def __init__(self, model, verbose=False):
        self.model = deepcopy(model)
        self.verbose = verbose

        # locate the biomass reaction
        biomass_reactions = [r for r in self.model.reactions
                             if r.objective_coefficient != 0]
        if len(biomass_reactions) != 1:
            raise Exception('There should be only one single biomass reaction')
        self.r_biomass = biomass_reactions[0]
        
        self.has_flux_as_variables = False
    
    def create_prob(self, sense=LpMaximize, use_glpk=False):
        # create the LP
        self.prob = LpProblem('OptKnock', sense=sense)
        if use_glpk:
            self.prob.solver = solvers.GLPK()
        else:
            self.prob.solver = solvers.CPLEX(msg=self.verbose)
        if not self.prob.solver.available():
            raise Exception("CPLEX not available")    

    def add_primal_variables_and_constraints(self):
        # create the continuous flux variables (can be positive or negative)
        self.var_v = {}
        for r in self.model.reactions:
            self.var_v[r] = LpVariable("v_%s" % r.id,
                                       lowBound=r.lower_bound,
                                       upBound=r.upper_bound,
                                       cat=LpContinuous)

        # this flag will be used later to know if to expect the flux
        # variables to exist
        self.has_flux_as_variables = True
        
        # add the mass-balance constraints to each of the metabolites (S*v = 0)
        for m in self.model.metabolites:
            S_times_v = LpAffineExpression([(self.var_v[r], r.get_coefficient(m))
                                            for r in m.reactions])
            self.prob.addConstraint(S_times_v == 0, 'mass_balance_%s' % m.id)
    
    def add_dual_variables_and_constraints(self):
        # create dual variables associated with stoichiometric constraints
        self.var_lambda = dict([(m, LpVariable("lambda_%s" % m.id, 
                                               lowBound=-M,
                                               upBound=M,
                                               cat=LpContinuous))
                                for m in self.model.metabolites])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_w_U = dict([(r, LpVariable("w_U_%s" % r.id, lowBound=0, upBound=M, cat=LpContinuous))
                             for r in self.model.reactions])
        self.var_w_L = dict([(r, LpVariable("w_L_%s" % r.id, lowBound=0, upBound=M, cat=LpContinuous))
                             for r in self.model.reactions])

        # add the dual constraints:
        #   S'*lambda + w_U - w_L = c_biomass
        for r in self.model.reactions:
            S_times_lambda = LpAffineExpression([(self.var_lambda[m], coeff)
                                                 for m, coeff in r._metabolites.iteritems()
                                                 if coeff != 0])
            row_sum = S_times_lambda + self.var_w_U[r] - self.var_w_L[r]
            self.prob.addConstraint(row_sum == r.objective_coefficient, 'dual_%s' % r.id)
                                   
    def prepare_FBA_primal(self, use_glpk=False):
        """
            Run standard FBA (primal)
        """
        self.create_prob(sense=LpMaximize, use_glpk=use_glpk)
        self.add_primal_variables_and_constraints()
        self.prob.setObjective(self.var_v[self.r_biomass])

    def prepare_FBA_dual(self, use_glpk=False):
        """
            Run shadow FBA (dual)
        """
        self.create_prob(sense=LpMinimize, use_glpk=use_glpk)
        self.add_dual_variables_and_constraints()
        
        w_sum = LpAffineExpression([(self.var_w_U[r], r.upper_bound)
                                    for r in self.model.reactions if r.upper_bound != 0] +
                                   [(self.var_w_L[r], -r.lower_bound)
                                    for r in self.model.reactions if r.lower_bound != 0])
        self.prob.setObjective(w_sum)
    
    def get_reaction_by_id(self, reaction_id):
        if reaction_id not in self.model.reactions:
            return None
        reaction_ind = self.model.reactions.index(reaction_id)
        reaction = self.model.reactions[reaction_ind]
        return reaction
        
    def add_optknock_variables_and_constraints(self):
        # create the binary variables indicating which reactions knocked out
        self.var_y = dict([(r, LpVariable("y_%s" % r.id, cat=LpBinary))
                           for r in self.model.reactions])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_mu = dict([(r, LpVariable("mu_%s" % r.id, cat=LpContinuous))
                             for r in self.model.reactions])

        # equate the objectives of the primal and the dual of the inner problem
        # to force its optimization:
        #   sum_j mu_j - v_biomass = 0
        constr = (lpSum(self.var_mu.values()) - self.var_v[self.r_biomass] == 0)
        self.prob.addConstraint(constr, 'daul_equals_primal')

        # add the knockout constraints (when y_j = 0, v_j has to be 0)
        for r in self.model.reactions:
            # L_jj * y_j <= v_j
            self.prob.addConstraint(r.lower_bound * self.var_y[r] <= self.var_v[r], 'v_lower_%s' % r.id)
            # v_j <= U_jj * y_j
            self.prob.addConstraint(self.var_v[r] <= r.upper_bound * self.var_y[r], 'v_upper_%s' % r.id)
            
        # set the constraints on the auxiliary variables (mu):
        #    mu_j == y_j * (U_jj * w_u_j - L_jj * w_l_j)
        for r in self.model.reactions:
            w_sum = LpAffineExpression([(self.var_w_U[r], r.upper_bound),
                                        (self.var_w_L[r], -r.lower_bound)])

            # mu_j + M*y_j >= 0
            self.prob.addConstraint(self.var_mu[r] + M*self.var_y[r] >= 0, 'aux_1_%s' % r.id)
            # -mu_j + M*y_j >= 0
            self.prob.addConstraint(-self.var_mu[r] + M*self.var_y[r] >= 0, 'aux_2_%s' % r.id)
            # mu_j - (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.addConstraint(self.var_mu[r] - w_sum + M*(1-self.var_y[r]) >= 0, 'aux_3_%s' % r.id)
            # -mu_j + (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.addConstraint(-self.var_mu[r] + w_sum + M*(1-self.var_y[r]) >= 0, 'aux_4_%s' % r.id)

    def add_knockout_bounds(self, ko_candidates=None, num_deletions=5):
        """ 
            construct the list of KO candidates and add a constraint that
            only K (num_deletians) of them can have a y_j = 0
        """
        ko_candidate_sum_y = []
        
        if ko_candidates is None:
            ko_candidates = [r for r in self.model.reactions if r != self.r_biomass]

        for r in set(self.model.reactions).difference(ko_candidates):
            # if 'r' is not a candidate constrain it to be 'active'
            # i.e.   y_j == 1
            self.prob.addConstraint(self.var_y[r] == 1, 'active_%s' % r.id)

        # set the upper bound on the number of knockouts (K)
        #   sum (1 - y_j) <= K
        ko_candidate_sum_y = [(self.var_y[r], 1) for r in ko_candidates]
        constr = (LpAffineExpression(ko_candidate_sum_y) >= len(ko_candidate_sum_y) - num_deletions)
        self.prob.addConstraint(constr, 'number_of_deletions')

    def prepare_optknock(self, target_reaction_id, ko_candidates=None, 
                         num_deletions=5, use_glpk=False):
        # find the target reaction
        self.r_target = self.get_reaction_by_id(target_reaction_id)

        self.create_prob(sense=LpMaximize, use_glpk=use_glpk)
        self.add_primal_variables_and_constraints()
        self.add_dual_variables_and_constraints()
        self.add_optknock_variables_and_constraints()

        # add the objective of maximizing the flux in the target reaction
        self.prob.setObjective(self.var_v[self.r_target])

        self.add_knockout_bounds(ko_candidates, num_deletions)

    def prepare_optslope(self, target_reaction_id, ko_candidates=None,
                         num_deletions=5, use_glpk=False):
        # add the objective of maximizing the flux in the target reaction
        self.r_target = self.get_reaction_by_id(target_reaction_id)

        # set biomass maximum to 0
        self.r_biomass.lower_bound = 0
        self.r_biomass.upper_bound = 0
        self.r_target.lower_bound = 0
        self.r_target.upper_bound = 0

        self.create_prob(sense=LpMaximize, use_glpk=use_glpk)
        self.add_primal_variables_and_constraints()
        self.add_dual_variables_and_constraints()
        self.add_optknock_variables_and_constraints()

        # set the objective as maximizing the shadow price of v_target upper bound
        self.prob.setObjective(self.var_w_U[self.r_target] - self.var_w_L[self.r_target])

        self.add_knockout_bounds(ko_candidates, num_deletions)

    def write_linear_problem(self, fname):
        self.prob.writeLP(fname)

    def solve(self):
        self.prob.solve()

        if self.prob.status != LpStatusOptimal:
            if self.verbose:
                print "LP was not solved because: ", LpStatus[self.prob.status]
            self.solution = Solution(None)
        else:        
            self.solution = Solution(self.prob.objective.value())
            if self.has_flux_as_variables:
                self.solution.x = [self.var_v[r].varValue for r in self.model.reactions]
        self.solution.status = self.prob.status
        
        return self.solution
    
    def get_objective_value(self):
        if self.solution.status != LpStatusOptimal:
            return None
        else:
            return self.prob.objective.value()

    def print_primal_results(self, short=True):
        obj = self.get_objective_value()
        if obj is None:
            return
        print "Objective : %6.3f" % obj
        if not short:
            print "List of reactions : "
            for r in self.model.reactions:
                print "%30s (%4g <= v <= %4g) : v = %6.3f" % \
                    (r.name, r.lower_bound, r.upper_bound, self.var_v[r].varValue)

    def print_dual_results(self, short=True):
        obj = self.get_objective_value()
        if obj is None:
            return
        print "Objective : %6.3f" % obj
        if not short:
            print "List of reactions : "
            for r in self.model.reactions:
                print "%30s (%4g <= v <= %4g) : w_L = %5.3f, w_U = %5.3f" % \
                    (r.id, r.lower_bound, r.upper_bound, 
                     self.var_w_L[r].varValue, self.var_w_U[r].varValue) 
            print "List of metabolites : "
            for m in self.model.metabolites:
                print "%30s : lambda = %5.3f, " % \
                    (m.id, self.var_lambda[m].varValue)
                
    def print_optknock_results(self, short=True):
        if self.solution.status != LpStatusOptimal:
            return
        print "Objective : %6.3f" % self.prob.objective.value()
        print "Biomass rate : %6.3f" % self.var_v[self.r_biomass].varValue
        print "Sum of mu : %6.3f" % np.sum([mu.varValue for mu in self.var_mu.values()])
        print "Knockouts : "
        print '   ;   '.join(['"%s" (%s)' % (r.name, r.id) for r, val in self.var_y.iteritems() if val.varValue == 0]) 
        if not short:
            print "List of reactions : "
            for r in self.model.reactions:
                print '%25s (%5s) : %4g  <=  v=%5g  <=  %4g ; y = %d ; mu = %g ; w_L = %5g ; w_U = %5g' % \
                    ('"' + r.name + '"', r.id,
                     r.lower_bound, self.var_v[r].varValue, r.upper_bound,
                     self.var_y[r].varValue, self.var_mu[r].varValue,
                     self.var_w_L[r].varValue, self.var_w_U[r].varValue) 
            print "List of metabolites : "
            for m in self.model.metabolites:
                print "%30s : lambda = %6.3f" % \
                    (m.id, self.var_lambda[m].varValue)

    def get_optknock_knockouts(self):
        return ','.join([r.id for r, val in self.var_y.iteritems() if val.varValue == 0])
    
    def get_optknock_model(self):
        if self.solution.status != LpStatusOptimal:
            raise Exception('OptKnock failed, cannot generate a KO model')
        
        optknock_model = deepcopy(self.model)
        knockout_reactions = [r for r, val in self.var_y.iteritems() if val.varValue == 0]
        for r in knockout_reactions:
            new_r = optknock_model.reactions[optknock_model.reactions.index(r.id)]
            new_r.lower_bound = 0
            new_r.upper_bound = 0
        return optknock_model

    def solve_FBA(self):
        self.create_prob(sense=LpMaximize)
        self.add_primal_variables_and_constraints()
        self.prob.setObjective(self.var_v[self.r_biomass])
        self.solve()
        max_biomass = self.get_objective_value()
        return max_biomass

    def solve_FVA(self, reaction_id):
        """
            Run Flux Variability Analysis on the provided reaction
        """
        self.create_prob(sense=LpMaximize)
        self.add_primal_variables_and_constraints()
        self.prob.setObjective(self.var_v[self.r_biomass])
        self.solve()
        max_biomass = self.get_objective_value()
        if max_biomass is None:
            raise Exception("Cannot run FVA because the model is infeasible")
        self.var_v[self.r_biomass].lowBound = max_biomass - 1e-5

        r_target = self.get_reaction_by_id(reaction_id)
        self.prob.setObjective(self.var_v[r_target])
        self.prob.sense = LpMaximize
        self.solve()
        max_v_target = self.get_objective_value()

        self.prob.sense = LpMinimize
        self.solve()
        min_v_target = self.get_objective_value()
        
        return min_v_target, max_v_target

    def get_PPP_data(self, reaction_id, bm_range=None):
        """
            Run FVA on a gradient of biomass lower bounds and generate
            the data needed for creating the Phenotype Phase Plane
        """
        self.create_prob(sense=LpMaximize)
        self.add_primal_variables_and_constraints()
        self.prob.setObjective(self.var_v[self.r_biomass])
        r_target = self.get_reaction_by_id(reaction_id)
        if r_target is None:
            return None

        self.solve()

        if bm_range is None:
            max_biomass = self.get_objective_value()
            if max_biomass is None:
                return None
            bm_range = np.linspace(1e-5, max_biomass - 1e-5, 50)

        self.prob.setObjective(self.var_v[r_target])

        data = []
        for bm_lb in bm_range:
            self.var_v[self.r_biomass].lowBound = bm_lb - 1e-3
            self.var_v[self.r_biomass].upBound = bm_lb + 1e-3

            self.prob.sense = LpMaximize
            self.solve()
            max_v_target = self.get_objective_value()

            self.prob.sense = LpMinimize
            self.solve()
            min_v_target = self.get_objective_value()
            
            data.append((bm_lb, min_v_target, max_v_target))
            
        return np.matrix(data)
        
    def get_slope(self, reaction_id, epsilon_bm=0.01):
        data = self.get_PPP_data(reaction_id, bm_range=[epsilon_bm])
        if data is None:
            return None
        else:
            return data[0, 1] / epsilon_bm
        
    def model_summary(self, html):
        import analysis_toolbox
        analysis_toolbox.model_summary(self.model, self.solution, html)
        
    def draw_svg(self, html):
        # Parse the SVG file of central metabolism
        drawer = DrawFlux('data/CentralMetabolism.svg')
        #drawer = DrawFlux('data/EcoliMetabolism.svg')
        drawer.ToSVG(self.model, self.solution, html)

