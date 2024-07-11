#For parameters related to code (e.g. timestep) and calculated 

# import sys
# import path as Path

# # directory reach
# directory = Path.Path(__file__).abspath()
 
# # setting path
# sys.path.append(directory.parent.parent)
 

from index_names import Species, Compartments, Treatments
import json
import numpy as np
#from parameter_values import *
import math

class model_params(object):
    def __init__(self, **kwargs):
        """
        Merge hardcoded default values (from json file) with user-specified parameter values
        :param kwargs:
        """
        # simulation parameters
        self.number_events = 12
        self.number_repeats = 10 ###### number of sims
        self.number_offset = 10 #offset for simulation seed number
        self.number_pathogens = 2
        self.number_compartments = 7

        # set time related parameters
        self.days_in_year = 365.25
        self.burnin_years = 5
        duration_years = 10 #excludes burn-in
        self.time_day_start = 0
        self.time_day_end = int((duration_years + self.burnin_years) * self.days_in_year) #Set to 1 years (previously 10 years). Add +1 or similar when used to make inclusive
        if 'time_day_step' in kwargs:
            self.time_day_step = kwargs['time_day_step']
        else:
            self.time_day_step =1.0  # looks like this needs to be about ~0.1 for vivax: notable differences between stochastic and deterministic solutions if 0.5, less but still some with 0.2
        self.time_vec = np.arange(start=self.time_day_start, stop=self.time_day_end, step=self.time_day_step)


        # #Added for p. vivax-only research. Changes implemented at time specified by time_treatment_changes, eg [0, 730, 1461]
        # if "treatment_changes_year" in kwargs:
        #     treatment_changes_year = kwargs["treatment_changes_year"]
        # else:
        #     treatment_changes_year = [duration_years]
        # # self.time_day_treatment_changes = [int(year * days_in_year) for year in treatment_changes_year]
        # self.time_treatment_changes = [int(year * days_in_year / self.time_day_step) for year in treatment_changes_year]

        # convert time to units of time_day_step
        self.time_start = int(self.time_day_start / self.time_day_step) #REMOVED ROUNDING, just floor it through int()
        self.time_end = int(self.time_day_end / self.time_day_step)

        # calculate time vectors needed for deterministic model
        self.t_det = np.arange(start=self.time_start, stop=self.time_end, step=0.5 / self.time_day_step)  # don't need to solve the odes with such a fine timestep
        self.time_vec_det = self.t_det * self.time_day_step

        self.human_population = 1000
        self.mozzie_human_pop_ratio = 1.35  # i.e. number of mosquitoes for each human
        self.period = self.days_in_year
        self.G6PD_band_ends = [0.3, 0.7, 1.0]
        self.hyp_snapshot_years = range(int(self.burnin_years),int(duration_years+self.burnin_years)+1)
        self.hyp_snapshot_days = [int(year*self.days_in_year) for year in self.hyp_snapshot_years] 
        self.hyp_snapshot_days[-1] += -1

        self.policy = {}
        self.policy["G6PD_maxes"] = [1.0]
        self.policy["treatments"] = [Treatments.Baseline]
        self.treatment_params = {}

        # transmission model parameters

        # update default values with preference to user-specified ones
        args, self.param_mins, self.param_maxes = self.defaults()
        args.update(kwargs)

        for key, value in args.items():
            self.__setattr__(key, value)
        self.param_mins['period'] = self.period - 0.01 * self.period
        self.param_maxes['period'] = self.period + 0.01 * self.period

        self.mozzie_pop = self.human_population * self.mozzie_human_pop_ratio  # not sure if this belongs here or elsewhere

        # setting the falciparum values to zero to prevent the `L` compartment being used
        not_falciparum_params = ['kappa', 'nu', 'pA', 'ph', 'pL', 'pP', 'alpha_hyp', 'mu_hyp', 'nu_hyp']
        for name in not_falciparum_params:
            self.__setattr__(name, [0.0, self.__getattribute__(name)])
            self.param_mins[name] = [0.0, self.param_mins[name]]
            self.param_maxes[name] = [0.0, self.param_maxes[name]]

        self.calculate_interactions()  # make sure default interactions are calculated

        # calculate new transmission probabilities
        self.epsilon_x = [None] * self.number_pathogens
        self.hat_epsilon_x = [None] * self.number_pathogens
        for sp in [Species.falciparum, Species.vivax]:  # or range(number_pathogens) if this doesn't work
            self.epsilon_x[sp] = [self.epsilonH[sp] * self.zetaI[sp], self.epsilonH[sp] * self.zetaA[sp], self.epsilonH[sp] * self.zetaT[sp], self.epsilonH[sp] * self.zetaG[sp]]
            self.hat_epsilon_x[sp] = [self.hatepsilonH[sp] * self.zetaI[sp], self.hatepsilonH[sp] * self.zetaA[sp], self.hatepsilonH[sp] * self.zetaT[sp], self.hatepsilonH[sp] * self.zetaG[sp]]

        if type(self.delta0) == list:
            self.delta0 = sum(self.delta0) / len(self.delta0)  # use arithematic average

    def defaults(self, filename=r"parameter_ranges.json"):
        """
        Pull out default values from a file in json format.
        :param filename: json file containing default parameter values, which can be overridden by user specified values
        :return:
        """
        param_min = dict()
        param_max = dict()
        with open(filename) as json_file:
            json_data = json.load(json_file)
        for key, value in json_data.items():
            # if this is a rate, multiply by `time_day_step`
            _tmp = value["exp"]
            _tmp_min = value["min"]
            _tmp_max = value["max"]

            if value["units"] == "per day":
                if type(_tmp) == list:
                    _tmp = [_tmp[idx] * self.time_day_step for idx in range(len(_tmp))]
                    _tmp_min = [_tmp_min[idx] * self.time_day_step for idx in range(len(_tmp_min))]
                    _tmp_max = [_tmp_max[idx] * self.time_day_step for idx in range(len(_tmp_max))]
                else:
                    _tmp = _tmp * self.time_day_step
                    _tmp_min = _tmp_min * self.time_day_step
                    _tmp_max = _tmp_max * self.time_day_step
            elif value["units"] == "day":
                if type(_tmp) == list:
                    _tmp = [_tmp[idx] / self.time_day_step for idx in range(len(_tmp))]
                    _tmp_min = [_tmp_min[idx] / self.time_day_step for idx in range(len(_tmp_min))]
                    _tmp_max = [_tmp_max[idx] / self.time_day_step for idx in range(len(_tmp_max))]
                else:
                    _tmp = _tmp / self.time_day_step
                    _tmp_min = _tmp_min / self.time_day_step
                    _tmp_max = _tmp_max / self.time_day_step

            json_data[key] = _tmp
            param_min[key] = _tmp_min
            param_max[key] = _tmp_max
        return json_data, param_min, param_max

    def interactions_off(self):
        """
        Set for the interaction parameters between falciparum and vivax such that there is none.
        :return:
        """
        interaction_params = {"a": [1.0, 1.0], "g": 1.0, "fv": 0.0, "ell": 1.0, "pLam": [1.0, 1.0],
                              "zf": 1.0, "zv": 1.0, "pZf": 0.0, "pZv": 0.0,
                              "qT": [1.0, 1.0], "qG": [1.0, 1.0], "qrho": 1.0, "qpsi": 1.0,
                              "hv": 1.0,
                              "s": [1.0, 1.0], "j": [1.0, 1.0], "k": [1.0, 1.0], "n": [1.0, 1.0], "u": [0.0, 0.0], "vW": [1.0, 1.0], "vY": [1.0, 1.0], "flag_entangled_treatment": False, "flag_triggering_by_pf": False}

        return interaction_params

    def interactions_on(self):
        """
        Set for the interaction parameters between falciparum and vivax such that there is both entanglement and triggering.
        :return:
        """
        interaction_params = {"a": [1.0, 1.0], "g": 1.0, "fv": 0.0, "ell": 1.0, "pLam": [1.0, 1.0],
                              "zf": 3.5, "zv": 1.0, "pZf": 0.0, "pZv": 0.0,
                              "qT": [1.0, 1.0], "qG": [1.0, 1.0], "qrho": 1.0, "qpsi": 1.0,
                              "hv": 1.0,
                              "s": [1.0, 1.0], "j": [1.0, 1.0], "k": [1.0, 1.0], "n": [1.0, 1.0], "u": [0.0, 0.0], "vW": [1.0, 1.0], "vY": [1.0, 1.0], "flag_entangled_treatment": True, "flag_triggering_by_pf": True}

        return interaction_params

    def calculate_interactions(self, **interaction_params):

        interaction_args = self.interactions_on() # turn on entanglement and triggering by pf
        interaction_args.update(interaction_params)

        # calculate new args
        new_args = dict()

        # RBC competition within-humans
        new_args['hatepsilonH'] = [a * b for a, b in
                                   zip(interaction_args['pLam'], self.epsilonH)]  # affects transmission to mozzies

        # competition between species for the vector
        new_args['u'] = interaction_args['u']
        new_args['vW'] = interaction_args['vW']
        new_args['vY'] = interaction_args['vY']

        interaction_args.update(new_args)

        for key, value in interaction_args.items():
            self.__setattr__(key, value)

    def initial_conditions(self, num_cases=[10, 10]):

        # initial condition
        if len(num_cases) < 2:
            human_initial_inf = [10, 10]
            print('Warning: insufficient number of initial case counts in `model_params.initial_conditions`')
        else:
            human_initial_inf = num_cases
        human_initial_mixed_only = 0  # number mixed infections at t=0
        mozzie_initial_inf = [self.mozzie_human_pop_ratio * human_initial_inf[0], self.mozzie_human_pop_ratio * human_initial_inf[1], 0]  # [falciparum-only, vivax-only, mixed]
        assert human_initial_mixed_only <= human_initial_inf[Species.falciparum] and human_initial_mixed_only <= \
               human_initial_inf[Species.vivax]

        # use the above to determine compartment counts
        initial_mozzie = [self.mozzie_pop - sum(mozzie_initial_inf), 0.0, mozzie_initial_inf[Species.falciparum], 0.0, mozzie_initial_inf[Species.vivax], 0.0, 0.0, 0.0, mozzie_initial_inf[Species.mixed]]

        human_initial_inf_comp_x_only = [[human_initial_inf[Species.falciparum] - human_initial_mixed_only, 0, 0, 0],
                                         [human_initial_inf[Species.vivax] - human_initial_mixed_only, 0, 0, 0]]
        human_initial_mixed_all = [[human_initial_mixed_only, 0, 0, 0],
                                   [0, 0, 0, 0],
                                   [0, 0, 0, 0],
                                   [0, 0, 0, 0]]
        initial_human_population_counts = [[self.human_population - human_initial_inf[Species.falciparum],
                                            human_initial_inf[Species.falciparum], 0, 0, 0, 0, 0, 0],
                                           [self.human_population - human_initial_inf[Species.vivax],
                                            human_initial_inf[Species.vivax], 0, 0, 0, 0, 0, 0]]  # last one is for `just_died`

        # identifying population counts for all compartment combos
        initial_human_combo_counts = [[0 for i in range(self.number_compartments)] for j in range(self.number_compartments)]
        # setting up no mixed infections ICs
        initial_human_combo_counts[Compartments.S] = initial_human_population_counts[Species.vivax][:-1].copy()
        for idx in range(self.number_compartments):
            initial_human_combo_counts[idx][0] = initial_human_population_counts[Species.falciparum][idx]
        # now fix S count
        initial_human_combo_counts[0][0] = int(self.human_population - sum(human_initial_inf))

        return human_initial_inf, human_initial_mixed_only, mozzie_initial_inf, human_initial_inf_comp_x_only, human_initial_mixed_all, initial_human_population_counts, initial_mozzie, initial_human_combo_counts

    def initialise_dict():
        it_dict = {'flag_entangled_treatment': 1}
        it_dict["flag_triggering_by_pf"] = 1
        it_dict["zf"] = 3.5

        it_dict["eta"] = [0, 0]

        it_dict["etaFSAT"] = [0, 0]
        it_dict["etaMDA"] = [0, 0]

        return it_dict

    def update_dict(it_dict, i1=0, i2=0):
        #includes hard-coded parameter values: adjust

        it_dict = dict(it_dict)
        # if treatment_changes_year:
        #     it_dict["treatment_changes_year"] = treatment_changes_year


        p_mask = 0.50
        mda1 = [270.0]  # lower bound of when MDA occurs
        mda2 = [270+30] # upper bound of when MDA occurs (i.e. 30 days of MDA)
        c_vec = [[0.3, 0.3]] # coverage scenarios
        pP_vec = [0.217]
        pG_vec = [[0.0066, 0.000006]]
        pN_vec = [[1.0, 0.843]] # treatment scenarios
        mda_coverage = [0.5]
        pN_mda_vec = [[[1.0, 0.843]]]

        MDA_vec = [False]
        FSAT_vec = [False] #focused screening and treatment

        mask_prob = p_mask * pN_vec[i1][0] + (1 - p_mask) * pN_vec[i1][1]
        it_dict["pN"] = pN_vec[i1]
        it_dict["mask_prob"] = mask_prob
        it_dict["pP"] = pP_vec[i1]
        it_dict["pG"] = pG_vec[i1]

        #it_dict["etaFSAT"] = [0, 0]
        #it_dict["etaMDA"] = [0, 0]

        # coverage scenarios
        it_dict["scenario"]=i2
        it_dict["mda_t1"] = mda1[i2]
        it_dict["mda_t2"] = mda2[i2]

        if MDA_vec[i2] == True:
            it_dict["MDA"] = True
            it_dict["etaMDA"] = [-math.log(1 - mda_coverage[i2]) / (mda2[i2] - mda1[i2]), -math.log(1 - mda_coverage[i2]) / (mda2[i2] - mda1[i2])]
        else:
            it_dict["MDA"] = False
            it_dict["etaMDA"] = [0, 0]

        if FSAT_vec[i2] == True:
            it_dict["FSAT"]=True
            it_dict["FSAT_period"] = 7
            it_dict["FSAT_exp_find"] = 10
            it_dict["localisation"] = 0.4
        else:
            it_dict["FSAT"] = False
            it_dict["etaFSAT"] = [0, 0]
            it_dict["FSAT_period"] = 7
            it_dict["FSAT_exp_find"] = 0
            it_dict["localisation"] = 0.4


        it_dict["pN_mda"] = pN_mda_vec[i1][i2]
        it_dict["mask_prob_mda"] = p_mask * pN_mda_vec[i1][i2][0] + (1 - p_mask) * pN_mda_vec[i1][i2][1]
        it_dict["c"] = [c_vec[i2][0], c_vec[i2][1]]


        return it_dict

    # read parameters from calibrated values
    def use_calibrated_params(prov,prov_file,treatment_file=None):
        params = dict()
        params["treatment_params"]=dict()

        # update the default parameter values using parameter values stored in `./sorted_calibrated_params.json`, after `parameter-play.py` processes the values in `./stored/model_calibration_params.json`, which were generated from `calibrated_to_cambodia_data.py`,
        with open(prov_file) as prov_file:
            prov_data = json.load(prov_file)

        for keys in prov_data[prov]:
            params[keys] = prov_data[prov][keys]


        with open(treatment_file) as treat_file:
            treat_data = json.load(treat_file)
            
        for treatment in treat_data:
            params["treatment_params"][Treatments[treatment]]=dict()
            #Note: formatted to have value and description. Can edit treatment_params to not have "value"
            for key in treat_data[treatment]:
                params["treatment_params"][Treatments[treatment]][key] = treat_data[treatment][key]["value"] #RC treatment only applicable to p vivax
                
                #If key (e.g. adherence) also exists in scenario, override with scenario-specific value
                if key in prov_data[prov]:
                    if treatment in prov_data[prov][key]:
                        params["treatment_params"][Treatments[treatment]][key] = prov_data[prov][key][treatment]
            # print(treatment)
            # print(params["treatment_params"][Treatments[treatment]])

        ics = params['ics']
        del params['ics']

        return params, ics

