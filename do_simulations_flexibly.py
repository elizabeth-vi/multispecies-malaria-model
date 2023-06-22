import H_abm_Mcomixing_pop as sim_codes
import numpy as np
import random
import multiprocessing
import time
import concurrent.futures
import math
from scipy.integrate import solve_ivp
from enum import IntEnum
import json

# import classes I've written
from index_names import Species, Compartments, Mozzie_labels
from disease_model import Disease
from mosquito_model import Mozzies
from human_agents import Agent
from model_params import model_params
from parameter_values import c_vec, pN_vec

start = time.time()  # for timing the code
#FOR REFERENCE *****************
#fig_dir = "../../figures/code_generated/"
#outfilename = "./stored/all_things"
#tangled_filename = "./stored/entangled_treatment"

import gc

########################### for user to change #################################
prov_file= './stored/sorted_calibrated_params2.json'
prov_list_epidemics = ['Example_Province'] # names need to exist in sorted_calibrated_params2.json
treatment_options_list = ['Primaquine_Lowdose','Primaquine_Highdose','Tafenoquine'] #names need to exist in sorted_calibrated_params2.json 
#treatment_options_list = ['Tafenoquine'] #names need to exist in sorted_calibrated_params2.json 
treatment_file = './stored/treatment_params.json'

# p_mask = 0.50
# mda1 = [270.0]  # lower bound of when MDA occurs
# mda2 = [270+30] # upper bound of when MDA occurs (i.e. 30 days of MDA)
# c_vec = [[0.3, 0.3]] # coverage scenarios
# pP_vec = [0.217]
# pG_vec = [[0.0066, 0.000006]]
# pN_vec = [[1.0, 0.843]] # treatment rates
# mda_coverage = [0.5]
# pN_mda_vec = [[[1.0, 0.843]]]
# MDA_vec = [False]
# it_dict["eta"] = [0, 0]
################################################################################

# FSAT_vec = [False] #focused screening and treatment
# it_dict["etaFSAT"] = [0, 0]
# it_dict["etaMDA"] = [0, 0]

in_parallel = False

if __name__ == '__main__':

    for treatment_scenario in treatment_options_list:

        for prov_name in prov_list_epidemics:
            # calibrated_params, ics = sim_codes.use_calibrated_params(prov=prov_name,file=prov_file)
            #baseline parameters
            calibrated_params_base, ics = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file)
            #params_baseline = model_params(**calibrated_params_base)

            #parameters after changing radical cure treatment
            calibrated_params_change, _ = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file,treatment=treatment_scenario,treatment_file=treatment_file)
            #params_treatment_change = model_params(**calibrated_params_change)
            
            it_dict = model_params.initialise_dict()
            treatment_rates = [1]#pN_vec #assume unchanged when treatment changes - only used for size ************
            coverage_scenarios = [1]#c_vec #assume unchanged**********************************

            #print(params.time_treatment_changes)

            #for iterate1 in range(len(pN_vec)):
            for iterate_treat in range(len(treatment_rates)):
                # treatment scenarios

                # mask_prob = p_mask * pN_vec[iterate1][0] + (1 - p_mask) * pN_vec[iterate1][1]
                # it_dict["pN"] = pN_vec[iterate1]
                # it_dict["mask_prob"] = mask_prob
                # it_dict["pP"] = pP_vec[iterate1]
                # it_dict["pG"] = pG_vec[iterate1]

                #for iterate2 in range(len(c_vec)):
                for iterate_cov in range(len(coverage_scenarios)):
                    # it_dict["etaFSAT"] = [0, 0]
                    # it_dict["etaMDA"] = [0, 0]

                    # coverage scenarios
                    # it_dict["scenario"]=iterate2
                    # it_dict["mda_t1"] = mda1[iterate2]
                    # it_dict["mda_t2"] = mda2[iterate2]

                    # if MDA_vec[iterate2] == True:
                    #     it_dict["MDA"] = True
                    #     it_dict["etaMDA"] = [-math.log(1 - mda_coverage[iterate2]) / (mda2[iterate2] - mda1[iterate2]), -math.log(1 - mda_coverage[iterate2]) / (mda2[iterate2] - mda1[iterate2])]
                    # else:
                    #     it_dict["MDA"] = False
                    #     it_dict["etaMDA"] = [0, 0]

                    # if FSAT_vec[iterate2] == True:
                    #     it_dict["FSAT"]=True
                    #     it_dict["FSAT_period"] = 7
                    #     it_dict["FSAT_exp_find"] = 10
                    #     it_dict["localisation"] = 0.4
                    # else:
                    #     it_dict["FSAT"] = False
                    #     it_dict["etaFSAT"] = [0, 0]
                    #     it_dict["FSAT_period"] = 7
                    #     it_dict["FSAT_exp_find"] = 0
                    #     it_dict["localisation"] = 0.4

                    # it_dict["pN_mda"] = pN_mda_vec[iterate1][iterate2]
                    # it_dict["mask_prob_mda"] = p_mask * pN_mda_vec[iterate1][iterate2][0] + (1 - p_mask) * pN_mda_vec[iterate1][iterate2][1]
                    # it_dict["c"] = [c_vec[iterate2][0], c_vec[iterate2][1]]

                    print([iterate_treat,iterate_cov])


                    ####
                    it_dict_baseline = model_params.update_dict(it_dict, iterate_treat, iterate_cov)
                    it_dict_change = model_params.update_dict(it_dict, iterate_treat, iterate_cov)
                    ###

                    # print(it_dict_baseline)
                    # print("******************************************")
                    # print(it_dict_change)                    
                    
                    
                    # calibrated_params_base.update(**it_dict_baseline)
                    
                    # calibrated_params_change.update(**it_dict_change)

                    #sim_codes.do_iterate(params, it_dict, prov_name, prov_file, treatment_scenario, in_parallel)

                    #calibrated_params_base.update(it_dict_baseline)
                    it_dict_baseline.update(calibrated_params_base)
                    params_baseline = model_params(**it_dict_baseline)

                    

                    #calibrated_params_change.update(it_dict_change)
                    it_dict_change.update(calibrated_params_change)

                    params_treatment_change = model_params(**it_dict_change)


                    # print(it_dict_baseline)
                    # print("******************************************")
                    # print(it_dict_change)

                    sim_codes.do_iterate([params_baseline, params_treatment_change], [it_dict_baseline, it_dict_change], ics, prov_name, prov_file, treatment_scenario, in_parallel)
                    gc.collect()
