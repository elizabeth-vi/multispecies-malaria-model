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
#from parameter_values import c_vec, pN_vec

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
#Added for p. vivax-only research. Changes implemented at the start (year 0), year 4, year 8.
#time_treatment_changes = [int(change_year * days_in_year / self.time_day_step) for change_year in [0, 2, 4]]
treatment_changes_year = [0, 0.1, 0.2]
treatment_file = './stored/treatment_params.json'

in_parallel = False

if __name__ == '__main__':

    for prov_name in prov_list_epidemics:
    
        #baseline parameters
        calibrated_params_base, ics = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file)
        it_dict = model_params.initialise_dict()
        it_dict_baseline = model_params.update_dict(it_dict) #uses default first value for treatment rates / coverage
        it_dict_baseline.update(calibrated_params_base)
        params_baseline = model_params(**it_dict_baseline)

        #baseline run
        sim_codes.do_iterate([params_baseline, params_baseline], [it_dict_baseline, it_dict_baseline], ics, prov_name, "Baseline", in_parallel)
        gc.collect()
        # duration = str(int(params_baseline.time_day_end))
        # baseline_file = "./stored/results_variables/duration_"+duration+"/baseline_timechange"+duration+"_duration"+duration+".json"

        for treatment_scenario in treatment_options_list:
        
            
            
            treatment_rates = [0]#pN_vec #assume unchanged when treatment changes - only used for size ************
            coverage_scenarios = [0]#c_vec #assume unchanged**********************************

            # for iterate1 in range(len(pN_vec)):
            # treatment scenarios
            for iterate_treat in range(len(treatment_rates)):
                #for iterate2 in range(len(c_vec)):
                for iterate_cov in range(len(coverage_scenarios)):

                    #parameters after changing radical cure treatment
                    calibrated_params_change, _ = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file,treatment=treatment_scenario,treatment_file=treatment_file)
                    it_dict_change = model_params.update_dict(it_dict, treatment_changes_year, iterate_treat, iterate_cov)
                    it_dict_change.update(calibrated_params_change)
                    params_treatment_change = model_params(**it_dict_change)


                    sim_codes.do_iterate([params_baseline, params_treatment_change], [it_dict_baseline, it_dict_change], ics, prov_name, treatment_scenario, in_parallel, baseline_file=True)
                    gc.collect()
