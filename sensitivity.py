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
from index_names import Species, Compartments, Mozzie_labels, Treatments
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
treatment_file = './stored/treatment_params.json'
scenario_list = ['BaseTreat_BaseAdherence','HighTreat_BaseAdherence']
treatment_policy_list = [Treatments.PLD, Treatments.Taf] #Run for these primary policies
adherence_sensitivity = {Treatments.PLD: np.linspace(0.5,1.0,11), Treatments.Taf: [1]}
ptreat_sensitivity = np.linspace(0.4,1.0,4)

# MtoH_sensitivity = np.linspace(1.2,1.6,9)
# treatment_policy_list = [Treatments.Baseline]
nsims = 2
offset = 3
#Added for p. vivax-only research. Changes implemented at the start (year 0), year 4, year 8.
scenario_changes_year = [0]

policies = {}
policies[Treatments.Baseline] = {"G6PD_maxes": [[0.7, 1.0],[0.3,1.0]], "treatments": [Treatments.ASMQ, Treatments.PLD], "blood_stage": Treatments.ASMQ} #Treatment regimens by G6PD test output

policies[Treatments.PLD] = {"G6PD_maxes": [[0.3, 1.0],[0.3, 1.0]], "treatments": [Treatments.PG6PD, Treatments.PLD], "blood_stage": Treatments.ASMQ}
# policies[Treatments.PHD] = {"G6PD_maxes": [0.3, 0.7, 1.0], "treatments": [Treatments.PG6PD, Treatments.PLD, Treatments.PHD]}
policies[Treatments.Taf] = {"G6PD_maxes": [[0.3, 0.7, 1.0],[0.3, 0.7, 1.0]], "treatments": [Treatments.PG6PD, Treatments.PLD, Treatments.Taf], "blood_stage": Treatments.ASMQ}

in_parallel = True

if __name__ == '__main__':

        #baseline parameters
        calibrated_params_base, ics = model_params.use_calibrated_params(prov=scenario_list[0],prov_file=prov_file,treatment_file=treatment_file)
        it_dict = model_params.initialise_dict()
        it_dict_baseline = model_params.update_dict(it_dict) #uses default first value for treatment rates / coverage
        it_dict_baseline.update(calibrated_params_base)
        params_baseline = model_params(**it_dict_baseline, **{"policy": policies[Treatments.Baseline]})
        params_baseline.number_repeats = nsims #Override model_params repeats
        params_baseline.number_offset = offset

        time_changes_year = [[year + params_baseline.burnin_years] for year in scenario_changes_year]

        #Set up multiprocessing
        # nsims = params_baseline.number_repeats
        # max_workers = int(multiprocessing.cpu_count()/nsims)-1  # Keep at least one thread free on personal computer
        # print(max_workers)
        max_workers=None

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:

            
            params_dict = {Treatments.Baseline: params_baseline}
            it_dict_dict = {Treatments.Baseline: it_dict_baseline}

            #Get parameters for each treatment option / scenario
            for primary_treatment in treatment_policy_list:
                print("primary treatment="+primary_treatment.name)

                for adh in adherence_sensitivity[primary_treatment]:
                    adh = round(adh, 2)
                    for ptreat_vivax in ptreat_sensitivity:
                        ptreat_vivax = round(ptreat_vivax, 2)
                        print("adh = "+str(adh)+", ptreat = "+str(ptreat_vivax))

                        for scenario in scenario_list[1:]:

                            calibrated_params_change, _ = model_params.use_calibrated_params(prov=scenario,prov_file=prov_file,treatment_file=treatment_file)
                            it_dict_change = model_params.update_dict(it_dict_baseline)
                            it_dict_change.update(calibrated_params_change)
                            params_treatment_change = model_params(**it_dict_change, **{"policy": policies[primary_treatment]})

                            #update params for sensitivity
                            params_treatment_change.pTreat = [params_treatment_change.pTreat[0], ptreat_vivax]
                            params_treatment_change.treatment_params[primary_treatment]["adherence"] = adh

                            #update corresponding dict
                            it_dict_change.update(pTreat = params_treatment_change.pTreat, treatment_params = params_treatment_change.treatment_params)
                            params_dict[primary_treatment] = params_treatment_change
                            it_dict_dict[primary_treatment] = it_dict_change

                            params_list = [params_dict[Treatments.Baseline],params_dict[primary_treatment]]
                            it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_dict[primary_treatment]]

                            #for file saving
                            file_id = "pTreat_"+str(ptreat_vivax)+"_adh_"+str(adh)
                            for time_change in time_changes_year:
                                time_changes_day = [int(year * params_baseline.days_in_year / params_baseline.time_day_step) for year in time_change]
                                if in_parallel:
                                    executor.submit(
                                        sim_codes.do_iterate, params_list, it_dict_list, ics, file_id, time_changes_day, in_parallel, False)
                                else:
                                    sim_codes.do_iterate(params_list, it_dict_list, ics, file_id, time_changes_day, in_parallel, False)
                                    #  sim_codes.do_iterate(params_list = [params_dict[Treatments.Baseline],params_dict[primary_treatment]], it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_dict[primary_treatment]], ics=ics, scenario=file_id, time_changes=time_changes_day, in_parallel=in_parallel, baseline_file=False)
                                # gc.collect()


# if __name__ == '__main__':

#         #baseline parameters
#         calibrated_params_base, ics = model_params.use_calibrated_params(prov=scenario_list[0],prov_file=prov_file,treatment_file=treatment_file)
#         it_dict = model_params.initialise_dict()
#         it_dict_baseline = model_params.update_dict(it_dict) #uses default first value for treatment rates / coverage
#         it_dict_baseline.update(calibrated_params_base)
#         params_baseline = model_params(**it_dict_baseline, **{"policy": policies[Treatments.Baseline]})
#         params_baseline.number_repeats = nsims #Override model_params repeats

#         time_changes_year = [[year + params_baseline.burnin_years] for year in scenario_changes_year]

#         #Set up multiprocessing
#         # nsims = params_baseline.number_repeats
#         # max_workers = int(multiprocessing.cpu_count()/nsims)-1  # Keep at least one thread free on personal computer
#         # print(max_workers)
#         max_workers=None

#         with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
#             # for adh in adherence_sensitivity:
#             #     for ptreat_vivax in ptreat_sensitivity:
#                 for MtoH in MtoH_sensitivity:
            
#                     params_dict = {Treatments.Baseline: params_baseline}
#                     it_dict_dict = {Treatments.Baseline: it_dict_baseline}

#                     #Get parameters for each treatment option / scenario
#                     for primary_treatment in treatment_policy_list:
#                         print("primary treatment="+primary_treatment.name)
#                         for scenario in scenario_list[0:1]:

#                             calibrated_params_change, _ = model_params.use_calibrated_params(prov=scenario,prov_file=prov_file,treatment_file=treatment_file)
#                             it_dict_change = model_params.update_dict(it_dict_baseline)
#                             it_dict_change.update(calibrated_params_change)
#                             params_treatment_change = model_params(**it_dict_change, **{"policy": policies[primary_treatment]})

#                             #update params for sensitivity
#                             # params_treatment_change.pTreat = [params_treatment_change.pTreat[0], ptreat_vivax]
#                             # params_treatment_change.treatment_params[primary_treatment]["adherence"] = adh
#                             #__________
#                             params_treatment_change.mozzie_human_pop_ratio = MtoH
#                             params_treatment_change.mozzie_pop = params_treatment_change.human_population * params_treatment_change.mozzie_human_pop_ratio
#                             file_id = "MtoH"+str(round(MtoH,3))
#                             #--------

#                             #update corresponding dict
#                             # it_dict_change.update(pTreat = params_treatment_change.pTreat, treatment_params = params_treatment_change.treatment_params)
#                             params_dict[primary_treatment] = params_treatment_change
#                             it_dict_dict[primary_treatment] = it_dict_change

#                             # params_list = [params_dict[Treatments.Baseline],params_dict[primary_treatment]]
#                             # it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_dict[primary_treatment]]
#                             params_list = [params_dict[Treatments.Baseline]]
#                             it_dict_list = [it_dict_dict[Treatments.Baseline]]

#                             #for file saving
#                             # file_id = "pTreat_"+str(round(ptreat_vivax,2))+"_adh_"+str(round(adh,2))
#                             # print(file_id)
#                             for time_change in time_changes_year:
#                                 time_changes_day = [int(year * params_baseline.days_in_year / params_baseline.time_day_step) for year in time_change]
#                                 #______
#                                 time_changes_day = []
#                                 #-------
#                                 if in_parallel:
#                                     executor.submit(
#                                         sim_codes.do_iterate, params_list, it_dict_list, ics, file_id, time_changes_day, in_parallel, False)
#                                 else:
#                                     sim_codes.do_iterate(params_list, it_dict_list, ics, file_id, time_changes_day, in_parallel, False)
#                                     #  sim_codes.do_iterate(params_list = [params_dict[Treatments.Baseline],params_dict[primary_treatment]], it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_dict[primary_treatment]], ics=ics, scenario=file_id, time_changes=time_changes_day, in_parallel=in_parallel, baseline_file=False)
#                                 # gc.collect()