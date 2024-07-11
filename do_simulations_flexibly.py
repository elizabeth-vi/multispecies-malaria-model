import H_abm_Mcomixing_pop as sim_codes
import numpy as np
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
# prov_list_epidemics = ['Example_Province'] # names need to exist in sorted_calibrated_params2.json
scenario_list = ['BaseTreat_BaseAdherence','HighTreat_BaseAdherence']
# scenario_list = ['BaseTreat_BaseAdherence','HighTreat_BaseAdherence']
# treatment_options_list = [Treatments.Baseline,Treatments.PLD,Treatments.PHD,Treatments.PG6PD,Treatments.Taf] #names need to exist in sorted_calibrated_params2.json
treatment_policy_list = [Treatments.PLD, Treatments.Taf] #Run for these primary policies
#Added for p. vivax-only research. Changes implemented at the start (year 0), year 4, year 8.
scenario_changes_year = [0, 4, 8]
# scenario_changes_year = [0]
# self.time_treatment_changes = [int(year * days_in_year / self.time_day_step) for year in treatment_changes_year]
#Scenarios [before timechange, after timechange]
# ********scenarios = [[Treatments.Baseline, Treatments.PLD], [Treatments.Baseline, Treatments.PHD], [Treatments.Baseline, Treatments.Taf]]
policies = {}
policies[Treatments.Baseline] = {"G6PD_maxes": [[0.7, 1.0],[0.3,1.0]], "treatments": [Treatments.ASMQ, Treatments.PLD], "blood_stage": Treatments.ASMQ} #Treatment regimens by G6PD test output

policies[Treatments.PLD] = {"G6PD_maxes": [[0.3, 1.0],[0.3, 1.0]], "treatments": [Treatments.PG6PD, Treatments.PLD], "blood_stage": Treatments.ASMQ}
# policies[Treatments.PHD] = {"G6PD_maxes": [0.3, 0.7, 1.0], "treatments": [Treatments.PG6PD, Treatments.PLD, Treatments.PHD]}
policies[Treatments.Taf] = {"G6PD_maxes": [[0.3, 0.7, 1.0],[0.3, 0.7, 1.0]], "treatments": [Treatments.PG6PD, Treatments.PLD, Treatments.Taf], "blood_stage": Treatments.ASMQ}

in_parallel = True

if __name__ == '__main__':

    # for prov_name in prov_list_epidemics:
    
        # num_processes = multiprocessing.cpu_count()  # Number of CPU cores
        
        #Set up multiprocessing
        with concurrent.futures.ThreadPoolExecutor() as executor:

            # #baseline parameters
            calibrated_params_base, ics = model_params.use_calibrated_params(prov=scenario_list[0],prov_file=prov_file,treatment_file=treatment_file)
            it_dict = model_params.initialise_dict()
            it_dict_baseline = model_params.update_dict(it_dict) #uses default first value for treatment rates / coverage
            it_dict_baseline.update(calibrated_params_base)
            params_baseline = model_params(**it_dict_baseline, **{"policy": policies[Treatments.Baseline]})

            time_changes_year = [[year + params_baseline.burnin_years] for year in scenario_changes_year]

            # #baseline run
            if in_parallel:
                executor.submit(
                    sim_codes.do_iterate, [params_baseline], [it_dict_baseline], ics, scenario_list[0], [], in_parallel)
            else:
                sim_codes.do_iterate([params_baseline], [it_dict_baseline], ics, scenario_list[0], [], in_parallel)
            
            # gc.collect()

            treatment_rates = [0]#pN_vec #assume unchanged when treatment changes - only used for size ************
            coverage_scenarios = [0]#c_vec #assume unchanged**********************************

            for iterate_treat in range(len(treatment_rates)):
                for iterate_cov in range(len(coverage_scenarios)):
            
                    params_dict = {Treatments.Baseline: params_baseline}
                    it_dict_dict = {Treatments.Baseline: it_dict_baseline}

                    #Get parameters for each treatment option / scenario
                    for primary_treatment in treatment_policy_list:
                        print("primary treatment=")
                        print(primary_treatment)
                        for scenario in scenario_list[1:]:

                            calibrated_params_change, _ = model_params.use_calibrated_params(prov=scenario,prov_file=prov_file,treatment_file=treatment_file)
                            it_dict_change = model_params.update_dict(it_dict_baseline, iterate_treat, iterate_cov)
                            it_dict_change.update(calibrated_params_change)
                            params_treatment_change = model_params(**it_dict_change, **{"policy": policies[primary_treatment]})

                            params_dict[primary_treatment] = params_treatment_change
                            it_dict_dict[primary_treatment] = it_dict_change

                            params_list = [params_dict[Treatments.Baseline],params_dict[primary_treatment]]
                            it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_dict[primary_treatment]]

                            for time_change in time_changes_year:
                                time_changes_day = [int(year * params_baseline.days_in_year / params_baseline.time_day_step) for year in time_change]
                                
                                if in_parallel:
                                    executor.submit(
                                        sim_codes.do_iterate, params_list, it_dict_list, ics, scenario, time_changes_day, in_parallel, False)
                                else:
                                    sim_codes.do_iterate(params_list, it_dict_list, ics, scenario, time_changes_day, in_parallel, False)
                                # gc.collect()
                
            
            # Close the pool and wait for all processes to complete
            # gc.collect()
            # pool.close()
            # pool.join()

                # #Additional baseline with higher adherence
                # calibrated_params_change, _ = model_params.use_calibrated_params(prov="BaseTreat_HighAdherence",prov_file=prov_file,treatment_file=treatment_file)
                # it_dict_change = model_params.update_dict(it_dict_baseline, iterate_treat, iterate_cov)
                # it_dict_change.update(calibrated_params_change)
                # params_treatment_change = model_params(**it_dict_change, **{"policy": policies[Treatments.Baseline]})

                # # params_dict[primary_treatment] = params_treatment_change
                # # it_dict_dict[primary_treatment] = it_dict_change
                # for time_change in time_changes_year:
                #     time_changes_day = [int(year * params_baseline.days_in_year / params_baseline.time_day_step) for year in time_change]
                #     sim_codes.do_iterate(params_list = [params_dict[Treatments.Baseline],params_treatment_change], it_dict_list = [it_dict_dict[Treatments.Baseline],it_dict_change], ics=ics, scenario="BaseTreat_HighAdherence", time_changes=time_changes_day, in_parallel=in_parallel, baseline_file=False)
                #     gc.collect()

                #     treatment_change_year = scenario_changes_year[scenario_i]
                #     treatments_used = scenarios[scenario_i]
                #     calibrated_params_change, _ = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file,scenario_changes_year = treatment_change_year, scenario=treatments_used,treatment_file=treatment_file)
                #     it_dict_change = model_params.update_dict(it_dict_baseline, treatment_change_year, treatments_used, iterate_treat, iterate_cov)
                #     it_dict_change.update(calibrated_params_change)
                #     params_treatment_change = model_params(**it_dict_change)

                #     params_list.append(params_treatment_change)
                #     it_dict_list.append(it_dict_change)                    


                # for scenario_i in range(scenario_changes_year):
                #     treatment_change_year = scenario_changes_year[scenario_i]
                #     treatments_used = scenarios[scenario_i]
                #     calibrated_params_change, _ = model_params.use_calibrated_params(prov=prov_name,prov_file=prov_file,scenario_changes_year = treatment_change_year, scenario=treatments_used,treatment_file=treatment_file)
                #     it_dict_change = model_params.update_dict(it_dict_baseline, treatment_change_year, treatments_used, iterate_treat, iterate_cov)
                #     it_dict_change.update(calibrated_params_change)
                #     params_treatment_change = model_params(**it_dict_change)

                #     params_list.append(params_treatment_change)
                #     it_dict_list.append(it_dict_change)

                # # # treatment policy scenarios
                # # for treatment_policy in treatment_policy_list:

                #     # sim_codes.do_iterate(params_list, it_dict_list, ics, prov_name, treatment_policy, in_parallel, baseline_file=False)
                #     sim_codes.do_iterate(params_list, it_dict_list, ics, time_changes=scenario_changes_year[scenario_i], in_parallel=in_parallel, baseline_file=False)
                #     gc.collect()
