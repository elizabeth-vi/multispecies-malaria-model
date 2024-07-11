import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.patches as mpatches
import json
import numpy as np
from model_params import model_params
import itertools
from index_names import Compartments, Mozzie_labels, Species, Treatments
import os
from collections import Counter
import time
import math

def get_filepath(folder, treatment, scenario, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    filepath = folder + treatment + "_" + scenario + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filepath

def get_filename(treatment, pTreat, adh, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    filename = treatment + "_pTreat_" + str(pTreat) + "_adh_" + str(adh) + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filename

def load_data(filepath, category):
    #Load and return data
    try:
        data_all = json.load(open(filepath, 'r'))
        data_cat = data_all.get(category)
        return data_cat
    except:
        print("Filepath "+filepath+" does not exist.")
        return
    
def get_plot_data(data, t_start, t_end, t_step, burnin_days):

    start_index = int(t_start/t_step)
    end_index = int(t_end/t_step)

    plot_times = np.arange(start=start_index*t_step-burnin_days, stop=end_index*t_step-burnin_days, step=t_step)
    plot_data = data[start_index:end_index]

    return plot_times, plot_data

def get_baseline_clinical(folders_list, max_files, plot_start, plot_end, t_step, burnin_days, treatment, scenario, timechange, duration):
    clinical_infect_i = 1
    count = 0
    sum_infect = 0
    for folder in folders_list:
        filepath = get_filepath(folder, treatment, scenario, timechange, duration)
        data_pv_outcomes = load_data(filepath, "pv_outcomes")
        for sim in range(len(data_pv_outcomes)):
            if count+1 > max_files:
                break
            count += 1
            data_infect_dict = Counter(data_pv_outcomes[sim][clinical_infect_i])
            data_infect = [data_infect_dict[t] for t in range(duration)]
            plot_times, plot_infect = get_plot_data(data_infect, plot_start, plot_end, t_step, burnin_days = burnin_days)
            sum_infect += sum(plot_infect)
            print(sum(plot_infect))
    avg_infect = round(sum_infect/count)
    # median_infect ??
    return avg_infect
            



params = model_params()
days_in_year = 365.25
treatments = [Treatments.Taf, Treatments.PLD] #Must match file names. Options: ["PLD","PHD","Taf"]]
# scenarios = ["HighTreat_BaseAdherence"]
# adherence_sensitivity = np.linspace(0.4,1.0,4)
ptreat_sensitivity = np.linspace(0.4,1.0,4)
adherence_sensitivity_treats = {Treatments.PLD: np.linspace(0.5,1.0,11), Treatments.Taf: [1]}
timechanges_year = [0]
timechanges = [[int((year + params.burnin_years)*days_in_year)] for year in timechanges_year]

duration = int(params.time_day_end) #Simulation duration
baseline = ["Baseline", "BaseTreat_BaseAdherence", [], duration]

burnin_years = params.burnin_years #years
burnin_days = burnin_years * days_in_year

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start+burnin_days, duration, params.time_day_step] #Time of start/end you want to plot, and timestep used
plt.rcParams["figure.figsize"] = (11,8.5)
plt.rcParams.update({'font.size': 16})
title_sensitivity = {Treatments.PLD: "Low-dose Primaquine Based Policy", Treatments.Taf: "Tafenoquine Based Policy"}

clinical_infect_i = 1
# folders_list = ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration5478_sims1_137579")]#"duration"+str(duration)+"_sims"+str(nsims)+"_")]
folders_list = ["./stored/results_variables/sensitivity_results/" + folder + "/" for folder in os.listdir('./stored/results_variables/sensitivity_results') if folder.startswith("duration")]
folders_list_baseline = ["./stored/results_variables/treatment_results/" + folder + "/" for folder in os.listdir('./stored/results_variables/treatment_results') if folder.startswith("duration")]

plot_ptreat_adh = True
if plot_ptreat_adh:
    # adherence_sensitivity = np.flip(adherence_sensitivity) #reverse for plotting; start with 1 and move down
    max_files = 10
    max_val = 22000
    timechange = timechanges[0]
    fig,ax = plt.subplots(len(treatments), sharex=True, gridspec_kw={"height_ratios": [len(adherence_sensitivity_treats[treatment]) for treatment in treatments]})
    fignum = 0

    # avg_baseline_infections = get_baseline_clinical(folders_list_baseline, max_files, plot_start, plot_end, t_step, burnin_days, *baseline)
    
    for treatment in treatments:
        adherence_sensitivity = np.round(np.flip(adherence_sensitivity_treats[treatment]), 2)
        
        data_sensitivity_arr = np.zeros(shape=[len(adherence_sensitivity),len(ptreat_sensitivity)])

        #iterate sensitivities
        for col in range(len(ptreat_sensitivity)):
            pTreat = ptreat_sensitivity[col]
            for row in range(len(adherence_sensitivity)):
                adh = round(adherence_sensitivity[row],2)

                sum_infect = 0

                #get data
                count=0
                for folder in folders_list:
                    filepath = folder + get_filename(treatment.name, pTreat, adh, timechange, duration)
                    data_pv_outcomes = load_data(filepath, "pv_outcomes")
                    for sim in range(len(data_pv_outcomes)):
                        if count+1 > max_files:
                            break
                        count += 1
                        data_infect_dict = Counter(data_pv_outcomes[sim][clinical_infect_i])
                        data_infect = [data_infect_dict[t] for t in range(duration)]
                        plot_times, plot_infect = get_plot_data(data_infect, timechange[0], plot_end, t_step, burnin_days = burnin_days)
                        # print(sum(plot_infect))
                        sum_infect += sum(plot_infect)

                #collate data, add text label
                avg_infect = round(sum_infect/count)
                data_sensitivity_arr[row][col] = avg_infect #avg_baseline_infections - avg_infect
                ax[fignum].text(col, row, round(data_sensitivity_arr[row][col]), ha="center", va="center", color="w", path_effects=[pe.withStroke(linewidth=4, foreground="k")])
        # max_val = max(max_val, np.max(data_sensitivity_arr)) #edit if needed
        im = ax[fignum].imshow(data_sensitivity_arr, aspect="auto", cmap = "RdBu_r", vmax=max_val)
        
        # Show all ticks and label them with the respective list entries
        ax[fignum].set_xticks(np.arange(len(ptreat_sensitivity)), labels=ptreat_sensitivity)
        plt.xlabel("Probability a symptomatic case receives treatment $(p_{treat})$")
        ax[fignum].set_yticks(np.arange(len(adherence_sensitivity)), labels=adherence_sensitivity)
        # ax[fignum].set_ylabel("Adherence")
        fig.text(.05, .5, 'Adherence', ha='center', va='center', rotation='vertical')

        ax[fignum].set_title(title_sensitivity[treatment])#treatment.name+", cumul clinical infections over "+str(round(duration/365.25 - burnin_years))+" years")
        fignum += 1
    # fig.savefig("./stored/figs/sensitivity_ptreat_adh_"+treatment.name+".png",bbox_inches='tight')
    # fig.subplots_adjust(right=1)
    fig.colorbar(im, ax=ax)
    # fig.tight_layout()
    fig.savefig("./stored/figs/sensitivity_ptreat_adh.png",bbox_inches='tight')
    

MtoH_sensitivity = np.linspace(1.2,1.6,9)


plot_MtoH = False
if plot_MtoH:
    # adherence_sensitivity = np.flip(adherence_sensitivity) #reverse for plotting; start with 1 and move down
    max_files = 1
    timechange = []
    clinical_per_year = 0
    
    for treatment in ["Baseline"]:
        fig,ax = plt.subplots()
        baseline_data_num = 100
        data_sensitivity_arr = np.zeros(shape=[len(MtoH_sensitivity)])
        # print(data_sensitivity_arr)
        print("\n"+treatment+" after "+str(burnin_years)+" years burn-in")
        print("Mosquitoes per human     |      Infections over "+str(round((duration-burnin_days)/365))+" years")

        #iterate sensitivities
        for col in range(len(MtoH_sensitivity)):
            MtoH = round(MtoH_sensitivity[col],3)
            print("\n"+str(MtoH))
            for row in range(1):
                adh = 1

                sum_infect = 0
                #get data
                count=0
                folders_list = ["./stored/results_variables/duration5478_sims1_223956/"]
                for folder in folders_list:
                    # filepath = folder + get_filename(treatment, pTreat, adh, timechange, duration)
                    filepath = folder + "Baseline_MtoH"+str(MtoH)+"_timechange[]_duration5478.json"
                    data_pv_outcomes = load_data(filepath, "pv_outcomes")
                    for sim in range(len(data_pv_outcomes)):
                        if count+1 > max_files:
                            break
                        count += 1
                        data_infect_dict = Counter(data_pv_outcomes[sim][clinical_infect_i])
                        data_infect = [data_infect_dict[t] for t in range(duration)]
                        plot_times, plot_infect = get_plot_data(data_infect, plot_start, plot_end, t_step, burnin_days = burnin_days)
                        sum_infect += sum(plot_infect)

                        #Check yearly infections for stability
                        for year in range(round((duration-burnin_days)/365)):
                            day_start = round(year*365.25 + burnin_days)
                            day_end = day_start+365
                            _, yearly_infect = get_plot_data(data_infect, day_start, day_end, t_step, burnin_days = burnin_days)
                            # print(str(year)+"____"+str(sum(yearly_infect)))


                #collate data, add text label
                avg_infect = sum_infect/count
                data_sensitivity_arr[col] = avg_infect
                ax.text(col, row, avg_infect, ha="center", va="center", color="w")
            print('{:.2f}'.format(MtoH)+"                             "+'{:.0f}'.format(avg_infect))
        # print(data_sensitivity_arr)
        
        data_sensitivity_arr = np.expand_dims(data_sensitivity_arr,axis=0)
        im = ax.imshow(data_sensitivity_arr)

        # Show all ticks and label them with the respective list entries
        # ax.set_xticks(np.arange(len(ptreat_sensitivity)), labels=ptreat_sensitivity)
        # ax.set_xlabel("pTreat")
        # ax.set_yticks(np.arange(len(adherence_sensitivity)), labels=adherence_sensitivity)
        # ax.set_ylabel("Adherence")
        ax.set_title(treatment+", cumul clinical infections over "+str(round(duration/365.25 - burnin_years))+" years")

        
        fig.tight_layout()


print("Count = "+str(count))
plt.show(block=False)
print("Plots finished. Press ENTER to continue.")
input()
plt.close()