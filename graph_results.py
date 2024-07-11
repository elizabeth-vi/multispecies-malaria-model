import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json
import numpy as np
from model_params import model_params
import itertools
from index_names import Compartments, Mozzie_labels
import os
from collections import Counter
import time
import math



def get_filepath(folder, treatment, scenario, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    filepath = folder + treatment + "_" + scenario + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filepath

def get_filename(treatment, scenario, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    filename = treatment + "_" + scenario + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filename

def load_pv_compartment_data(filepath):

    #Load and return data
    try:
        data = json.load(open(filepath, 'r'))
        vals_pv = data.get("human_pop_pv_history")[0]
        return vals_pv
    except:
        print("Filepath "+filepath+" does not exist.")
        return
    
def load_relapses(filepath):
    #Load and return data
    try:
        data = json.load(open(filepath, 'r'))
        vals_pv = data.get("pv_relapses")[0]
        return vals_pv
    except:
        print("Filepath "+filepath+" does not exist.")
        return

def load_data(filepath, result_type):
    #Load and return data
    try:
        data_all = json.load(open(filepath, 'r'))
        if any(cat in result_type for cat in ["daily","cumulative"]):
            category = "pv_outcomes"
        elif "compartment" in result_type:
            category = "human_pop_pv_history"
        elif "hypnozoite" in result_type:
            category = "hypnozoite_dist"
        else:
            print("Undefined load results type")
        data_cat = data_all.get(category)
        return data_cat
    except:
        print("Filepath "+filepath+" does not exist.")
        return

def process_data(data, results_type, duration):
    if "clinical" in results_type:
        data_infect_dict = Counter(data[1])
        data_infect = [data_infect_dict[t] for t in range(duration)] 
        if "cumulative" in results_type:
            data_infect = np.cumsum(data_infect)
    elif "relapses" in results_type:
        data_infect_dict = Counter(data[2])
        data_infect = [data_infect_dict[t] for t in range(duration)]
        if "cumulative" in results_type:
            data_infect = np.cumsum(data_infect)
    elif "compartment" in results_type:
        if "I" in results_type:
            for t in range(duration):
                data_infect = [data[t][1] for t in range(duration)]
    else:
        print("Undefined process results type")
    return data_infect

def get_plot_data(data, t_start, t_end, t_step, burnin_days, results_type):
    start_index = int(t_start/t_step)
    end_index = int(t_end/t_step)

    plot_times = np.arange(start=start_index*t_step-burnin_days, stop=end_index*t_step-burnin_days, step=t_step)
    plot_data = data[start_index:end_index]

    if "cumulative" in results_type:
        plot_data = plot_data - data[int(burnin_days/t_step)]

    return plot_times, plot_data

def moving_average(arr, n=1):
  if n == 1:
    return arr
  ret = np.cumsum(arr, dtype=float)
  ret[n:] = ret[n:] - ret[:-n]
  divisor = np.ones((len(arr)))
  divisor[:] = n
  divisor [:n] = np.arange(1, n+1)
  return ret / divisor

def getlims(plot_type):
    xlims = None
    ylims = None
    if "cumulative" in plot_type:
        if "clinical" in plot_type:
            ylims = [0, 4000]
        elif "relapses" in plot_type:
            ylims = [0, 8000]
    else:
        ylims = [0, None]

    return xlims, ylims

def print_cumul_values(plot_type, baseline_data, treatment_data, treatment, timechange, CIs):
    baseline_data = [data[-1] for data in baseline_data]
    treatment_data = [data[-1] for data in treatment_data]
    cases_averted = [baseline_data[sim] - treatment_data[sim] for sim in range(len(baseline_data))]
    percentage_cases_averted = [cases_averted[sim]/baseline_data[sim] for sim in range(len(baseline_data))]
    CI = 0.5-0.5*CIs[0], 0.5+0.5*CIs[0]

    # print("\n"+treatment+" at t="+str(round(timechange[0]/days_in_year-burnin_years))+" years")
    year = round(timechange[0]/days_in_year-burnin_years)
    if year == 0:
        if treatment == "PLD":
            print("\multirow{3}{*}{\makecell{Low-dose \\\\primaquine}}")
        elif treatment == "Taf":
            print("\multirow{3}{*}{Tafenoquine}")
    print("& "+str(year))
    result_list = [treatment_data, cases_averted, [100*pc for pc in percentage_cases_averted]]
    cat = "symptomatic" if "clinical" in plot_type else "relapses" if "relapse" in plot_type else plot_type
    result_name = ["_cases_","_averted_","_pcaverted_"]
    for result_i in range(len(result_list)):
        result = result_list[result_i]
        name = cat + result_name[result_i]
        precision = 1 if "pc" in name else None
        median = round(np.median(result), precision)
        lower_CI = round(np.quantile(result, CI[0]),precision)
        upper_CI = round(np.quantile(result, CI[1]),precision)
        # print(str(median)+" ("+str(lower_CI)+"-"+str(upper_CI)+")")
        print("& "+str(median).replace('.','·')+"\\numlabel{"+str(median).replace('.','·')+"}{num:"+name+treatment+str(year)+"} ("+str(lower_CI).replace('.','·')+"-"+str(upper_CI).replace('.','·')+")")
    if year == 8:
        print("\\\\\\hline")
    else:
        print("\\\\\\cline{2-5}")
    # print("Numbers: "+str(round(np.median(treatment_data)))+" ("+str(round(np.quantile(treatment_data, 0.025)))+"-"+str(round(np.quantile(treatment_data, 0.975)))+")")
    # print("Reduction: "+str(round(np.median(cases_averted)))+" ("+str(round(np.quantile(cases_averted, 0.025)))+"-"+str(round(np.quantile(cases_averted, 0.975)))+")")
    # print("Percentage reduction: "+str(round(100*np.median(percentage_cases_averted),1))+" ("+str(round(100*np.quantile(percentage_cases_averted, 0.025),1))+"-"+str(round(100*np.quantile(percentage_cases_averted, 0.975),1))+")")


#*********************************SCRIPT START*********************************


#Filepath name info
#[filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]


params = model_params()
days_in_year = 365.25
# moving_average_num = 7
treatments = ["PLD","Taf"] #Must match file names. Options: ["PLD","PHD","Taf"]]
scenarios = ["HighTreat_BaseAdherence"]
# scenarios = ["HighTreat_BaseAdherence"]
# timechanges = [[1826],[3287],[4748]] #[[0], [1461], [2922]] #which time changes you want to plot, in days. e.g.: [0, 730, 1461]
# timechanges = [[365]] #[[0], [1461], [2922]] #which time changes you want to plot, in days. e.g.: [0, 730, 1461]
timechanges_year = [0, 4, 8]
timechanges = [[int((year + params.burnin_years)*days_in_year)] for year in timechanges_year]
print(timechanges)
treatments_titles = ["Change to Greater Rates of Tafenoquine Treatment"] #["Change to Policy 1", "Change to Policy 2", "Change to Policy 3"]
treatments_legend = ["Primaquine","Tafenoquine"] #["Change to Policy 1", "Change to Policy 2", "Change to Policy 3"]
y_labels = {"cumulative_clinical": "Cumulative Symptomatic Cases", "cumulative_relapses": "Cumulative Relapses", "hypnozoites_mean_all": "Mean number of hypnozoites"}
# colours = ['purple','blue','red', 'k'] #for plotting
colours_list_timechange = ['purple','blue','red', 'k'] #for plotting
# colours_map_timechange = {timechanges[0]:'purple',timechanges[1]:'blue',timechanges[2]:'red', "Baseline": 'k'} #for plotting
colours_map_treatment = {"Baseline":"black", "PLD":"red", "Taf":"blue"}
colours_map_scenario = {"Baseline":"black", scenarios[0]:"blue"}#, scenarios[1]:"red", scenarios[2]:"purple"}

# timechanges_year = [[5],[9],[13]]
duration = int(params.time_day_end) #Simulation duration
baseline = ["Baseline", "BaseTreat_BaseAdherence", [], duration]

burnin_years = params.burnin_years #years
burnin_days = burnin_years * days_in_year

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start+burnin_days, duration, params.time_day_step] #Time of start/end you want to plot, and timestep used
plt.rcParams["figure.figsize"] = (11,8.5)
plt.rcParams.update({'font.size': 16})

legend_treat = ["Baseline (current practice)"]
legend_treat.extend(treatments_legend)
legend_timechange = ["Baseline (current practice)"]
legend_timechange.extend(["Change policy at t = "+str(round((x[0]-burnin_days)/days_in_year))+" years" for x in timechanges])
legend_scenario = ["Baseline (current practice)"]
legend_scenario.extend(scenarios)
folders_list = ["./stored/results_variables/treatment_results/" + folder + "/" for folder in os.listdir('./stored/results_variables/treatment_results') if folder.startswith("duration")]
print(folders_list)

# plot with median and quantiles
plot_multiple = False
if plot_multiple:
    plot_types = ["cumulative_clinical", "cumulative_relapses"]#["daily_clinical","compartment_I", "daily_relapses"]

    alpha = 0.4
    max_files = 30
    treatments_plot = treatments 
    scenarios_plot = scenarios
    CIs_plot = [0.5]
    CIs_table = [0.9]
    clinical_infect_i = 1

    for plot in plot_types:
        baseline_data_list = []
        treatment_data_list = []

        count=0
        for folder in folders_list:
            #baseline plot data
            filepath = folder+get_filename(*baseline) #Baseline
            data_pv_outcomes = load_data(filepath, plot)
            for sim in range(len(data_pv_outcomes)):
                if count+1 > max_files:
                    break
                count += 1
                data_infect = process_data(data_pv_outcomes[sim],plot, duration)
                plot_times_baseline, plot_infect = get_plot_data(data_infect, plot_start, plot_end, t_step, burnin_days = burnin_days, results_type=plot)
                baseline_data_list.append(plot_infect)
            
        baseline_data = np.array(baseline_data_list)
        median_baseline = np.quantile(baseline_data, 0.5, axis=0)
        print_cumul_values(plot, baseline_data, baseline_data, "Baseline", [0], CIs_table)

        #Plot for treatments
        for treatment in treatments_plot:
            treatment_i = treatments.index(treatment)
            for scenario in scenarios_plot:
                #Baseline plotting
                fig_infect, ax_infect = plt.subplots(1) #set up figures
                trans = ax_infect.get_xaxis_transform() #Transformation for labelling line

                #Plot median
                ax_infect.plot(plot_times_baseline/days_in_year, moving_average(median_baseline), color = colours_map_treatment["Baseline"], alpha=1, linewidth = 1, zorder = 100) 
                #Plot quantiles
                for CI in CIs_plot:
                    quants = [0.5-0.5*CI, 0.5+0.5*CI]
                    low_quant = np.quantile(baseline_data, quants[0], axis=0)
                    high_quant = np.quantile(baseline_data, quants[1], axis=0)
                    plt.fill_between(plot_times_baseline/days_in_year, moving_average(low_quant), moving_average(high_quant), color = colours_map_treatment["Baseline"], alpha=alpha, label='_nolegend_')

                #Treatment / scenario plotting for all timechanges
                for timechange_i in range(len(timechanges)):
                    timechange = timechanges[timechange_i]
                    treatment_data_list = [] #reset list
                    count = 0
                    for folder in folders_list:
                        filepath = folder + get_filename(treatment, scenario, timechange, duration)
                        data_pv_outcomes = load_data(filepath, plot)

                        for sim in range(len(data_pv_outcomes)):
                            if count+1 > max_files:
                                break
                            count += 1
                            # data_infect = [data_pv_outcomes[1].count(i) for i in range(duration)]
                            
                            data_infect = process_data(data_pv_outcomes[sim],plot, duration)
                            plot_times, plot_infect = get_plot_data(data_infect, timechange[0], plot_end, t_step, burnin_days = burnin_days, results_type=plot)
                            treatment_data_list.append(plot_infect)

                    treatment_data = np.array(treatment_data_list)
                    median = np.quantile(treatment_data, 0.5, axis=0)
                    if "cumulative" in plot:
                        print_cumul_values(plot, baseline_data, treatment_data, treatment, timechange, CIs_table)

                    ax_infect.plot(plot_times/days_in_year, moving_average(median), color = colours_list_timechange[timechange_i], alpha=1, linewidth = 1, label = legend_treat[treatment_i+1], zorder=100)

                    for CI in CIs_plot:
                        quants = [0.5-0.5*CI, 0.5+0.5*CI]
                        low_quant = np.quantile(treatment_data, quants[0], axis=0)
                        high_quant = np.quantile(treatment_data, quants[1], axis=0)
                        plt.fill_between(plot_times/days_in_year, moving_average(low_quant), moving_average(high_quant), color = colours_list_timechange[timechange_i], alpha=alpha, label='_nolegend_')

                    _timechange = round((timechange[0]-burnin_days)/days_in_year)
                    ax_infect.axvline(x=_timechange,c=colours_list_timechange[timechange_i],ls='--',lw=1,label='_nolegend_') #Plot vertical line
                    ax_infect.text(_timechange, 0.10, " t="+str(_timechange), transform=trans,c=colours_list_timechange[timechange_i])
                    
                # if plot in ["daily_clinical", "daily_relapses"]:    
                #     ax_infect.set_ylim(bottom=0)
                # else:
                #     ax_infect.set_ylim(bottom=0)
                xlims, ylims = getlims(plot)
                if xlims:
                    ax_infect.set_xlim(xlims)
                if ylims:
                    ax_infect.set_ylim(ylims)
                ax_infect.yaxis.grid()
                # ax_infect.set_title(treatments_legend[treatment_i]+scenario)

                leg = fig_infect.legend(legend_timechange, loc="upper left", bbox_to_anchor = (0.13, 0.96))
                for line in leg.get_lines():
                    line.set_linewidth(10.0)

                fig_infect.supxlabel("Time (years)")
                fig_infect.supylabel(y_labels.get(plot, plot))
                    # ax_relapse[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN
                fig_infect.tight_layout()
                fig_infect.savefig("./stored/figs/"+plot+"_"+treatment+"_"+scenario+".png",bbox_inches='tight')

    print(str(count)+" simulations plotted")


#plot histogram of hypnozoite distribution at specified timepoints
hypnozoite_dist = False
if hypnozoite_dist:

    time_plots_years = [0,1,2,10]
    time_plots_days = [int((year+burnin_years)*days_in_year -1*(year==10)) for year in time_plots_years]
    print(time_plots_days)
    subplot_labels = ["a","b","c","d"]
    # assert all(time in params.hyp_snapshot_years for time in time_plots)
    n_subplots = len(time_plots_years)#2#len(params.hyp_snapshot_days)
    max_files = 30
    xlim_lower = 1

    for treatment in treatments:
        fig_hyp_dist, ax_hyp_dist = plt.subplots(math.ceil(n_subplots/2), 2)
        data_hyp_dist_all = []

        count = 0
        
        #load data
        for folder in folders_list:
            filepath = folder+get_filename(treatment, scenarios[0], timechanges[0], duration)
            data_hyp_dist = load_data(filepath, "hypnozoite_dist")
            for sim in range(len(data_hyp_dist)):
                if count+1 > max_files:
                    break
                count += 1
                data_hyp_dist_all.append(data_hyp_dist[sim])

        #plot for each time specified
        for subplot in range(n_subplots):
            time_index = params.hyp_snapshot_days.index(time_plots_days[subplot])

            max_hyp = 0
            #get length info
            for data in data_hyp_dist_all:
                print(time_index)
                data_time = data[time_index][1]
                print(data_time)
                max_hyp = max(max_hyp, len(data_time))

            #combine the lists together to get overall pmf
            np_data = np.zeros(max_hyp)
            for data in data_hyp_dist_all:
                data_time = data[time_index][1]
                #add data from this iteration to existing numpy list
                np_data = np.add(np_data, np.pad(data_time, (0, max_hyp-len(data_time)), mode='constant'))

            #add in grid
            ax_hyp_dist[(int(subplot/2), subplot%2)].set_axisbelow(True)
            ax_hyp_dist[(int(subplot/2), subplot%2)].grid(zorder=0)

            #plot
            ax_hyp_dist[(int(subplot/2), subplot%2)].bar(range(xlim_lower, max_hyp),np_data[xlim_lower:]/sum(np_data))
            ax_hyp_dist[(int(subplot/2), subplot%2)].set_title(subplot_labels[subplot] + ") Distribution at t = "+str(time_plots_years[subplot])+" years")#round(params.hyp_snapshot_days[subplot]/days_in_year)-burnin_years)+" years")
            ax_hyp_dist[(int(subplot/2), subplot%2)].set_xlim(right=30)
            ax_hyp_dist[(int(subplot/2), subplot%2)].set_ylim([0, 0.002])

            #set line for mean
            mean = round(sum([np_data[i]*i for i in range(max_hyp)])/sum(np_data), 3)
            trans = ax_hyp_dist[(int(subplot/2), subplot%2)].get_xaxis_transform()
            ax_hyp_dist[(int(subplot/2), subplot%2)].axvline(x=mean,c='black',ls='--',lw=2,label='_nolegend_') #Plot vertical line
            ax_hyp_dist[(int(subplot/2), subplot%2)].text(mean, 0.5, " mean ="+str(mean), transform=trans,c='black')

            print("Time = "+str(time_plots_years[subplot])+" years")
            print("Number with 0 hypnozoites = "+str(np_data[0]/sum(np_data)))
            print("Max hypnozoites = "+str(max_hyp-1))
            # print("Median = "+str(np.quantile(np_data/sum(np_data), 0.5)))
            print("Median = "+str((np_data/sum(np_data)).cumsum().searchsorted(0.5)))
            print("Mean = "+str(sum([np_data[i]*i for i in range(max_hyp)])/sum(np_data)))
            print("********************")
            
            #ax_hyp_dist[(subplot%2, int(subplot/2))]

            # count += 1
        fig_hyp_dist.suptitle(treatment)
        fig_hyp_dist.supxlabel("Number of hypnozoites")
        fig_hyp_dist.supylabel("Probability")
        fig_hyp_dist.tight_layout()
        fig_hyp_dist.savefig("./stored/figs/hyp_dist"+treatment+scenarios[0]+".png",bbox_inches='tight')

    print("filecount = "+str(count))

hypnozoite_dist_mean_all = True
if hypnozoite_dist_mean_all:


    time_plots = params.hyp_snapshot_days
    # time_plots = [params.hyp_snapshot_days[0]]
    n_times = len(time_plots)#2#len(params.hyp_snapshot_days)
    max_files = 30
    data_baseline_hyp = []
    data_treatment_hyp = [[[] for _ in range(len(timechanges))] for _ in range(len(timechanges))]
    distance_between = 0.15
    width = 2
    markerwidth = 4
    linetypes = {"Baseline": '-o', treatments[0]: '-^', treatments[1]: '-s'}
    plot_types = ["hypnozoites_mean_all", "hypnozoites_mean_nonzero", "hypnozoites_prevalence_pc"]

    count = 0
    
    #load data
    for folder in folders_list:
        filepath = folder+get_filename(*baseline)
        data_baseline = load_data(filepath, "hypnozoite_dist")

        for sim in range(len(data_baseline)):
            if count+1 > max_files:
                break
            count += 1
            data_baseline_hyp.append(data_baseline[sim])
        
        for treatment_i in range(len(treatments)):
            for timechange_i in range(len(timechanges)):
                filepath = folder+get_filename(treatments[treatment_i], scenarios[0],timechanges[timechange_i],duration)
                data_treat = load_data(filepath, "hypnozoite_dist")
                for sim in range(len(data_treat)):
                    data_treatment_hyp[treatment_i][timechange_i].append(data_treat[sim])
                # print(sum(data_treatment_hyp[treatment_i][timechange_i][0][0][1]))

    for plot in plot_types:
        fig_hyp_dist, ax_hyp_dist = plt.subplots(1)
        trans = ax_hyp_dist.get_xaxis_transform() #Transformation for labelling line

        #plot for each time specified
        means = []
        for i in range(n_times):
            time_index = params.hyp_snapshot_days.index(int(time_plots[i]))

            #baseline plot
            max_hyp = 0
            #get length info
            for data in data_baseline_hyp:
                data_time = data[time_index][1]
                max_hyp = max(max_hyp, len(data_time))
            # print("Baseline max hypnozoites at t = "+str((round(time_plots[i]/days_in_year - burnin_years)))+" = "+str(max_hyp))

            #combine the lists together to get overall pmf
            np_data = np.zeros(max_hyp)
            for data in data_baseline_hyp:
                data_time = data[time_index][1]
                #add data from this iteration to existing numpy list
                np_data = np.add(np_data, np.pad(data_time, (0, max_hyp-len(data_time)), mode='constant'))

            #plot
            raw_data = []#[hyp for _ in range(np_data[hyp])] for hyp in range(max_hyp)]

            if "hypnozoites_mean_all" in plot:
                i0 = 0
            elif plot == "hypnozoites_mean_nonzero":
                i0 = 1
            elif plot == "hypnozoites_prevalence_pc":
                i0 = 1

            if "mean" in plot:
                for nhyp in range(i0, max_hyp):
                    for _ in range(int(np_data[nhyp])):
                        raw_data.append(nhyp)
                if "per_thousand" in plot:
                    means.append(100*np.mean(raw_data))
                else: 
                    means.append(np.mean(raw_data))
                    
            elif "prevalence" in plot:
                means.append(100*sum(np_data[i0:])/sum(np_data))

        print("Baseline")
        print(means[-1])

        h_baseline = ax_hyp_dist.plot([(time_plots[i]/days_in_year)-burnin_years-distance_between for i in range(n_times)], means, linetypes["Baseline"], color = colours_list_timechange[-1], linewidth = width, markeredgewidth = markerwidth)
            # vp_baseline = ax_hyp_dist.violinplot(raw_data, positions = [(time_plots[i]/days_in_year)-burnin_years-distance_between], showmeans=True, widths = [width])
            # for pc in vp_baseline["bodies"]:
            #     pc.set_facecolor(colours_map_treatment["Baseline"])
            #     pc.set_edgecolor(colours_map_treatment["Baseline"])
            
            # for part in ["cmeans", "cmins", "cmaxes", "cbars"]:
            #     vp_baseline[part].set_color(colours_map_treatment["Baseline"])

            #plot for treatments

            # vp_treat=[0,0]
        h_treatment=[0]*len(treatments)
        for treatment_i in range(len(treatments)):
            scenario = treatments[treatment_i]
            for timechange_i in range(len(timechanges)):
                timechange = timechanges[timechange_i]

                times = []
                means = []
                for i in range(n_times):
                    time_index = params.hyp_snapshot_days.index(int(time_plots[i]))

                    if time_plots[i] >= timechange[0]:
                        max_hyp = 0
                        #get length info
                        for data in data_treatment_hyp[treatment_i][timechange_i]:
                            data_time = data[time_index][1]
                            max_hyp = max(max_hyp, len(data_time))

                        #combine the lists together to get overall pmf
                        np_data = np.zeros(max_hyp)
                        for data in data_treatment_hyp[treatment_i][timechange_i]:
                            data_time = data[time_index][1]
                            #add data from this iteration to existing numpy list
                            np_data = np.add(np_data, np.pad(data_time, (0, max_hyp-len(data_time)), mode='constant'))

                        #get means
                        raw_data = []#[hyp for _ in range(np_data[hyp])] for hyp in range(max_hyp)]

                        # if plot == "hypnozoites_mean_all":
                        #     i0 = 0
                        # elif plot == "hypnozoites_mean_nonzero":
                        #     i0 = 1
                        # elif plot == "hypnozoites_prevalence_pc":
                        #     i0 = 1

                        if "mean" in plot:
                            for nhyp in range(i0, max_hyp):
                                for _ in range(int(np_data[nhyp])):
                                    raw_data.append(nhyp)#-distance_between*.2*(treatment_i+1))   
                            if "per_thousand" in plot:
                                means.append(100*np.mean(raw_data))
                            else: 
                                means.append(np.mean(raw_data))
                        elif "prevalence" in plot:
                            means.append(100*sum(np_data[i0:])/sum(np_data))

                        times.append((time_plots[i]/days_in_year)-burnin_years+treatment_i*distance_between)

                        # if i==0:
                        #     print("Times: ")
                        #     print(times)
                        #     print("Means: ")
                        #     print(means)

                # vp_treat[treatment_i] = 
                print("\n"+scenario+str(timechange))
                print(means[-1])
                h_treatment[treatment_i] = ax_hyp_dist.plot(times, means, linetypes[scenario], color = colours_list_timechange[timechange_i], linewidth=width, markeredgewidth = markerwidth)
    
        ax_hyp_dist.set_ylim(bottom=0)

        # #Create legend
        # bl_handle = [mpatches.Patch(color=colours[-1])]
        
        for timechange_i in range(len(timechanges_year)):
            timechange = timechanges_year[timechange_i]
            # _timechange = round((timechange[0]-burnin_days)/days_in_year)
            ax_hyp_dist.axvline(x=timechange,c=colours_list_timechange[timechange_i],ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_hyp_dist.text(timechange, 0.05, " t="+str(timechange), transform=trans,c=colours_list_timechange[timechange_i])

        bl_handle = [plt.scatter([],[],marker = linetypes[scenario][1], color='k',s=80) for scenario in ["Baseline"]]
        scenario_handle = [plt.scatter([],[],marker = linetypes[scenario][1], color='k',facecolors='none',s=80) for scenario in treatments]
        time_handle = [mpatches.Patch(color=colours_list_timechange[timechange_i]) for timechange_i in range(len(timechanges))]

        fig_hyp_dist.legend(bl_handle+scenario_handle+time_handle, legend_treat+legend_timechange[1:], loc="upper right", bbox_to_anchor = (0.97, 0.97),framealpha=0.5)
        # fig_hyp_dist.legend([bl_handle, treat_handle[0], treat_handle[1]], timechange, loc="upper right", bbox_to_anchor = (0.97, 0.97),framealpha=1)

            # print("Time = "+str(time_plots[subplot])+" years")
            # print("Number with 0 hypnozoites = "+str(np_data[0]/sum(np_data)))
            # print("Max hypnozoites = "+str(max_hyp-1))
            # # print("Median = "+str(np.quantile(np_data/sum(np_data), 0.5)))
            # print("Median = "+str((np_data/sum(np_data)).cumsum().searchsorted(0.5)))
            # print("Mean = "+str(sum([np_data[i]*i for i in range(max_hyp)])/sum(np_data)))
            # print("********************")
            
            #ax_hyp_dist[(subplot%2, int(subplot/2))]

            # count += 1
        # fig_hyp_dist.suptitle("Hypnozoite distribution before and after a policy change [low-dose PQ]")
        fig_hyp_dist.supxlabel("Time (years)")
        # fig_hyp_dist.supylabel("Mean number of hypnozoites")
        fig_hyp_dist.supylabel(y_labels.get(plot, plot))
        
        xticks = [round(time/days_in_year-burnin_years) for time in time_plots]
        ax_hyp_dist.set_xticks(xticks)
        fig_hyp_dist.tight_layout()
        fig_hyp_dist.savefig("./stored/figs/"+plot+".png",bbox_inches='tight')

    print("filecount = "+str(count))


plt.show(block=False)
print("Plots finished. Press ENTER to continue.")
input()
plt.close()

