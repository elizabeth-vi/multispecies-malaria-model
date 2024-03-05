import matplotlib.pyplot as plt
import json
import numpy as np
from model_params import model_params
import itertools
from index_names import Compartments, Mozzie_labels
import os


def plot_all_human_pv_compartments(treatment, timechange, data, t_start, t_end, t_step):
    #Set up fig
    human_fig, human_axs = plt.subplots(len(Compartments))
    human_fig.supxlabel('time (days)')
    human_fig.suptitle('Human p.v. Compartments for '+treatment+" introduced at time "+timechange)

    # variables = [*data.keys()]
    
    #Plot by compartments
    # pv_hist = variables[1]
    # vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = data[t_start : t_start+n_vals]

    for compartment in Compartments:
        pv_compartment_data = [val[compartment.value] for val in vals_pv]

        #Plot human p.v. infections
        human_axs[compartment.value].plot(times,pv_compartment_data)
        human_axs[compartment.value].set_title(compartment.name,x=1.03,y=0)
        human_axs[compartment.value].grid()
        human_axs[compartment.value].set_ylim(bottom=0)

def plot_human_select_pv_compartments_together(pv_data, t_start, t_end, t_step, compartments_str, fig=None, colour='b'):
    #Set up fig
    fig = plt.figure(fig)

    # variables = [*data.keys()]
    compartments = [Compartments[comp] for comp in compartments_str]
    
    #Plot by compartments
    # pv_hist = variables[1]
    # vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = pv_data[t_start : t_start+n_vals]

    for index in range(len(compartments)):
        compartment = compartments[index]
        #compartment = Compartments(vars(compartments[index]))
        pv_compartment_data = [val[compartment.value] for val in vals_pv]

        #Plot human p.v. infections
        plt.plot(times,pv_compartment_data)#,colour)
    
    legend = [compartment.name for compartment in compartments]
    plt.xlabel('time (days)')
    plt.title("Human p.v. Compartments "+str(legend))#+" for "+treatment+" introduced at time "+timechange)
    plt.ylim(bottom=0)
    #plt.legend(legend)
    plt.grid()

    return fig

def plot_human_pv_infections(data, t_start, t_end, t_step, fig=None):
    
    #Compartments you want to sum and graph as "infected"
    infect_comps_str = ["I","A"] 
    infect_comps = [Compartments[comp] for comp in infect_comps_str]
    
    #Set up fig
    fig = plt.figure(fig)

    #Get the pv data
    variables = [*data.keys()]
    pv_hist = variables[1] #human_pop_pv_history
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)

    vals_pv = vals_pv[t_start : t_start+n_vals] #All compartments data
    pv_compartment_data = [(sum(val[comp.value] for comp in infect_comps)) for val in vals_pv] #Summed values

    #Plot human p.v. infections
    plt.plot(times,pv_compartment_data)
    plt.xlabel('time (days)')
    plt.ylabel('Infections')
    plt.title('Human p.v. infections: sum(' + str(infect_comps_str))# + ") for "+treatment+" introduced at time "+timechange)
    plt.ylim(bottom=0) #set lower bound for plotting
    plt.grid()
    
    return fig

#WORK IN PROGRESS. Not particularly informative
def plot_all_tracked_variables(data, t_start, t_end, t_step):

    #variables = [*data.keys()] #For all variables including pf
    variables = ["human_pop_pv_history","mozzie_pv_infectious","mozzie_pop_history"]#,"pv_outcomes","num_TGD"] #Only plot PV-related
    for variable in variables:

        vals = data.get(variable)[0]
        times = np.arange(start=t_start, stop=t_end, step=t_step)
        n_vals = len(times)

        plt.figure()
        plt.grid(True)

        if isinstance(vals[0],list):
            try:
                for i in range(len(vals[0])):
                    yvals = [val[i] for val in vals]

                    ## LINE GRAPH ##
                    plt.plot(times,yvals[t_start : t_start+n_vals])
                    
                #print(Compartments._member_names_)  
                if variable == "human_pop_pv_history":
                    plt.legend(Compartments._member_names_)
                elif variable == "mozzie_pop_history":
                    plt.legend(Mozzie_labels._member_names_) 
            except Exception as e:
                print("Breaking, variable = "+variable)
                print(repr(e))
                plt.close()
                break
        else:
            yvals = vals

            ## LINE GRAPH ##
            plt.plot(times,yvals[t_start : t_start+n_vals])
        plt.xlabel('time')
        plt.ylabel(variable)
        plt.title(variable)

def get_filepath(treatment, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    #duration = str(int(params.time_day_end)) #Simulation duration. Used in filename
    filepath_loc = "./stored/results_variables/duration_" + str(duration) + "_48632" + "/"
    filepath = filepath_loc + treatment + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filepath

def get_filename(treatment, timechange, duration):
    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    filename = treatment + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
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

def load_data(filepath, category):
    #Load and return data
    try:
        data_all = json.load(open(filepath, 'r'))
        data_cat = data_all.get(category)[0]
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

#*********************************SCRIPT START*********************************


#Filepath name info
#[filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]


params = model_params()
treatments = ["PLD","Taf"] #Must match file names. Options: ["PLD","PHD","Taf"]]
treatments_titles = ["Change to Greater Rates of low-dose Primaquine Treatment","Change to Greater Rates of Tafenoquine Treatment"] #["Change to Policy 1", "Change to Policy 2", "Change to Policy 3"]
treatments_legend = ["Low-dose Primaquine","Tafenoquine"] #["Change to Policy 1", "Change to Policy 2", "Change to Policy 3"]
colours = ['r','b','g', 'k'] #for plotting
colours_map_treatment = {"Baseline":"black", "PLD":"red", "Taf":"blue"}
timechanges = [[1826],[3287],[4748]] #[[0], [1461], [2922]] #which time changes you want to plot, in days. e.g.: [0, 730, 1461]
duration = int(params.time_day_end)-1 #Simulation duration
days_in_year = 365.25
burnin_years = 5 #years
burnin_days = burnin_years * days_in_year

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start+burnin_days, duration, params.time_day_step] #Time of start/end you want to plot, and timestep used
plt.rcParams["figure.figsize"] = (11,8.5)
# plt.tight_layout()

#Set up figures
# fig_infect = plt.figure()
# fig_IAL = plt.figure()
# fig_TG = plt.figure()
legend_treat = ["Baseline"]
legend_treat.extend(treatments_legend)
legend_timechange = ["Baseline"]
legend_timechange.extend(["Change policy at t = "+str(round(x[0]-burnin_days))+" days" for x in timechanges])

#baseline plotting
# if (treatment := "baseline"):
#     timechange = duration
#     data = get_file_data("baseline", timechange, duration)
#     #filepath = filepath_loc + treatment + filepath_timechange + timechange + filepath_duration + duration + filepath_type
#     #data = json.load(open(filepath, 'r'))
#     legend.append("Baseline")

#     fig_infect, ax_infect = plt.subplots(len(treatments))
#     #Figures to plot
#     baseline_infect = plot_human_pv_infections(data, plot_start, plot_end, t_step)
#     # fig_IAL = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["A"], fig=fig_IAL, colour = colours[-1])
#     # fig_TG = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["G"], fig=fig_TG, colour = colours[-1])

#Set up fig 1
#Infections for all treatment options
plot_infections_1 = False
if plot_infections_1:
    #set up fig
    fig_infect, ax_infect = plt.subplots(len(treatments))
    # fig_infect.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_infect.suptitle("Clinical infections (I) when treatment policy is changed at different times")
    #plot
    for treatment_i in range(len(treatments)):
        #baseline plotting
        baseline_filepath = get_filepath("Baseline",[],duration)
        compartment_data_baseline = load_data(baseline_filepath, "human_pop_pv_history")
        infections_baseline = [val[1] for val in compartment_data_baseline]
        times_baseline, data_baseline_I = get_plot_data(infections_baseline, plot_start, plot_end, t_step)
        # plt.plot(times_baseline,data_baseline)
        # baseline_infect = plot_human_pv_infections(data_baseline, plot_start, plot_end, t_step, ax_infect[treatment_i])

        # print(data_baseline_I)

        ax_infect[treatment_i].plot(times_baseline,data_baseline_I,colours[-1],zorder=4)

        #fig_infect[treatment_i] = plt.figure(baseline_infect)
        # fig_infect[treatment_i].legend(legend,loc=7)
        # print(treatment_i)

        #Transformation for labelling line
        trans = ax_infect[treatment_i].get_xaxis_transform()
        

        
        #Plot infections
        for timechange_i in range(len(timechanges)):
            colour = colours[timechange_i]
            treatment = treatments[treatment_i]
            timechange = timechanges[timechange_i]
            filepath = get_filepath(treatment,timechange,duration)
            data = load_data(filepath, "human_pop_pv_history")
            infections = [val[1] for val in data]
            times, data_I = get_plot_data(infections, plot_start, plot_end, t_step)
            
            #Plot
            #For each treatment, plot I for each timechange option
            ax_infect[treatment_i].plot(times,data_I,colour)
            ax_infect[treatment_i].axvline(x=timechange[0],c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_infect[treatment_i].text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colour)
        
        ax_infect[treatment_i].set_ylim(bottom=0)
        ax_infect[treatment_i].set_title(treatments_titles[treatment_i])

    fig_infect.legend(legend_timechange)
    fig_infect.supxlabel("Time (days)")
    fig_infect.supylabel("Clinical Infections (I)")
        # ax_infect[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN

    fig_infect.tight_layout()
    fig_infect.savefig("./stored/figs/clinical_infections.png",bbox_inches='tight')
        # plot_all_human_pv_compartments(treatment, "0", get_compartment_data(treatment,0,duration), plot_start, plot_end, t_step)

plot_clinical_2 = False
if plot_clinical_2:
    fig_infect, ax_infect = plt.subplots(1)
    alpha = 0.3
    # fig_relapse.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_relapse.suptitle("Relapses when treatment policy is changed at different times")
    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:
        print(folder)
        
    # for treatment_i in range(len(treatments)):
        #baseline plotting
        filepath = folder+get_filename("Baseline",[],duration)
        data_pv_outcomes = load_data(filepath, "pv_outcomes")
        # baseline_filepath = get_filepath("Baseline",[],duration)
        # relapse_data_baseline = load_data(baseline_filepath, "pv_recorded_relapses")
        #slice data based on plot start/end
        data_infect = [data_pv_outcomes[1].count(i) for i in range(duration)]
        # data_cumul_relapse = np.cumsum(data_relapse)
        plot_times, plot_infect = get_plot_data(data_infect, plot_start, plot_end, t_step, burnin_days = burnin_days)
        
        # times_baseline, relapse_baseline = get_plot_data(relapse_data_baseline, plot_start, plot_end, t_step)
        # ax_relapse[treatment_i].plot(times_baseline,relapse_baseline,colours[-1],zorder=4)
        ax_infect.plot(plot_times, plot_infect, color = colours_map_treatment["Baseline"], alpha=alpha, zorder = 4)


        #Transformation for labelling line
        trans = ax_infect.get_xaxis_transform()

        
        #Plot infections

        for timechange_i in range(len(timechanges)):
            timechange = timechanges[timechange_i]
            treatment_i = 1
            if True:
            # for treatment_i in range(len(treatments)):
                treatment = treatments[treatment_i]
                filepath = folder + get_filename(treatment,timechange,duration)
                data_pv_outcomes = load_data(filepath, "pv_outcomes")
                data_infect = [data_pv_outcomes[1].count(i) for i in range(duration)]
                plot_times, plot_infect = get_plot_data(data_infect, plot_start, plot_end, t_step, burnin_days = burnin_days)

                ax_infect.plot(plot_times,plot_infect, color = colours[timechange_i], alpha=alpha, label = legend_treat[treatment_i+1])

            _timechange = round(timechange[0]-burnin_days)
            ax_infect.axvline(x=_timechange,c=colours[timechange_i],ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_infect.text(_timechange, 0.15, " t="+str(_timechange), transform=trans,c=colours[timechange_i])
        
    ax_infect.set_ylim(bottom=0)
    ax_infect.set_title("(Fig 1a) Daily clinical infections for baseline vs tafenoquine treatment changes at different times")
    fig_infect.legend(legend_timechange)
    fig_infect.supxlabel("Time (days)")
    fig_infect.supylabel("Daily clinical infections")
        # ax_relapse[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN
    fig_infect.tight_layout()
    fig_infect.savefig("./stored/figs/clinical_daily.png",bbox_inches='tight')

#plot relapses
plot_relapses_1 = False
if plot_relapses_1:
    fig_relapse, ax_relapse = plt.subplots(len(treatments))
    # fig_relapse.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_relapse.suptitle("Relapses when treatment policy is changed at different times")
    
    for treatment_i in range(len(treatments)):
        #baseline plotting
        baseline_filepath = get_filepath("Baseline",[],duration)
        relapse_data_baseline = load_data(baseline_filepath, "pv_recorded_relapses")
        #slice data based on plot start/end
        times_baseline, relapse_baseline = get_plot_data(relapse_data_baseline, plot_start, plot_end, t_step)
        ax_relapse[treatment_i].plot(times_baseline,relapse_baseline,colours[-1],zorder=4)

        #Transformation for labelling line
        trans = ax_relapse[treatment_i].get_xaxis_transform()

        
        #Plot infections
        for timechange_i in range(len(timechanges)):
            colour = colours[timechange_i]
            treatment = treatments[treatment_i]
            timechange = timechanges[timechange_i]
            filepath = get_filepath(treatment,timechange,duration)
            data = load_data(filepath, "pv_recorded_relapses")
            times, data_plot = get_plot_data(data, plot_start, plot_end, t_step)
            
            #Plot
            #For each treatment, plot I for each timechange option
            ax_relapse[treatment_i].plot(times,data_plot,colour)
            ax_relapse[treatment_i].axvline(x=timechange[0],c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_relapse[treatment_i].text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colour)

            # data = load_data(filepath, "pv_recorded_relapses")
            # times, data_plot = get_plot_data(data, plot_start, plot_end, t_step)
            # ax_relapse[treatment_i].plot(times,data_plot,colour, marker='x')

        
        ax_relapse[treatment_i].set_ylim(bottom=0)
        ax_relapse[treatment_i].set_title(treatments_titles[treatment_i])
    fig_relapse.legend(legend_timechange)
    fig_relapse.supxlabel("Time (days)")
    fig_relapse.supylabel("daily recorded relapses")
        # ax_relapse[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN
    fig_relapse.tight_layout()
    fig_relapse.savefig("./stored/figs/recorded_relapses.png",bbox_inches='tight')

plot_relapses_2 = False
if plot_relapses_2:
    fig_relapse, ax_relapse = plt.subplots(1)
    alpha = 0.3
    # fig_relapse.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_relapse.suptitle("Relapses when treatment policy is changed at different times")
    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:
        print(folder)
        
    # for treatment_i in range(len(treatments)):
        #baseline plotting
        filepath = folder+get_filename("Baseline",[],duration)
        # data_pv_outcomes = load_data(filepath, "pv_outcomes")
        data_relapse = load_data(filepath, "pv_actual_relapses")
        # baseline_filepath = get_filepath("Baseline",[],duration)
        # relapse_data_baseline = load_data(baseline_filepath, "pv_recorded_relapses")
        #slice data based on plot start/end
        # data_relapse = [data_pv_outcomes[2].count(i) for i in range(duration)]
        # data_cumul_relapse = np.cumsum(data_relapse)
        plot_times, plot_relapse = get_plot_data(data_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)
        
        # times_baseline, relapse_baseline = get_plot_data(relapse_data_baseline, plot_start, plot_end, t_step)
        # ax_relapse[treatment_i].plot(times_baseline,relapse_baseline,colours[-1],zorder=4)
        ax_relapse.plot(plot_times, plot_relapse, color = colours_map_treatment["Baseline"], alpha=alpha, zorder = 4)


        #Transformation for labelling line
        trans = ax_relapse.get_xaxis_transform()

        
        #Plot infections

        for timechange_i in range(len(timechanges)):
            timechange = timechanges[timechange_i]
            treatment_i = 1
            if True:
            # for treatment_i in range(len(treatments)):
                treatment = treatments[treatment_i]
                filepath = folder + get_filename(treatment,timechange,duration)
                # data_pv_outcomes = load_data(filepath, "pv_outcomes")
                data_relapse = load_data(filepath, "pv_actual_relapses")
                # data_relapse = [data_pv_outcomes[2].count(i) for i in range(duration)]
                plot_times, plot_relapse = get_plot_data(data_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)

                ax_relapse.plot(plot_times,plot_relapse, color = colours[timechange_i], alpha=alpha, label = legend_treat[treatment_i+1])

            _timechange = round(timechange[0]-burnin_days)
            ax_relapse.axvline(x=_timechange,c=colours[timechange_i],ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_relapse.text(_timechange, 0.15, " t="+str(_timechange), transform=trans,c=colours[timechange_i])
        
    ax_relapse.set_ylim(bottom=0)
    ax_relapse.set_title("(Fig 2a) Daily relapses for baseline vs tafenoquine treatment changes at different times")
    fig_relapse.legend(legend_timechange, loc = "upper right")
    fig_relapse.supxlabel("Time (days)")
    fig_relapse.supylabel("daily relapses")
        # ax_relapse[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN
    fig_relapse.tight_layout()
    fig_relapse.savefig("./stored/figs/daily_relapses.png",bbox_inches='tight')

plot_IA_TG = False
if plot_IA_TG:
    #setup IA fig
    fig_IA, ax_IA = plt.subplots(len(treatments), figsize = (11,8.5))
    # fig_IA.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_IA.suptitle("Blood-stage infections (I+A) when treatment policy is changed at different times")
    #setup TG fig
    fig_TG, ax_TG = plt.subplots(len(treatments), figsize = (11,8.5))
    # fig_TG.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_TG.suptitle("People undergoing treatment (T+G) when treatment policy is changed at different times")
    #setup TG new treatments fig
    fig_treat, ax_treat = plt.subplots(len(treatments), figsize = (11,8.5))
    # fig_treat.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat.suptitle("New treatments (T+G) each day")

    fig_treat_T, ax_treat_T = plt.subplots(len(treatments), figsize = (11,8.5))
    # fig_treat_T.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat_T.suptitle("New treatments (T) each day")
    fig_treat_G, ax_treat_G = plt.subplots(len(treatments), figsize = (11,8.5))
    # fig_treat_G.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat_G.suptitle("New treatments (G) each day")

    for treatment_i in range(len(treatments)):
        #baseline plotting
        baseline_filepath = get_filepath("Baseline",[],duration)
        compartment_data_baseline = load_data(baseline_filepath, "human_pop_pv_history")
        
        #IA
        fig1_comps = [2]
        IA_baseline = [sum([val[comp] for comp in fig1_comps]) for val in compartment_data_baseline]
        times_baseline, data_baseline_IA = get_plot_data(IA_baseline, plot_start, plot_end, t_step)
        ax_IA[treatment_i].plot(times_baseline,data_baseline_IA,colours[-1],zorder=4)

        #TG
        fig2_comps = [5,6]
        TG_baseline = [sum([val[comp] for comp in fig2_comps]) for val in compartment_data_baseline]
        times_baseline, data_baseline_TG = get_plot_data(TG_baseline, plot_start, plot_end, t_step)
        ax_TG[treatment_i].plot(times_baseline,data_baseline_TG,colours[-1],zorder=4)

        #TG
        fig3_comps = [3,4]
        compartment_data_baseline_3 = load_data(baseline_filepath, "pv_outcomes")

        treat_baseline_T = [sum([compartment_data_baseline_3[comp].count(i) for comp in [3]]) for i in range(duration)]
        times_baseline, data_baseline_treat_T = get_plot_data(treat_baseline_T, plot_start, plot_end, t_step)
        ax_treat_T[treatment_i].plot(times_baseline,data_baseline_treat_T,colours[-1],zorder=4)

        treat_baseline_G = [sum([compartment_data_baseline_3[comp].count(i) for comp in [4]]) for i in range(duration)]
        times_baseline, data_baseline_treat_G = get_plot_data(treat_baseline_G, plot_start, plot_end, t_step)
        ax_treat_G[treatment_i].plot(times_baseline,data_baseline_treat_G,colours[-1],zorder=4)

        # treat_baseline = [sum([compartment_data_baseline_3[comp].count(i) for comp in fig3_comps]) for i in range(duration)]
        treat_baseline = [treat_baseline_G[i] + treat_baseline_T[i] for i in range(len(treat_baseline_G))]
        times_baseline, data_baseline_treat = get_plot_data(treat_baseline, plot_start, plot_end, t_step)
        ax_treat[treatment_i].plot(times_baseline,data_baseline_treat,colours[-1],zorder=4)

        #Transformation for labelling line
        trans = ax_treat[treatment_i].get_xaxis_transform()

        #Plot
        for timechange_i in range(len(timechanges)):
            colour = colours[timechange_i]
            treatment = treatments[treatment_i]
            timechange = timechanges[timechange_i]
            filepath = get_filepath(treatment,timechange,duration)
            data = load_data(filepath, "human_pop_pv_history")
            IA = [sum([val[comp] for comp in fig1_comps]) for val in data]
            times, data_IA = get_plot_data(IA, plot_start, plot_end, t_step)
            TG = [sum([val[comp] for comp in fig2_comps]) for val in data]
            times, data_TG = get_plot_data(TG, plot_start, plot_end, t_step)
            
            data_2 = load_data(filepath, "pv_outcomes")
            treat_all = [sum([data_2[comp].count(i) for comp in fig3_comps]) for i in range(duration)]
            times, data_treat_all = get_plot_data(treat_all, plot_start, plot_end, t_step)
            treat_T = [sum([data_2[comp].count(i) for comp in [3]]) for i in range(duration)]
            times, data_treat_T = get_plot_data(treat_T, plot_start, plot_end, t_step)
            treat_G = [sum([data_2[comp].count(i) for comp in [4]]) for i in range(duration)]
            times, data_treat_G = get_plot_data(treat_G, plot_start, plot_end, t_step)

            #Plot
            #For each treatment, plot IA for each timechange option
            ax_IA[treatment_i].plot(times,data_IA,colour)
            ax_IA[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_IA[treatment_i].text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colour)
            #For each treatment, plot TG for each timechange option
            ax_TG[treatment_i].plot(times,data_TG,colour)
            ax_TG[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_TG[treatment_i].text(timechange[0], 0.05, " t="+str(timechange[0]), transform=trans,c=colour)
            #For each treatment, plot new treatments for each timechange option
            ax_treat[treatment_i].plot(times,data_treat_all,colour)
            ax_treat[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_treat[treatment_i].text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colour)
            #For each treatment, plot new treatments T for each timechange option
            ax_treat_T[treatment_i].plot(times,data_treat_T,colour)
            ax_treat_T[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_treat_T[treatment_i].text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colour)
            #For each treatment, plot new treatments G for each timechange option
            ax_treat_G[treatment_i].plot(times,data_treat_G,colour)
            ax_treat_G[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
            ax_treat_G[treatment_i].text(timechange[0]+0.05, 0.83, " t="+str(timechange[0]), transform=trans,c=colour)
        
        ax_IA[treatment_i].set_ylim(bottom=0)
        ax_IA[treatment_i].set_title(treatments_titles[treatment_i]) #adjust to be exact pN
        ax_TG[treatment_i].set_ylim((0,2100))
        ax_TG[treatment_i].set_title(treatments_titles[treatment_i]) #adjust to be exact pN
        ax_treat[treatment_i].set_title(treatments_titles[treatment_i]) #adjust to be exact pN
        ax_treat[treatment_i].set_ylim((0,200))
        ax_treat_T[treatment_i].set_title(treatments_titles[treatment_i]) #adjust to be exact pN
        ax_treat_T[treatment_i].set_ylim((0,150))
        ax_treat_G[treatment_i].set_title(treatments_titles[treatment_i]) #adjust to be exact pN
        ax_treat_G[treatment_i].set_ylim((0,150))
    
    
    fig_IA.legend(legend_timechange)
    fig_IA.supxlabel("Time (days)")
    fig_IA.supylabel("Asymptomatic Infections (A)")
    # ax_IA[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN
    # plt.tight_layout()
    # fig_IA.savefig("./stored/figs/IA.png",bbox_inches='tight')
    fig_IA.tight_layout()
    fig_IA.savefig("./stored/figs/IA.png", bbox_inches='tight')
    
    fig_TG.legend(legend_timechange)
    fig_TG.supxlabel("Time (days)")
    fig_TG.supylabel("People in treatment (T+G)")
    # plt.tight_layout()
    # fig_TG.savefig("./stored/figs/TG.png",bbox_inches='tight')
    fig_TG.tight_layout()
    fig_TG.savefig("./stored/figs/TG.png", bbox_inches='tight')

    
    fig_treat.legend(legend_timechange)
    fig_treat.supxlabel("Time (days)")
    fig_treat.supylabel("New treatments (per day)")
    # plt.tight_layout()
    # fig_treat.savefig("./stored/figs/new_treatments.png",bbox_inches='tight')
    fig_treat.tight_layout()
    fig_treat.savefig("./stored/figs/new_treatments.png", bbox_inches='tight')

    
    fig_treat_T.legend(legend_timechange)
    fig_treat_T.supxlabel("Time (days)")
    fig_treat_T.supylabel("New blood-stage treatments (per day)")
    # plt.tight_layout()
    # fig_treat_T.savefig("./stored/figs/new_T.png",bbox_inches='tight')
    fig_treat_T.tight_layout()
    fig_treat_T.savefig("./stored/figs/new_T.png", bbox_inches='tight')

    
    fig_treat_G.legend(legend_timechange)
    fig_treat_G.supxlabel("Time (days)")
    fig_treat_G.supylabel("New radical cure treatments (per day)")

    # # fig_IA.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_TG.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat_T.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_treat_G.tight_layout(rect=[0, 0.03, 1, 0.95])


    fig_treat_G.tight_layout()
    fig_treat_G.savefig("./stored/figs/new_G.png", bbox_inches='tight')


#plot relapses
plot_relapses_G6PD_1 = False
G6PD_statuses = ["Ineligible", "0-30%", "30-70%","70-100%"]
if plot_relapses_G6PD_1:
    fig_relapseG6PD, ax_relapseG6PD = plt.subplots()
    # fig_relapseG6PD.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig_relapseG6PD.suptitle("Relapses by G6PD status")
    G6PD_i = 1
    timechange = timechanges[1]
    #baseline plotting
    # filepath = get_filepath("Baseline","[]",duration)
    filepath = get_filepath(treatments[-1],timechange,duration)
    relapse_data = load_data(filepath, "G6PD_hypnozoites")
    num_G6PD = load_data(filepath, "num_G6PD")
    # relapse_G6PD = [relapse_data[i][G6PD_i]/num_G6PD[G6PD_i] for i in range(len(relapse_data))]
    # times_G6PD, relapse_G6PD = get_plot_data(relapse_G6PD, plot_start, plot_end, t_step)
    # ax_relapseG6PD.plot(times_G6PD,relapse_G6PD,colours[-1],zorder=4)

    trans = ax_relapseG6PD.get_xaxis_transform()



    for G6PD_i in range(len(G6PD_statuses)):
        relapse_G6PD = [relapse_data[i][G6PD_i]/num_G6PD[G6PD_i] for i in range(len(relapse_data))]
        # print(relapse_G6PD)
        #slice data based on plot start/end
        times_G6PD, relapse_G6PD = get_plot_data(relapse_G6PD, plot_start, plot_end, t_step)
        ax_relapseG6PD.plot(times_G6PD,relapse_G6PD,colours[G6PD_i],zorder=4)

    # for treatment_i in range(len(treatments)):
    #     filepath = get_filepath(treatments[treatment_i],"[0]",duration)
    #     relapse_data = load_data(filepath, "G6PD_hypnozoites")
    #     num_G6PD = load_data(filepath, "num_G6PD")
    #     relapse_G6PD = [relapse_data[i][G6PD_i]/num_G6PD[G6PD_i] for i in range(len(relapse_data))]
    #     # print(relapse_G6PD)
    #     #slice data based on plot start/end
    #     times_G6PD, relapse_G6PD = get_plot_data(relapse_G6PD, plot_start, plot_end, t_step)
    #     ax_relapseG6PD.plot(times_G6PD,relapse_G6PD,colours[treatment_i],zorder=4)
        
    ax_relapseG6PD.axvline(x=timechange[0],c=colours[2],ls='--',lw=1,label='_nolegend_') #Plot vertical line
    ax_relapseG6PD.text(timechange[0], 0.15, " t="+str(timechange[0]), transform=trans,c=colours[2])
    ax_relapseG6PD.set_ylim(bottom=0)
    ax_relapseG6PD.legend(G6PD_statuses)
    # ax_relapseG6PD.legend(["Baseline"] + treatments_titles)
    ax_relapseG6PD.set_ylabel("proportion with hypnozoites")
    ax_relapseG6PD.set_xlabel("Time (days)")
    # ax_relapseG6PD.set_title("Proportion of people who have hypnozoites, by G6PD status") #adjust to be exact pN
    plt.tight_layout()
    # fig_relapseG6PD.savefig("./stored/figs/hypnozoites_severe_0.png",bbox_inches='tight')
    fig_relapseG6PD.savefig("./stored/figs/hypnozoites_G6PD_Taf_1461.png",bbox_inches='tight')

plot_relapses_G6PD_2 = False
G6PD_statuses = ["Ineligible", "0-30%", "30-70%","70-100%"]
if plot_relapses_G6PD_2:
    fig_relapseG6PD, ax_relapseG6PD = plt.subplots(1)
    alpha=0.3

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:

        #baseline plotting
        filepath = folder+get_filename("Baseline",[],duration)
        data_G6PD_relapse = load_data(filepath, "G6PD_hypnozoites")
        data_num_G6PD = load_data(filepath, "num_G6PD")
        for G6PD_i in range(len(G6PD_statuses)):
            relapse_G6PD = [data_G6PD_relapse[i][G6PD_i]/data_num_G6PD[G6PD_i] for i in range(len(data_G6PD_relapse))]
            plot_times_G6PD, plot_relapse_G6PD = get_plot_data(relapse_G6PD, plot_start, plot_end, t_step,burnin_days=burnin_days)
            # total_relapses = sum(plot_relapse_G6PD)
            end_relapses = plot_relapse_G6PD[-1]
        # ax_relapseG6PD.plot(times_G6PD,relapse_G6PD,colours[-1],zorder=4)
            ax_relapseG6PD.plot(round(duration/days_in_year-burnin_years),end_relapses, color=colours[G6PD_i], marker="o", alpha=alpha, label = legend_treat[0])

        #Tafenoquine plotting
        treatment_i = 1
        for timechange in timechanges:
            treatment = treatments[treatment_i]
            filepath = folder + get_filename(treatment,timechange,duration)
            data_G6PD_relapse = load_data(filepath, "G6PD_hypnozoites")
            data_num_G6PD = load_data(filepath, "num_G6PD")
            
            for G6PD_i in range(len(G6PD_statuses)):
                relapse_G6PD = [data_G6PD_relapse[i][G6PD_i]/data_num_G6PD[G6PD_i] for i in range(len(data_G6PD_relapse))]
                plot_times_G6PD, plot_relapse_G6PD = get_plot_data(relapse_G6PD, plot_start, plot_end, t_step, burnin_days=burnin_days)
                # total_relapses = sum(plot_relapse_G6PD)
                end_relapses = plot_relapse_G6PD[-1]
                ax_relapseG6PD.plot(round(timechange[0]/days_in_year-burnin_years),end_relapses, color=colours[G6PD_i], marker="o", alpha=alpha, label = legend_treat[0])

    # ax_relapseG6PD.set_ylim(bottom=0)
    ax_relapseG6PD.legend(G6PD_statuses)
    # ax_relapseG6PD.legend()
    ax_relapseG6PD.set_ylabel("Proportion with hypnozoites")
    ax_relapseG6PD.set_xlabel("Time of policy change to predominantly tafenoquine (years)")
    ax_relapseG6PD.set_title("(Fig 3) Proportion of people who have hypnozoites at the end of 10 years, by G6PD status")
    
    timechanges_years_minus_burnin = [round(timechange[0]/days_in_year-burnin_years) for timechange in timechanges]
    xaxis_ticks = timechanges_years_minus_burnin
    xaxis_ticks.append(round(duration/days_in_year-burnin_years))
    
    ax_relapseG6PD.set_xlim([0-2,round(duration/days_in_year-burnin_years)+2])
    ax_relapseG6PD.set_xticks(xaxis_ticks)
    
    plt.tight_layout()
    # fig_relapseG6PD.savefig("./stored/figs/hypnozoites_severe_0.png",bbox_inches='tight')
    fig_relapseG6PD.savefig("./stored/figs/hypnozoites_G6PD_Taf.png",bbox_inches='tight')

plot_cumulative = False
#WIP nearly done
if plot_cumulative:
    fig_cumul_I, ax_cumul_I = plt.subplots(1)
    fig_cumul_relapse, ax_cumul_relapse = plt.subplots(1)

    alpha = 0.3
#     days_in_year = 365.25
#     fig_all_cumul_relapse, ax_all_cumul_relapse = plt.subplots(1)
#     baseline_final_cum_relapse = []
#     handles = [0]*(len(treatments)+1)

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:

        #baseline first
        filepath = folder+get_filename("Baseline",[],duration)
        data_pv_outcomes = load_data(filepath, "pv_outcomes")

        data_I = [data_pv_outcomes[1].count(i) for i in range(duration)]
        data_cumul_I = np.cumsum(data_I)
        plot_times, plot_cumul_I = get_plot_data(data_cumul_I, plot_start, plot_end, t_step, burnin_days = burnin_days)
        ax_cumul_I.plot(plot_times, plot_cumul_I, color = colours_map_treatment["Baseline"], alpha=alpha)#label = legend_treat[0]

        #fill cumulative relapse array
        data_relapse = [data_pv_outcomes[2].count(i) for i in range(duration)]
        data_cumul_relapse = np.cumsum(data_relapse)
        plot_times, plot_cumul_relapse = get_plot_data(data_cumul_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)
        
        ax_cumul_relapse.plot(plot_times, plot_cumul_relapse, color = colours_map_treatment["Baseline"], alpha=alpha)

        print("processing...") 

        #One treatment: Tafenoquine changing at middle timestep
        treatment_i = 1 #Taf
        timechange_i = 1 #Midway
        filepath = folder+get_filename(treatments[treatment_i],timechanges[timechange_i],duration)
        data_pv_outcomes = load_data(filepath, "pv_outcomes")

        #plot I
        data_I = [data_pv_outcomes[1].count(i) for i in range(duration)]
        data_cumul_I = np.cumsum(data_I)
        plot_times, plot_cumul_I = get_plot_data(data_cumul_I, plot_start, plot_end, t_step, burnin_days = burnin_days)
        ax_cumul_I.plot(plot_times, plot_cumul_I, color = colours_map_treatment[treatments[treatment_i]], alpha=alpha)

        #plot relapses
        data_relapse = [data_pv_outcomes[2].count(i) for i in range(duration)]
        data_cumul_relapse = np.cumsum(data_relapse)
        plot_times, plot_cumul_relapse = get_plot_data(data_cumul_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)
        ax_cumul_relapse.plot(plot_times, plot_cumul_relapse, color = colours_map_treatment[treatments[treatment_i]], alpha=alpha)
        print("simulation plotted") 
        break


    # fig_all_cumul_I.legend(handles=handles, loc="center right")
    # ax_all_cumul_I.set_xlabel("Time policy is changed (years)")
    # timechanges_years_minus_burnin = [round(timechange[0]/days_in_year-burnin_years) for timechange in timechanges]
    # xaxis_ticks = timechanges_years_minus_burnin
    # xaxis_ticks.append(round(duration/days_in_year-burnin_years))
    
    # ax_all_cumul_I.set_xlim([0-2,round(duration/days_in_year-burnin_years)+2])
    # ax_all_cumul_I.set_xticks(xaxis_ticks)
    # fig_all_cumul_I.supylabel("Cumulative clinical infections (I)")
    # fig_all_cumul_I.tight_layout()
    
    #Vertical lines
    trans = ax_cumul_I.get_xaxis_transform()
    timechange = round(timechanges[timechange_i][0]-burnin_days)
    ax_cumul_I.axvline(x=timechange,c=colours_map_treatment[treatments[treatment_i]],ls='--',lw=1,label='_nolegend_') #Plot vertical line
    ax_cumul_I.text(timechange, 0.15, " t="+str(timechange), transform=trans,c=colours_map_treatment[treatments[treatment_i]])

    trans = ax_cumul_relapse.get_xaxis_transform()
    ax_cumul_relapse.axvline(x=timechange,c=colours_map_treatment[treatments[treatment_i]],ls='--',lw=1,label='_nolegend_') #Plot vertical line
    ax_cumul_relapse.text(timechange, 0.15, " t="+str(timechange), transform=trans,c=colours_map_treatment[treatments[treatment_i]])


    ax_cumul_I.set_ylim(bottom=0)
    ax_cumul_I.set_title("(Fig 1a) Baseline ("+colours_map_treatment["Baseline"]+") vs "+treatments_titles[treatment_i]+" ("+colours_map_treatment[treatments[treatment_i]]+")")
    ax_cumul_I.set_title("All scenarios")

    fig_cumul_I.legend(legend_treat)
    ax_cumul_I.set_xlabel("Time (days)")
    ax_cumul_I.set_xlabel("Time policy is changed (days)")
    # ax_cumul_I.set_xlim([-2,10])
    # ax_cumul_I.set_xticks([0,4,8])
    fig_cumul_I.supylabel("Cumulative clinical infections (I)")
    fig_cumul_I.tight_layout()
    fig_cumul_I.savefig("./stored/figs/cumul_I.png",bbox_inches='tight')

    ax_cumul_relapse.set_ylim(bottom=0)
    ax_cumul_relapse.set_title("(Fig 2a) Baseline ("+colours_map_treatment["Baseline"]+") vs "+treatments_titles[-1]+" ("+colours_map_treatment[treatments[treatment_i]]+")")
    ax_cumul_relapse.set_title("All scenarios")

    fig_cumul_relapse.legend(legend_treat)
    ax_cumul_relapse.set_xlabel("Time (days)")
    ax_cumul_relapse.set_xlabel("Time policy is changed (days)")
    # ax_cumul_relapse.set_xlim([-2,10])
    # ax_cumul_relapse.set_xticks([0,4,8])
    fig_cumul_relapse.supylabel("Cumulative relapses")
    fig_cumul_relapse.tight_layout()
    fig_cumul_relapse.savefig("./stored/figs/cumul_relapse.png",bbox_inches='tight')

        # relapses
        # total_relapse = len(data_pv_outcomes[2])
        # baseline_final_cum_relapse.append(total_relapse)
        
    # median_baseline_relapse = np.median(baseline_final_cum_relapse)
    # scale_factor = 1/median_baseline_relapse

    # for total_relapse in baseline_final_cum_relapse:
        # handles[0], = ax_all_cumul_relapse.plot(round(duration/days_in_year-burnin_years),total_relapse*scale_factor, color=colours_map_treatment["Baseline"], marker="o", alpha=alpha, label = legend_treat[0])


    # #plot a

    # filepaths = []
    # filepaths.append(get_filepath("Baseline",[],duration))
    # filepaths.append(get_filepath(treatments[-1],timechanges[0],duration))
    # fig_relapse.tight_layout(rect=[0, 0.03, 1, 0.95])
    # for filepath_i in range(len(filepaths)):

        # #pv data
        # data_pv_outcomes = load_data(filepaths[filepath_i], "pv_outcomes")
        # # treat_baseline_T = [sum([compartment_data_baseline_3[comp].count(i) for comp in [3]]) for i in range(duration)]
        # data_I = [data_pv_outcomes[1].count(i) for i in range(duration)]
        # data_cumul_I = np.cumsum(data_I)
        
        # #plot cumulative I
        # plot_times, plot_cumul_I = get_plot_data(data_cumul_I, plot_start, plot_end, t_step)
        # ax_cumul_I[0].plot(plot_times, plot_cumul_I, colours[filepath_i])



        # #relapse data
        # data_relapse = [data_pv_outcomes[2].count(i) for i in range(duration)]
        # #fill cumulative relapse array
        # data_cumul_relapse = np.cumsum(data_relapse)
        
        # #plot cumulative relapses
        # plot_times, plot_cumul_relapse = get_plot_data(data_cumul_relapse, plot_start, plot_end, t_step)
        # ax_cumul_relapse[0].plot(plot_times, plot_cumul_relapse, colours[filepath_i])

        
    # ax_cumul_I[0].set_ylim(bottom=0)
    # ax_cumul_I[0].set_title("(Fig 1a) Baseline ("+colours[0]+") vs "+treatments_titles[-1]+" ("+colours[1]+")")
    # ax_cumul_I[1].set_title("All scenarios")

    # fig_cumul_I.legend(legend_timechange)
    # ax_cumul_I[0].set_xlabel("Time (days)")
    # ax_cumul_I[1].set_xlabel("Time policy is changed (days)")
    # ax_cumul_I[1].set_xlim([-2,10])
    # ax_cumul_I[1].set_xticks([0,4,8])
    # fig_cumul_I.supylabel("Cumulative clinical infections (I)")
    # fig_cumul_I.tight_layout()
    # fig_cumul_I.savefig("./stored/figs/cumul_I.png",bbox_inches='tight')

    # ax_cumul_relapse[0].set_ylim(bottom=0)
    # ax_cumul_relapse[0].set_title("(Fig 2a) Baseline ("+colours[0]+") vs "+treatments_titles[-1]+" ("+colours[1]+")")
    # ax_cumul_relapse[1].set_title("All scenarios")

    # fig_cumul_relapse.legend(legend_timechange)
    # ax_cumul_relapse[0].set_xlabel("Time (days)")
    # ax_cumul_relapse[1].set_xlabel("Time policy is changed (days)")
    # ax_cumul_relapse[1].set_xlim([-2,10])
    # ax_cumul_relapse[1].set_xticks([0,4,8])
    # fig_cumul_relapse.supylabel("Cumulative relapses")
    # fig_cumul_relapse.tight_layout()
    # fig_cumul_relapse.savefig("./stored/figs/cumul_relapse.png",bbox_inches='tight')


##### Ending value of cumulative clinical infections, all policies and simulations
all_policies_cum_I = False
if all_policies_cum_I:

    alpha = 0.3
    fig_all_cumul_I, ax_all_cumul_I = plt.subplots(1)
    baseline_final_cum_I = []
    handles = [0]*(len(treatments)+1)

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:
        #baseline first
        filepath = folder+get_filename("Baseline",[],duration)
        data_pv_outcomes = load_data(filepath, "pv_outcomes")
        data_I = [data_pv_outcomes[1].count(i) for i in range(duration)]
        plot_times, plot_I = get_plot_data(data_I, plot_start, plot_end, t_step, burnin_days = burnin_days)
        total_I = sum(plot_I)
        baseline_final_cum_I.append(total_I)
        
    median_baseline_I = np.median(baseline_final_cum_I)
    scale_factor = 1/median_baseline_I
    for total_I in baseline_final_cum_I:
        handles[0], = ax_all_cumul_I.plot(round(duration/days_in_year-burnin_years),total_I*scale_factor, color=colours_map_treatment["Baseline"], marker="o", alpha=alpha, label = legend_treat[0])

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:# else:
        for timechange in timechanges:
            for treatment_i in range(len(treatments)):
                treatment = treatments[treatment_i]
                filepath = folder + get_filename(treatment,timechange,duration)
                data_pv_outcomes = load_data(filepath, "pv_outcomes")
                data_I = [data_pv_outcomes[1].count(i) for i in range(duration)]
                plot_times, plot_I = get_plot_data(data_I, plot_start, plot_end, t_step, burnin_days = burnin_days)
                total_I = sum(plot_I)
                handles[treatment_i+1], = ax_all_cumul_I.plot(round(timechange[0]/days_in_year-burnin_years),total_I*scale_factor, color = colours_map_treatment[treatment], marker="o", alpha=alpha, label = legend_treat[treatment_i+1])

    ax_all_cumul_I.set_title("(Fig 1b) Cumulative infections as a proportion of baseline median, with policy changed at different time points")
    # fig_all_cumul_I.legend(legend_treat, loc="center right")
    fig_all_cumul_I.legend(handles=handles, loc="center right")
    ax_all_cumul_I.set_xlabel("Time policy is changed (years)")
    timechanges_years_minus_burnin = [round(timechange[0]/days_in_year-burnin_years) for timechange in timechanges]
    xaxis_ticks = timechanges_years_minus_burnin
    xaxis_ticks.append(round(duration/days_in_year-burnin_years))
    
    ax_all_cumul_I.set_xlim([0-2,round(duration/days_in_year-burnin_years)+2])
    ax_all_cumul_I.set_xticks(xaxis_ticks)
    fig_all_cumul_I.supylabel("Cumulative clinical infections (I)")
    fig_all_cumul_I.tight_layout()

    fig_all_cumul_I.savefig("./stored/figs/cumul_I_all_policies.png",bbox_inches='tight')

# ##### Time series cumulative relapses, Taf policy and all simulations
# cum_relapse = True
# if cum_relapse:

#     alpha = 0.3
#     days_in_year = 365.25
#     fig_all_cumul_relapse, ax_all_cumul_relapse = plt.subplots(1)
#     baseline_final_cum_relapse = []
#     handles = [0]*(len(treatments)+1)

#     for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:
#         #baseline first
#         filepath = folder+get_filename("Baseline",[],duration)
#         data_pv_outcomes = load_data(filepath, "pv_outcomes")
#         total_relapse = len(data_pv_outcomes[2])
#         baseline_final_cum_relapse.append(total_relapse)
        
#     median_baseline_relapse = np.median(baseline_final_cum_relapse)
#     scale_factor = 1/median_baseline_relapse

#     for total_relapse in baseline_final_cum_relapse:
#         handles[0], = ax_all_cumul_relapse.plot(round(duration/days_in_year-burnin_years),total_relapse*scale_factor, color=colours_map_treatment["Baseline"], marker="o", alpha=alpha, label = legend_treat[0])

#     for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:# else:
#         for timechange in timechanges:
#             for treatment_i in range(len(treatments)):
#                 treatment = treatments[treatment_i]
#                 filepath = folder + get_filename(treatment,timechange,duration)
#                 data_pv_outcomes = load_data(filepath, "pv_outcomes")
#                 total_relapse = len(data_pv_outcomes[2])
#                 handles[treatment_i+1], = ax_all_cumul_relapse.plot(round(timechange[0]/days_in_year-burnin_years),total_relapse*scale_factor, color = colours_map_treatment[treatment], marker="o", alpha=alpha, label = legend_treat[treatment_i+1])

#     ax_all_cumul_relapse.set_title("(Fig 2b) Cumulative relapses as a proportion of baseline median, with policy changed at different time points")
#     # fig_all_cumul_relapse.legend(legend_treat, loc="center right")
#     fig_all_cumul_relapse.legend(handles=handles, loc="center right")
#     ax_all_cumul_relapse.set_xlabel("Time policy is changed (years)")
#     timechanges_years_minus_burnin = [round(timechange[0]/days_in_year-burnin_years) for timechange in timechanges]
#     xaxis_ticks = timechanges_years_minus_burnin
#     xaxis_ticks.append(round(duration/days_in_year-burnin_years))
    
#     ax_all_cumul_relapse.set_xlim([0-2,round(duration/days_in_year-burnin_years)+2])
#     ax_all_cumul_relapse.set_xticks(xaxis_ticks)
#     fig_all_cumul_relapse.supylabel("Cumulative relapses")
#     fig_all_cumul_relapse.tight_layout()

#     fig_all_cumul_relapse.savefig("./stored/figs/cumul_relapse_all_policies.png",bbox_inches='tight')


##### Ending value of cumulative relapses, all policies and simulations
all_policies_cum_relapse = False
if all_policies_cum_relapse:

    alpha = 0.3
    fig_all_cumul_relapse, ax_all_cumul_relapse = plt.subplots(1)
    baseline_final_cum_relapse = []
    handles = [0]*(len(treatments)+1)

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:
        #baseline first
        filepath = folder+get_filename("Baseline",[],duration)
        data_relapse = load_data(filepath, "pv_actual_relapses")
        plot_times, plot_relapse = get_plot_data(data_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)
        total_relapse = sum(plot_relapse)
        baseline_final_cum_relapse.append(total_relapse)
        
    median_baseline_relapse = np.median(baseline_final_cum_relapse)
    scale_factor = 1/median_baseline_relapse

    for total_relapse in baseline_final_cum_relapse:
        handles[0], = ax_all_cumul_relapse.plot(round(duration/days_in_year-burnin_years),total_relapse*scale_factor, color=colours_map_treatment["Baseline"], marker="o", alpha=alpha, label = legend_treat[0])

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration))]:# else:
        for timechange in timechanges:
            for treatment_i in range(len(treatments)):
                treatment = treatments[treatment_i]
                filepath = folder + get_filename(treatment,timechange,duration)
                # data_pv_outcomes = load_data(filepath, "pv_outcomes")
                data_relapse = load_data(filepath, "pv_actual_relapses")
                plot_times, plot_relapse = get_plot_data(data_relapse, plot_start, plot_end, t_step, burnin_days = burnin_days)
                # total_relapse = len(data_pv_outcomes[2])
                total_relapse = sum(plot_relapse)
                handles[treatment_i+1], = ax_all_cumul_relapse.plot(round(timechange[0]/days_in_year-burnin_years),total_relapse*scale_factor, color = colours_map_treatment[treatment], marker="o", alpha=alpha, label = legend_treat[treatment_i+1])

    ax_all_cumul_relapse.set_title("(Fig 2b) Cumulative relapses as a proportion of baseline median, with policy changed at different time points")
    # fig_all_cumul_relapse.legend(legend_treat, loc="center right")
    fig_all_cumul_relapse.legend(handles=handles, loc="center right")
    ax_all_cumul_relapse.set_xlabel("Time policy is changed (years)")
    timechanges_years_minus_burnin = [round(timechange[0]/days_in_year-burnin_years) for timechange in timechanges]
    xaxis_ticks = timechanges_years_minus_burnin
    xaxis_ticks.append(round(duration/days_in_year-burnin_years))
    
    ax_all_cumul_relapse.set_xlim([0-2,round(duration/days_in_year-burnin_years)+2])
    ax_all_cumul_relapse.set_xticks(xaxis_ticks)
    fig_all_cumul_relapse.supylabel("Cumulative relapses")
    fig_all_cumul_relapse.tight_layout()

    fig_all_cumul_relapse.savefig("./stored/figs/cumul_relapse_all_policies.png",bbox_inches='tight')

#plot histogram of hypnozoite distribution at specified timepoints
hypnozoite_dist = True
if hypnozoite_dist:

    time_plots = params.hyp_snapshot_days
    # assert all(time in params.hyp_snapshot_years for time in time_plots)
    n_subplots = len(time_plots)#2#len(params.hyp_snapshot_days)
    fig_hyp_dist, ax_hyp_dist = plt.subplots(n_subplots)

    for folder in ["./stored/results_variables/" + folder + "/" for folder in os.listdir('./stored/results_variables/') if folder.startswith("duration_"+str(duration)+"_24")]:
        
        print("****")
        #baseline first
        # filepath = folder+get_filename(treatments[0],timechanges[0],duration)
        filepath = folder+get_filename(treatments[0],timechanges[0],duration)
        data_hyp_dist = load_data(filepath, "hypnozoite_dist")
        
        #Test data
        # data_hyp_dist = [[1,[30792, 16067, 10799, 7960, 6015, 4706, 3713, 2957, 2380, 2042, 1714, 1449, 1272, 1014, 890, 753, 637, 520, 452, 404, 391, 304, 245, 249, 199, 189, 139, 144, 114, 105, 103, 66, 72, 51, 51, 42, 48, 37, 25, 30, 23, 19, 15, 21, 16, 13, 10, 5, 3, 5, 11, 6, 2, 6, 7, 4, 5, 3, 3, 2, 2, 1, 2, 0, 4, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [2,[30792, 16067, 10799, 7960, 6015, 4706, 3713, 2957, 2380, 2042, 1714, 1449, 1272, 1014, 890, 753, 637, 520, 452, 404, 391, 304, 245, 249, 199, 189, 139, 144, 114, 105, 103, 66, 72, 51, 51, 42, 48, 37, 25, 30, 23, 19, 15, 21, 16, 13, 10, 5, 3, 5, 11, 6, 2, 6, 7, 4, 5, 3, 3, 2, 2, 1, 2, 0, 4, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], [3,[30792, 16067, 10799, 7960, 6015, 4706, 3713, 2957, 2380, 2042, 1714, 1449, 1272, 1014, 890, 753, 637, 520, 452, 404, 391, 304, 245, 249, 199, 189, 139, 144, 114, 105, 103, 66, 72, 51, 51, 42, 48, 37, 25, 30, 23, 19, 15, 21, 16, 13, 10, 5, 3, 5, 11, 6, 2, 6, 7, 4, 5, 3, 3, 2, 2, 1, 2, 0, 4, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], [4,[30792, 16067, 10799, 7960, 6015, 4706, 3713, 2957, 2380, 2042, 1714, 1449, 1272, 1014, 890, 753, 637, 520, 452, 404, 391, 304, 245, 249, 199, 189, 139, 144, 114, 105, 103, 66, 72, 51, 51, 42, 48, 37, 25, 30, 23, 19, 15, 21, 16, 13, 10, 5, 3, 5, 11, 6, 2, 6, 7, 4, 5, 3, 3, 2, 2, 1, 2, 0, 4, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]]

        # for time in params.hyp_snapshot_days:
        count=0
        # for subplot in [0,2]:
        for subplot in range(n_subplots):
            vals = data_hyp_dist[subplot][1]
            bins = range(len(vals))
            ax_hyp_dist[count].bar(bins,vals)
            ax_hyp_dist[count].set_title("Snapshot at time t = "+str(time_plots[count])+" years")#round(params.hyp_snapshot_days[subplot]/days_in_year)-burnin_years)+" years")

            count += 1
    # fig_hyp_dist.suptitle("Hypnozoite distribution before and after a policy change [low-dose PQ]")
    fig_hyp_dist.supxlabel("Number of hypnozoites")
    fig_hyp_dist.supylabel("Number of people")
    fig_hyp_dist.savefig("./stored/figs/hyp_dist.png",bbox_inches='tight')




plt.show(block=False)
print("Plots finished. Press ENTER to continue.")
input()
plt.close()

