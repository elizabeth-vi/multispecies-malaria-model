import matplotlib.pyplot as plt
import json
import numpy as np
from model_params import model_params
import itertools
from index_names import Compartments, Mozzie_labels


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
    filepath_loc = "./stored/results_variables/duration_" + str(duration) + "/"
    filepath = filepath_loc + treatment + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type
    return filepath

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
    
def get_plot_data(data, t_start, t_end, t_step):

    start_index = int(t_start/t_step)
    end_index = int(t_end/t_step)

    plot_times = np.arange(start=start_index*t_step, stop=end_index*t_step, step=t_step)
    plot_data = data[start_index:end_index]

    return plot_times, plot_data

#*********************************SCRIPT START*********************************


#Filepath name info
#[filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]


params = model_params()
treatments = ["PLD","PHD","Taf"] #Must match file names. Options: ["PLD","PHD","Taf"]]
treatments_titles = ["Change to Policy 1", "Change to Policy 2", "Change to Policy 3"]
colours = ['r','b','g', 'k'] #for plotting
timechanges = [[0], [1461], [2922]] #which time changes you want to plot, in days. e.g.: [0, 730, 1461]
duration = int(params.time_day_end) #Simulation duration

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start, duration, params.time_day_step] #Time of start/end you want to plot, and timestep used
plt.rcParams["figure.figsize"] = (11,8.5)
# plt.tight_layout()

#Set up figures
# fig_infect = plt.figure()
# fig_IAL = plt.figure()
# fig_TG = plt.figure()
legend_treat = ["Baseline"]
legend_treat.extend(treatments)
legend_timechange = ["Baseline"]
legend_timechange.extend(["Change policy at t = "+str(x[0])+" days" for x in timechanges])

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
plot_infections = False
if plot_infections:
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

#plot relapses
plot_relapses = True
if plot_relapses:
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
plot_relapses_G6PD = True
G6PD_statuses = ["Ineligible", "0-30%", "30-70%","70-100%"]
if plot_relapses_G6PD:
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

# data = load_data(get_filepath(treatment="Taf", timechange="120",duration=plot_end), "human_pop_pv_history")
# plot_all_human_pv_compartments(treatment="Taf", timechange="120", data=data, t_start=plot_start, t_end=plot_end, t_step=t_step)

# fig_infect.legend(legend,loc=7)
# fig_IAL.legend(legend,loc=7)
# fig_TG.legend(legend,loc=7)





plt.show(block=False)
input()
plt.close()
