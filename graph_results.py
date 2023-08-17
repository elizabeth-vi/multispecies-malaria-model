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

    variables = [*data.keys()]
    
    #Plot by compartments
    pv_hist = variables[1]
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = vals_pv[t_start : t_start+n_vals]

    for compartment in Compartments:
        pv_compartment_data = [val[compartment.value] for val in vals_pv]

        #Plot human p.v. infections
        human_axs[compartment.value].plot(times,pv_compartment_data)
        human_axs[compartment.value].set_title(compartment.name,x=1.03,y=0)
        human_axs[compartment.value].grid()
        human_axs[compartment.value].set_ylim(bottom=0)

def plot_human_select_pv_compartments_together(data, t_start, t_end, t_step, compartments_str, fig=None, colour='b'):
    #Set up fig
    fig = plt.figure(fig)

    variables = [*data.keys()]
    compartments = [Compartments[comp] for comp in compartments_str]
    
    #Plot by compartments
    pv_hist = variables[1]
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = vals_pv[t_start : t_start+n_vals]

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

def get_compartment_data(treatment, timechange, duration):

    #Filename / filepath format 
    [filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]
    #duration = str(int(params.time_day_end)) #Simulation duration. Used in filename
    filepath_loc = "./stored/results_variables/duration_" + str(duration) + "/"
    filepath = filepath_loc + treatment + filepath_timechange + str(timechange) + filepath_duration + str(duration) + filepath_type

    #Load and return data
    try:
        data = json.load(open(filepath, 'r'))
        vals_pv = data.get("human_pop_pv_history")[0]
        return vals_pv
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
treatments = ["Primaquine_Lowdose","Primaquine_Highdose","Tafenoquine"] #Must match file names. Options: ["Primaquine_Lowdose","Primaquine_Highdose","Tafenoquine"]]
colours = ['r','b','g', 'k'] #for plotting
timechanges = [0, 36, 73] #which time changes you want to plot, in days. e.g.: [0, 730, 1461]
duration = int(params.time_day_end) #Simulation duration

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start, duration, params.time_day_step] #Time of start/end you want to plot, and timestep used

#Set up figures
# fig_infect = plt.figure()
# fig_IAL = plt.figure()
# fig_TG = plt.figure()
legend_treat = ["baseline"]
legend_treat.extend(treatments)
legend_timechange = ["baseline"]
legend_timechange.extend(["Change at t = "+str(x)+" days" for x in timechanges])

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

#Set up fig 2
fig_infect, ax_infect = plt.subplots(len(treatments))
fig_infect.tight_layout(rect=[0, 0.03, 1, 0.95])
fig_infect.suptitle("Infections when treatment policy is changed at different times")

#Set up fig 1

for treatment_i in range(len(treatments)):
    
    #baseline plotting
    compartment_data_baseline = get_compartment_data("baseline",duration,duration)
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
        data = get_compartment_data(treatment,timechange,duration)
        infections = [val[1] for val in data]
        times, data_I = get_plot_data(infections, plot_start, plot_end, t_step)
        
        #Plot
        #For each treatment, plot I for each timechange option
        ax_infect[treatment_i].plot(times,data_I,colour)
        ax_infect[treatment_i].axvline(x=timechange,c=colour,ls='--',lw=1,label='_nolegend_') #Plot vertical line
        ax_infect[treatment_i].text(timechange, 0.5, " t="+str(timechange), transform=trans,c=colour)
    
    ax_infect[treatment_i].set_ylim(bottom=0)
    ax_infect[treatment_i].legend(legend_timechange)
    ax_infect[treatment_i].set_title("Treatment policy changes to "+treatments[treatment_i]+" with higher treatment rates") #adjust to be exact pN

    # plot_all_human_pv_compartments(treatment, "0", get_compartment_data(treatment,0,duration), plot_start, plot_end, t_step)

        


if False:
    for treatment_i, timechange_i in itertools.product(range(len(treatments)), range(len(timechanges))):

        # plot_number = timechange_i + treatment_i*len(timechanges)
        treatment = treatments[treatment_i]
        timechange = timechanges[timechange_i]

        filepath = filepath_loc + treatment + filepath_timechange + timechange + filepath_duration + duration + filepath_type
        data = json.load(open(filepath, 'r'))
        #variables = [key for key, value in data.items()]

        # plot_all_human_pv_compartments(data, plot_start, plot_end, t_step)
        fig_infect = plot_human_pv_infections(data, plot_start, plot_end, t_step, fig_infect)
        # fig_infect.legend(treatments)
        legend.append(treatment+" introduced at time "+timechange)


        # fig_IAL = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["I","A","L"], fig=fig_IAL)
        # fig_TG = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["T","G"], fig=fig_TG)

        fig_IAL = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["A"], fig=fig_IAL, colour = colours[treatment_i])
        fig_TG = plot_human_select_pv_compartments_together(data, plot_start, plot_end, t_step, compartments_str=["G"], fig=fig_TG, colour = colours[treatment_i])
        #plot_all_human_pv_compartments( data, plot_start, plot_end, t_step)
                    
        ###WIP
        #plot_all_tracked_variables(data, plot_start, plot_end, t_step)

# fig_infect.legend(legend,loc=7)
# fig_IAL.legend(legend,loc=7)
# fig_TG.legend(legend,loc=7)


plt.show(block=False)
input()
plt.close()
