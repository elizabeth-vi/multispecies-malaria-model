import matplotlib.pyplot as plt
import json
import numpy as np
from model_params import model_params
import itertools
import index_names


def plot_all_human_pv_compartments(treatment, timechange, data, t_start, t_end, t_step):
    #Set up fig
    human_fig, human_axs = plt.subplots(len(index_names.Compartments))
    human_fig.supxlabel('time (days)')
    human_fig.suptitle('Human p.v. Compartments for '+treatment+" introduced at time "+timechange)

    variables = [*data.keys()]
    
    #Plot by compartments
    pv_hist = variables[1]
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = vals_pv[t_start : t_start+n_vals]

    for compartment in index_names.Compartments:
        pv_compartment_data = [val[compartment.value] for val in vals_pv]

        #Plot human p.v. infections
        human_axs[compartment.value].plot(times,pv_compartment_data)
        human_axs[compartment.value].set_title(compartment.name,x=1.03,y=0)
        human_axs[compartment.value].grid()
        human_axs[compartment.value].set_ylim(bottom=0)

def plot_human_select_pv_compartments_together(treatment, timechange, data, t_start, t_end, t_step, compartments_str):
    #Set up fig
    plt.figure()

    variables = [*data.keys()]
    compartments = [index_names.Compartments[comp] for comp in compartments_str]
    
    #Plot by compartments
    pv_hist = variables[1]
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = vals_pv[t_start : t_start+n_vals]

    for index in range(len(compartments)):
        compartment = compartments[index]
        #compartment = index_names.Compartments(vars(compartments[index]))
        pv_compartment_data = [val[compartment.value] for val in vals_pv]

        #Plot human p.v. infections
        plt.plot(times,pv_compartment_data)
    
    legend = [compartment.name for compartment in compartments]
    plt.xlabel('time (days)')
    plt.title("Human p.v. Compartments "+str(legend)+" for "+treatment+" introduced at time "+timechange)
    plt.ylim(bottom=0)
    plt.legend(legend)
    plt.grid()

def plot_human_pv_infections(treatment, timechange, data, t_start, t_end, t_step):
    #Set up fig
    plt.figure()

    variables = [*data.keys()]
    
    #Plot by compartments
    pv_hist = variables[1] #human_pop_pv_history
    vals_pv = data.get(pv_hist)[0]
    
    times = np.arange(start=t_start, stop=t_end, step=t_step)
    n_vals = len(times)
    
    vals_pv = vals_pv[t_start : t_start+n_vals]

    pv_compartment_data = [val[index_names.Compartments['I']]+val[index_names.Compartments['A']] for val in vals_pv]

    #Plot human p.v. infections
    plt.plot(times,pv_compartment_data)
    plt.xlabel('time (days)')
    plt.ylabel('Infections')
    plt.title('Human p.v. infections (I + A)'+" for "+treatment+" introduced at time "+timechange)
    plt.ylim(bottom=0) #set lower bound for plotting
    plt.grid()

#WORK IN PROGRESS
def plot_all_tracked_variables(treatment, timechange, data, t_start, t_end, t_step):

    #variables = [*data.keys()] #For all variables
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
                    
                #print(index_names.Compartments._member_names_)  
                plt.legend(index_names.Compartments._member_names_)
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

#*********************************SCRIPT START*********************************

#File locations
filepath_loc = "./stored/results_variables/"
#Filepath name info
[filepath_timechange,filepath_duration,filepath_type] = ["_timechange","_duration",".json"]


params = model_params()
treatments = ["Tafenoquine"] #Must match file names. Options: ["Primaquine_Lowdose","Primaquine_Highdose","Tafenoquine"]]
timechanges = ["730"] #which time changes you want to plot, in days. Options: ["0", "730", "1461"]
duration = str(int(params.time_day_end)) #Simulation duration. Used in filename

plt.close("all")
nplots = len(treatments)*len(timechanges)
[plot_start, plot_end, t_step] = [params.time_day_start,params.time_day_end,params.time_day_step] #Time of start/end you want to plot, and timestep used

for treatment_i, timechange_i in itertools.product(range(len(treatments)), range(len(timechanges))):

    plot_number = timechange_i + treatment_i*len(timechanges)
    treatment = treatments[treatment_i]
    timechange = timechanges[timechange_i]

    filepath = filepath_loc + treatment + filepath_timechange + timechange + filepath_duration + duration+ filepath_type
    data = json.load(open(filepath, 'r'))
    variables = [key for key, value in data.items()]

    plot_all_human_pv_compartments(treatment, timechange, data, plot_start, plot_end, t_step)
    plot_human_pv_infections(treatment, timechange, data, plot_start, plot_end, t_step)
    plot_human_select_pv_compartments_together(treatment, timechange, data, plot_start, plot_end, t_step, compartments_str=["I","A","L"])
    plot_human_select_pv_compartments_together(treatment, timechange, data, plot_start, plot_end, t_step, compartments_str=["T","G"])
    
    ###WIP
    #plot_all_tracked_variables(treatment, timechange, data, plot_start, plot_end, t_step)


plt.show(block=False)
input()
plt.close()
