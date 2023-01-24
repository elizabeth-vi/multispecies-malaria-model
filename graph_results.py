from pyexpat import model
import matplotlib.pyplot as plt
import json
import numpy as np
from model_params import model_params
from operator import add

plt.close("all")
#plt.ion()

#old filepath = "/Users/elizabeth/Desktop/RA/multispecies-malaria-model/stored/entangled_treatmentExample_Province_all_pN1.0_0.843_scenario0.json"
#"C:\Users\eliza\Desktop\RA\multispecies-malaria-model\stored\entangled_treatment_Example_Province_Primaquine_Highdoseyear048.json"

#filepath = "./stored\entangled_treatment_all_Example_Province_Primaquine_Highdose_timechange_0.json"
filepaths = ["./stored/results\entangled_treatment_all_Example_Province_Tafenoquine_timechange_0.json","./stored/results\entangled_treatment_all_Example_Province_Tafenoquine_timechange_183.json","./stored/results\entangled_treatment_all_Example_Province_Tafenoquine_timechange_365.json"]
titles = ["Tafenoquine introduced at t=0","Tafenoquine introduced at t=183","Tafenoquine introduced at t=365"]

params = model_params()
#times = np.arange(params.time_day_start,params.time_day_end,params.time_day_step)
#print(times)


#print(len(variables[0]))

human_fig, human_axs = plt.subplots(3)
human_fig.supxlabel('time')
human_fig.supylabel('Infections')
human_fig.suptitle('Human p.v. infections (I and A)')
plt.subplots_adjust(hspace=0.5)


mosq_fig, mosq_axs = plt.subplots(3)
mosq_fig.supxlabel('time')
mosq_fig.supylabel('Infections')
mosq_fig.suptitle('Mozzie p.v. infections')
plt.subplots_adjust(hspace=0.5)

for plot_number in range(len( filepaths)):
    filepath = filepaths[plot_number]
    title = titles[plot_number]

    dictionary = json.load(open(filepath, 'r'))
    variables = [key for key, value in dictionary.items()]

    

    plot_selected = True
    if plot_selected:

        #Plot human p.v. infections
        yvals = []
        pf_hist = variables[0] #name of variable, 'human_pop_pf_history'
        pv_hist = variables[1]
        mixed_hist = variables[2]
        times = np.linspace((params.time_day_start),int(params.time_day_end),len(dictionary.get(pf_hist)[0]))
        
        vals_pf = dictionary.get(pf_hist)[0] #get values associated with first variable
        vals_pv = dictionary.get(pv_hist)[0]
        vals_mixed = dictionary.get(mixed_hist)[0]
        
        #yvals = yvals + vals[0]
        #yvals = list(map(add,[val[1]+val[2] for val in vals_pf],[val[1]+val[2] for val in vals_mixed]))#doesn't work
        pf_only_data = [val[1]+val[2] for val in vals_pf]
        pv_only_data = [val[1]+val[2] for val in vals_pv]
        #yvals = pv_only_data
        #mixed_data = [val[1]+val[2] for val in vals_mixed]
        #print(mixed_data)
        #yvals = list(map(add,pf_only_data,mixed_data))
        
        #print(yvals)
        
        #"""
        
        #Plot p.f.
        # plt.figure()
        # plt.grid(True)
        # plt.plot(times,pf_only_data, color='maroon')
        # plt.xlabel('time')
        # plt.ylabel('Infections')
        # plt.title('Human p.f. infections (I and A)')
        
        #Plot human p.v. infections
        
        human_axs[plot_number].plot(times,pv_only_data)
        human_axs[plot_number].set_title(title)
        human_axs[plot_number].grid()
        
        # plt.figure()
        # plt.grid(True)
        # plt.plot(times,pv_only_data, color='maroon')
        # plt.xlabel('time')
        # plt.ylabel('Infections')
        # plt.title(title + '\nHuman p.v. infections (I and A)')

        #"""


        #Plot mosquito p.v. infections
        mozzie_pv = variables[4]
        mozzie_pv_data = dictionary.get(variables[4])[0]

        #Plot mosquito p.v. infections
        mosq_axs[plot_number].plot(times, mozzie_pv_data)
        mosq_axs[plot_number].set_title(title)
        mosq_axs[plot_number].grid()
        #"""



    plot_all = False
    if plot_all:
        for variable in variables:
            #print(variable)
            vals = dictionary.get(variable)[0]
            #print(vals)
            times = np.linspace((params.time_day_start),int(params.time_day_end),len(vals))
            #print([time for time in times])
            #print(times)
            #print(np.size(vals))
            #print(len(vals))
            #print(vals[0])
            #print(vals)

            #length = len(vals[0])
            #print(length)

            plt.figure()
            plt.grid(True)

            if isinstance(vals[0],list):
                #print(len(vals[0]))
                #if len(vals[0])>20:
                    #plt.close()
                    #break

                try:
                    for i in range(len(vals[0])):
                        #plot
                        #print(i)
                        yvals = [val[i] for val in vals]
                        #print(yvals)
                        
                        #plt.figure(i)

                        ## LINE GRAPH ##
                        plt.plot(times,yvals, color='maroon')
                        #plt.setp(fig,linewidth=0.5) #set fig= plt.plot...
                        plt.plot()
                        #if i in [1,2]:
                        #    plt.legend('I or A')
                    plt.xlabel('time')
                    plt.ylabel(variable)
                    plt.title(variable)
                    #print(variable)
                except:
                    plt.close()
                    break
            else:
                yvals = vals
                #print(yvals)

                ## LINE GRAPH ##
                plt.plot(times,yvals, color='maroon')
                plt.xlabel('time')
                plt.ylabel(variable)
                plt.title(variable)

plt.show(block=False)
input()
plt.close()
