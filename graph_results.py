from pyexpat import model
import matplotlib.pyplot as plt
import json
import numpy as np
from model_parameters import model_parameters
from operator import add

plt.close("all")
#plt.ion()

filepath = "C:\\Users\eliza\Desktop\RA\multispecies-malaria-model\stored\entangled_treatmentExample_Province_all_pN1.0_0.843_scenario0.json"

dictionary = json.load(open(filepath, 'r'))

params = model_parameters()
#times = np.arange(params.time_day_start,params.time_day_end,params.time_day_step)
#print(times)

variables = [key for key, value in dictionary.items()]
print(variables)
#print(len(variables[0]))

plot_selected = True
if plot_selected:
    yvals = []
    pf_hist = variables[0]
    pv_hist = variables[1]
    mixed_hist = variables[2]
    times = np.linspace((params.time_day_start),int(params.time_day_end),len(dictionary.get(pf_hist)[0]))
    
    vals_pf = dictionary.get(pf_hist)[0]
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
    plt.figure()
    plt.grid(True)
    plt.plot(times,pf_only_data, color='maroon')
    plt.xlabel('time')
    plt.ylabel('Infections')
    plt.title('Human p.f. infections (I and A)')
    
    #Plot p.v.
    plt.figure()
    plt.grid(True)
    plt.plot(times,pv_only_data, color='maroon')
    plt.xlabel('time')
    plt.ylabel('Infections')
    plt.title('Human p.v. infections (I and A)')
    #"""



plot_all = True
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

#xAxis = [key for key, value in dictionary.items()]
#yAxis = [value for key, value in dictionary.items()]

"""
plt.grid(True)

## LINE GRAPH ##
plt.plot(xAxis,yAxis, color='maroon', marker='o')
plt.xlabel('variable')
plt.ylabel('value')
"""

"""
## BAR GRAPH ##
fig = plt.figure()
plt.bar(xAxis,yAxis, color='maroon')
plt.xlabel('variable')
plt.ylabel('value')
"""

plt.show(block=False)
input()
plt.close()
