'''
Created on 20 mar. 2020

@author: Miguel
'''
from disease import DiseaseSimulation
from copy import copy
import numpy as np

def stepOptimizer(h_max, tolerance=0.2, **diseaseKwargs):
    """ 
    This function iterates to find a value of the step for which the SIR simulation
    converges under a certain tolerance. The steps are divided by 2 in each 
    iteration (for a maximum of h_max/2^-6).
    Args:
    :h_max first step value.
    :tolerance (=0.2 by default), tolerance ratio on step (1.0 for 100%).
    :diseaseKwargs, are the parameters for the DiseaseSimulation.
    
    Return:
    :h optimized or not.
    """
    
    max_vals_prev = (0,0)
    max_values = {}
    optimiced = False
    DiseaseSimulation.stopWhenMaxInfectedReached(True)
    for i in range(7):
        h = h_max / (2**(i))
        
        ds_h = DiseaseSimulation(t_step= h, **diseaseKwargs)
        ds_h()
        max_vals = ds_h.max_infected
        
        max_values[h] = max_vals
        if ((abs(max_vals[0]-max_vals_prev[0])/max_vals[0] > tolerance)
             or (abs(max_vals[1]-max_vals_prev[1])/max_vals[1] > tolerance)):
            max_vals_prev = max_vals
        else:
            if i==0:
                max_vals_prev = max_vals
                continue 
#             print("convergence achieved for h[{}]/2^{} with tolerance: {}%"
#                   .format(h_max, i, 100*tolerance))
#             for h, max in max_values.items():
#                 print(f"max(h={h:6.6f}) = t:{max[0]:8.2f}  max_infected: {max[1]:8.1f}")
            optimiced = True
            #Return optimized step
            return h
            break
    
    if not optimiced:
        print("WARNING: Method has not reach the [{}]% tolerance criteria."
              .format(100*tolerance))
    
    return h
    #return max_values


#===============================================================================
#     MODEL PARAMETER FITTER
#===============================================================================
def toleranceAchieved(value_ini, value_post, tolerance = .05):
    return abs(value_ini - value_post)/value_post < tolerance

def paramStep(data_t, model_t, param_value):
    step = param_value * (data_t - model_t) / (max(data_t, model_t)**1.2)
    return param_value + step

MAX_STEP = 50

def modelOptimizerFromData(N_population, 
                           parameters0, 
                           data, 
                           h_max=0.01, 
                           h_tolerance=0.1, 
                           data_tolerance=0.1):
    """From a set of data, find the best parameters. Optimize h in each step
    Args:
    :h_max = 0.01
    :parameters0 <dict> ={contagious_rate, recovery_rate, mortality_rate}
    :data <list of tuples> = [(day, infected, recovered, dead)]
    
    Return:
    <tuple> The most upgraded parameter sets and time step achieved
    """
    
    # TODO: Many variables could be grouped, avoiding single purpose definitions
    # TODO: Refactor in simple functions, excessive extension and cumbersome
    # TODO: Generalize the model for a general number of recoveries
    
    aux_params = {**parameters0,
                  't_step': h_max,
                  'days': 200, 
                  'N_population': N_population}
    post_params = {}
    iniArgs = {'day_0' : data[0][0], 'infected_0' : data[0][1], 
               'dead_0': data[0][3], 'recovered_0': data[0][2]}
    # Set up the first elements for t, infect, ... with the first data row
    DiseaseSimulation.setInitializers(**iniArgs)
    DiseaseSimulation.stopWhenMaxInfectedReached(True)
    DiseaseSimulation.setLogsPrint(False)
    
    dataTolAchieved = [dict([(key, False) for key in parameters0])
                       for _ in data]
    tolAchieved = [False for _ in data]
    ITER = 0
    evolutionParams = []
    keys = ('recovery_rate', 'mortality_rate', 'contagious_rate')
    while ((ITER < MAX_STEP) and (False in tolAchieved)):
        print('ITER Opt:'+str(ITER))
        ITER += 1
        # Adapt the time step (stepOptimizer)
        del aux_params['t_step']
        t_step = stepOptimizer(h_max, tolerance=h_tolerance, **aux_params)
        # Calculate the difference between the data and the model(with t_step optimized)
        # in t_data time. For equations of
        #          Death->M,  Recovered->RR,  (RR, M ,N, Infected)->CR
        aux_params['t_step'] = t_step
        model = DiseaseSimulation(**aux_params)
        model()
        # Calculate the step for each parameter as a difference data normalized 
        # by the minimum(difference).
        result_zip = model.getResults()
        #result_zip = (day, SUSCEPTIBLE, INFECTED, RECOVERED, DEAD)
        aux_rates = [dict([(key, 0) for key in keys]) for i in data]
        for i in range(len(data)):
            tuple_data = data[i]
            day = tuple_data[0]
            infc, reco, dead = tuple_data[1], tuple_data[2], tuple_data[3]
            
            _elem = 0
            for element in result_zip:
                if round(abs(element[0] - day), 6) > t_step:
                    _elem += 1
                    continue
                post_params[keys[0]] = paramStep(reco, element[3], aux_params[keys[0]])
                post_params[keys[1]] = paramStep(dead, element[4], aux_params[keys[1]])
                post_params[keys[2]] = paramStep(infc, element[2], aux_params[keys[2]]) 
                
                for k in range(3):
                    aux_rates[i][keys[k]] = post_params[keys[k]]
                break
            # check an arbitrary value of tolerance, return the parameters and t_step 
            # if it's exceeded
            for key in parameters0:
                dataTolAchieved[i][key] = toleranceAchieved(aux_params[key],
                                                            post_params[key],
                                                            data_tolerance)
            
            if False in dataTolAchieved[i].values():
                aux_params = {**aux_params, **post_params}
                tolAchieved[i] = False
        # Refresh the values with the median
        for k in range(3):
            aux_params[keys[k]] = np.median([rate[keys[k]] for rate in aux_rates])
        
        evolutionParams.append(aux_params)
        
    DiseaseSimulation.stopWhenMaxInfectedReached(False)
    graphEvolutionAndResultantModel(evolutionParams, aux_params)
    
    return aux_params


def graphEvolutionAndResultantModel(evolutionParams, finalParams):
    cRates = [cc['contagious_rate'] for cc in evolutionParams]
    rRates = [cc['recovery_rate'] for cc in evolutionParams]
    mRates = [cc['mortality_rate'] for cc in evolutionParams]
    
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3)
    plt.suptitle('Evolution of parameters with modelOptimizerFromData')
      
    axs[0].set_ylabel('r. contagious')
    axs[0].plot(cRates, 'rv-')
    axs[1].set_ylabel('r. recovery')
    axs[1].plot(rRates, 'g^-')
    axs[2].set_ylabel('r. mortality')
    axs[2].plot(mRates, 'k*-')
     
    #===========================================================================
    #     PRINT MODEL WITH RESULTS
    #===========================================================================
    DiseaseSimulation.setLogsPrint(True)
    model = DiseaseSimulation(**finalParams)
    model()
    print(model)
    model.graph()
