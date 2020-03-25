'''
Created on 20 mar. 2020

@author: Miguel
'''
from disease import DiseaseSimulation

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
#     tol_ = 1 - tolerance
    max_vals_prev = (0,0)
    max_values = {}
    optimiced = False
    for i in range(7):
        h = h_max / (2**(i))
        
        ds_h = DiseaseSimulation(t_step= h, **diseaseKwargs)
        ds_h.stopWhenMaxInfectedReached(True)
        ds_h()
        max_vals = ds_h.max_infected
        
        max_values[h] = max_vals
        print('t:', round(abs(max_vals[0]-max_vals_prev[0])/max_vals[0], 4),
                  '\t i:', round(abs(max_vals[1]-max_vals_prev[1])/max_vals[1],4))
        if ((abs(max_vals[0]-max_vals_prev[0])/max_vals[0] > tolerance) or
            (abs(max_vals[1]-max_vals_prev[1])/max_vals[1] > tolerance)):
            max_vals_prev = max_vals
        else:
            if i==0:
                max_vals_prev = max_vals
                continue 
            print("convergence achieved for h[{}]/2^{} with tolerance: {}%"
                  .format(h_max, i, 100*tolerance))
            for h, max in max_values.items():
                print(f"max(h={h:6.6f}) = t:{max[0]:8.2f}  max_infected: {max[1]:8.1f}")
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
def modelOptimizerFromData(N_Population, parameters0, data, h_max=0.01, day0=0):
    """From a set of data, find the best parameters. Optimize h in each step
    Args:
    :h_max = 0.01
    :parameters0 <dict> ={CONTAGIOUS_RATE, RECOVERY_RATE, MORTALITY}
    :data <list of tuples> = [(day, infected, recovered, dead)]
    
    Return:
    <tuple> The most upgraded parameter sets and time step achieved"""
    
    t, infc, reco, dead = zip(data)
    

    
    
    