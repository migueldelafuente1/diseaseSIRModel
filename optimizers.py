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
#             print('t:', round(abs(max_vals[0]-max_vals_prev[0])/max_vals[0], 4),
#                   '\t i:', round(abs(max_vals[1]-max_vals_prev[1])/max_vals[1],4))
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
            break
    
    if not optimiced:
        print("WARNING: Method has not reach the [{}]% tolerance criteria."
              .format(100*tolerance))
    return max_values
                    
def modelOptimizerFromData():
    pass
    