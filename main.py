'''
Created on 18 mar. 2020

@author: Miguel
'''
from disease import DiseaseSimulation
from optimizers import stepOptimizer, modelOptimizerFromData

# =============================================================================
#   RUN
# =============================================================================

if __name__ == '__main__':
    #===========================================================================
    # BASIC RANGE
    #===========================================================================
    ## Let be the parameters: (contagious_rate, recovery_rate) [1/day*people]
#     params = [(0.265, 0.05), (2.5, 10.)]
#     params_r_const = [(0.1 + 0.3*r, 20.) for r in range(5)]
#     params_c_const = [(1.25, 1./(2.1 + 4.*c)) for c in range(0, 6)]
#     params_t_var   = [(0.05, 6.0, 0.01 + 0.01*t) for t in range(0, 10, 2)]
#      
#     print("CR\tRR\tmax_infected")
#     DiseaseSimulation.setInitializers(infected_0=10000, dead_0=2, recovered_0=80)
#     for param in params_c_const:
#         ds = DiseaseSimulation(contagious_rate=param[0], 
#                                recovery_rate=param[1]
#                                #,t_step=param[2]
#                                ,mortality_rate=0.04
#                                )
#         ds()
#         ds.graph(details=False)
#         #ds.getDetails()
#         print(f"{ds.CONT_RATE}\t{ds.RECO_RATE}\t{ds.max_infected}")
#         
    #===========================================================================
    # OPTIMIZATION OF STEP FOR FIXED PARAMETERS (Example)
    #===========================================================================
#    params = {'contagious_rate': 1.25, 'recovery_rate': 1/14, 'days': 200}
#    maximum_values = stepOptimizer(0.01, tolerance=0.0005, **params)
#    
    #===========================================================================
    # CALIBRATION FOR MADRID
    #===========================================================================
    N_Population = 6550000
    data = [(17, 1990, 1, 81),  #13-3
            (20, 4165, 474, 255),
            (23, 6777, 498, 941),
            (26, 9702, 1899, 1022)]
    params0 = {'CONT_RATE': 1.25, 'RECO_RATE': 1/7, 'MORTALITY': 0.05}
    # initial parameters
    modelOptimizerFromData(N_Population, params0, data, h_max=0.01, day0=17)
    