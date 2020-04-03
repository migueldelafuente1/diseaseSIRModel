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
    # BASIC RANGES AND EXAMPLES
    #===========================================================================
    ## Let be the parameters: (contagious_rate, recovery_rate) [1/day*people]
    params = [(0.265, 0.05), (2.5, 10.)]
    params_r_const = [(0.1 + 0.3*r, 20.) for r in range(5)]
    params_c_const = [(1.25, 1./(2.1 + 4.*c)) for c in range(0, 6)]
    
    print("CR\tRR\tmax_infected")
    DiseaseSimulation.setInitializers(infected_0=10000, dead_0=2, recovered_0=80)
    for param in params_c_const:
        ds = DiseaseSimulation(contagious_rate=param[0], 
                               recovery_rate=param[1],
                               t_step=0.005,
                               mortality_rate=0.04
                                )
        max_infected_analytical = ds.analyticalValueOfMaximumInfected()
        ds()
        ds.graph(details=False)
        ds.getDetails()
        print(f"{ds.CONT_RATE}\t{ds.RECO_RATE}\t{ds.max_infected}\t({max_infected_analytical})")
         
    #===========================================================================
    # OPTIMIZATION OF STEP FOR FIXED PARAMETERS (Example)
    #===========================================================================
    params = {'N_population': 6550000,
              'contagious_rate': 0.2759, 
              'recovery_rate': 0.0336, 
              'mortality_rate': 0.0268, 
              'days': 500,
              't_step': 0.01}
    DiseaseSimulation.setInitializers(day_0=17, infected_0=1990, 
                                      dead_0=81, recovered_0=1)
    #maximum_values = stepOptimizer(0.01, tolerance=0.005, **params)
    model = DiseaseSimulation(**params)
    model()
    model.graph(grid=True)
    #===========================================================================
    # CALIBRATION FOR MADRID
    #===========================================================================
    N_Population = 6550000
    data = [(17, 1990, 1, 81),  #13-3-2020
            (20, 4165, 474, 255),
            (23, 6777, 498, 941),
            (26, 9702, 1899, 1022)]
    params = {'contagious_rate': 1.25, 
              'recovery_rate': 1./7, 
              'mortality_rate': 0.05}
    # initial parameters
    params = modelOptimizerFromData(N_Population, params, data, h_max=0.01, 
                                    h_tolerance=0.05, data_tolerance=0.05)
    print(params)
    
        
    
    
    