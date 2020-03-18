'''
Created on 18 mar. 2020

@author: Miguel
'''
from disease import DiseaseSimulation

# =============================================================================
#   RUN
# =============================================================================

if __name__ == '__main__':
    ## Let be the parameters: (contagious_rate, recovery_rate) [people/day]
    params = [(0.265, 0.05), (2.5, 10.)]
    params_r_const = [(0.1 + 0.3*r, 20.) for r in range(5)]
    params_c_const = [(1.25, 1./(2.1 + 4.*c)) for c in range(0, 6)]
    params_t_var   = [(0.05, 6.0, 0.01 + 0.01*t) for t in range(0, 10, 2)]
    
    print("CR\tRR\tmax_infected")
    for param in params_c_const:
        DiseaseSimulation()
        ds = DiseaseSimulation(contagious_rate=param[0], 
                               recovery_rate=param[1]
                               #,t_step=param[2]
                               )
        ds()
        ds.graph(details=False)
#        ds.getDetails()
        print(f"{ds.CONT_RATE}\t{ds.RECO_RATE}\t{ds.max_infected}")