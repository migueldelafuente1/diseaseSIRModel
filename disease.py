# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 11:08:42 2020

@author: Miguel
"""

import numpy as np

# =============================================================================
#   ODE SOLVERS
# =============================================================================
def euler(x_i, f_i, h):
    return h*x_i + f_i

# Python program to implement Runge Kutta method 
# A sample differential equation "dy / dx = (x - y)/2" 
def dydx(x, y): 
    return ((x - y)/2) 

def rungeKutta4(x0, y0, x, h): 
    # Count number of iterations using step size or 
    # step height h 
    n = (int)((x - x0)/h)  
    # Iterate for number of iterations 
    y = y0 
    for i in range(1, n + 1): 
        k1 = h * dydx(x0, y) 
        k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1) 
        k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2) 
        k4 = h * dydx(x0 + h, y + k3) 
  
        # Update next value of y 
        y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) 
  
        # Update next value of x 
        x0 = x0 + h 
    return y 
    
        

class DiseaseSimulation(object):
    # =========================================================================
    #   DATA INPUT
    # =========================================================================
    
    GROUPS_OF_RECOVERY = 1
    INFECTEDS_COULD_DIE = False
    
    def __init__(self,
                 t_step  = 0.01,
                 days    = 100,
                 N_population    = 1000,
                 contagious_rate = 0.0,
                 recovery_rate   = 0.0):
        
        self.t_step  = t_step #min(t_step, days * contagious_rate * recovery_rate / 1e5)
        self.days    = days
        self.N_steps = int(days / t_step)
        self.N_population = N_population
        
        self.time = np.arange(.0, days, t_step)
        self.susceptible = np.zeros(self.N_steps)
        self.infected    = np.zeros(self.N_steps)
        self.recovered   = np.zeros(self.N_steps)
        
        self.d_susceptible = np.zeros(self.N_steps)
        self.d_infected    = np.zeros(self.N_steps)
        self.d_recovered   = np.zeros(self.N_steps)
        
        self.susceptible[0] = N_population
        self.infected[0]    = 1
        
        self.CONT_RATE = round(contagious_rate, 4)
        self.RECO_RATE = round(recovery_rate, 4)
        
        self.max_infected = None
    
    def __repr__(self):

        return f"""
  =================================================================
    ***         DISEASE SIMULATION, SIR MODEL: INPUTS       ***
    -----------------------------------------------------------
       time step:\t{self.t_step}\t [days]
       days:     \t{self.days}\t [days]
       N_population:\t{self.N_population}\t [persons]
       
       contagious rate:\t{self.CONT_RATE} [persons/day]
       recovery rate:  \t{self.RECO_RATE} [persons/day]
       
       Considering People Die: {self.INFECTEDS_COULD_DIE}
       Considering [{self.GROUPS_OF_RECOVERY}] groups of recovey
       
  =================================================================
"""
    # =========================================================================
    #   ODE SYSTEM
    # =========================================================================
    SUSCEPTIBLE = 'susceptible'
    INFECTED    = 'infected'
    RECOVERED   = 'recovered'
    
    DERIVATES_1st = {SUSCEPTIBLE : 'd_susceptible',
                     INFECTED    : 'd_infected',
                     RECOVERED   : 'd_recovered'}
    
    VAR_INDEX = {SUSCEPTIBLE: 0, 
                 INFECTED   : 1,
                 RECOVERED  : 2}
    
    VAR_COLORS= {SUSCEPTIBLE: 'b', 
                 INFECTED   : 'r',
                 RECOVERED  : 'g'}
    
    def getVariables(self, i):
        return (self.susceptible[i], self.infected[i], self.recovered[i])
    
    def getResults(self):
        """ Get all the set of results for other uses"""
        return (self.susceptible, self.infected, self.recovered)
    
    def _derivates(self, sus_i, inf_i, rec_i, variable=None):
        if variable == self.SUSCEPTIBLE:
            return -1 * self.CONT_RATE * sus_i * inf_i
        
        elif variable == self.INFECTED:
            return (self.CONT_RATE * sus_i * inf_i) - (self.RECO_RATE * inf_i)
        
        elif variable == self.RECOVERED:
            return self.RECO_RATE * inf_i
        
        else:
            raise Exception(f"Invalid variable given [{variable}]")
            
    def _euler(self, *vars_, variable=None):
        """ Euler method for ODE solving """
        i_step = (self.t_step**2) * self._derivates(*vars_, variable)
        val = max(min(self.N_population, 
                      i_step + vars_[self.VAR_INDEX[variable]]), 0)
        return val, i_step
    
    def __call__(self):
        """ run the execution for the object inputs """
        for i in range(self.N_steps-1):
            vars_ = self.getVariables(i)
            for _variable in self.VAR_INDEX.keys():
                value, step = self._euler(*vars_, variable=_variable)
                
                if _variable == self.SUSCEPTIBLE:
                    self.susceptible[i+1], self.d_susceptible[i] = value, step
                elif _variable == self.INFECTED:
                    self.infected[i+1], self.d_infected[i] = value, step
                    if (step < 0) and (not self.max_infected):
                        self.max_infected = (round(self.time[i],2), round(value))
                elif _variable == self.RECOVERED:
                    self.recovered[i+1], self.d_recovered[i] = value, step 
    
    GRAPH_LABEL = 0
    @classmethod
    def graphLabelIncrement(cls):
        cls.GRAPH_LABEL += 1
        
    def graph(self, details=True, logY=False):
        """ Graph the results using matplotlib, also prints the object inputs"""
        if details:
            print(self)
        
        import matplotlib.pyplot as plt
        #gn = f"Graph {self.GRAPH_LABEL}"
        self.graphLabelIncrement()
        
        fig, ax = plt.subplots()
        for _var in self.VAR_INDEX.keys():
            ax.plot(self.time, 
                    getattr(self, _var),
                    self.VAR_COLORS[_var],
                    label=_var)
        
        ax.set_title("Time evolution of Disease\n", loc='center', fontsize=18)
        ax.set_title(f"CR: {self.CONT_RATE} RR: {self.RECO_RATE} [persons/day] Max={self.max_infected}", 
                     loc='left', fontsize=13, color='grey')
        ax.set_xlabel("time [days]")
        ax.set_ylabel("People")
        ax.legend()
        
        if logY:
            plt.semilogy()
        plt.show()
        


# =============================================================================
#   OBTAIN CONTAGIOUS COEFICIENT (data since 12-2-20)
# =============================================================================
from datetime import datetime, timedelta 
def getContagiousRate(detected_list):
    t, detected = zip(*detected_list)
    x = np.vstack([np.array(t), np.ones(len(t))]).T
    detected = list(map(lambda d: np.log(d), detected))

    A, B = np.linalg.lstsq(x, detected, rcond=None)[0]
    return round(A, 4), round(np.exp(B), 4)

t_0 = datetime(2020, 2, 12)

detected_list = [(0,  2),    (16, 31), 
                 (23, 365),  (27, 1639), 
                 (28, 2140), (29, 2965), 
                 (30, 4231), (31, 5753)]

i_rate, inf_y_abs = getContagiousRate(detected_list)
t = 33
print(f"IR:{i_rate} abscI:{inf_y_abs}")
_infected = round(inf_y_abs*np.exp(i_rate*t))
print(f"{(t_0+timedelta(days=t)).date()} : infected= {_infected}")

# =============================================================================
#   RUN
# =============================================================================

if __name__ == '__main__':
    ## Let be the parameters: (contagious_rate, recovery_rate) [people/day]
    params = [(0.265, 0.05), (2.5, 10.)]
    params_r_const = [(0.1 + 0.3*r, 20.) for r in range(5)]
    params_c_const = [(0.05, 2.1 + 4.*c) for c in range(6)]
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
        print(f"{ds.CONT_RATE}\t{ds.RECO_RATE}\t{ds.max_infected}")
