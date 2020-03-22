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
    STOP_WHEN_MAX_INFECTED_ACHIEVED = False
    
    def __init__(self,
                 t_step  = 0.01, # days
                 days    = 200,  # days
                 N_population    = 200000, # persons
                 contagious_rate = 0.0,    # persons/day
                 recovery_rate   = 0.0):   # person/day (= time to recovery**-1)
        """ 
        Args:
        :contagious_rate = contacts/(day * person) *
            * transmission probability by contact
        :recovery_rate = 1 / days to recovery (period of infection)
        """
        
        self.t_step  = t_step #min(t_step, days * contagious_rate * recovery_rate / 1e5)
        self.days    = days
        self.N_steps = int(days / t_step)
        self.N_population = N_population
        
        self.time = np.arange(.0, days, t_step)
        self.susceptible = [N_population]
        self.infected    = [1]
        self.recovered   = [0]
        self._setCalculationVars()
        
        self.CONT_RATE = contagious_rate/N_population
        self.RECO_RATE = recovery_rate
        
        self.max_infected = None
        self.__defineDerivates()
    
    def _setCalculationVars(self):
        """ Setting attributes as lists"""
        for _, deriv in self.DERIVATES_1st.items():
            setattr(self, deriv, [])
        for _, err in self.ERROR_2ord.items():
            setattr(self, err, [])
    
    def __repr__(self):

        return f"""
  =================================================================
    ***         DISEASE SIMULATION, SIR MODEL: INPUTS       ***
    -----------------------------------------------------------
       time step:\t{self.t_step}\t [days]
       days:     \t{self.days}\t [days]
       N_population:\t{self.N_population}\t [persons]
       
       contagious rate:\t{self.CONT_RATE} [1/ persons day]
       recovery rate:  \t{self.RECO_RATE} [1/ days]
       
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
    
    VARS = [SUSCEPTIBLE, INFECTED, RECOVERED]
    
    VAR_COLORS= {SUSCEPTIBLE: 'b', 
                 INFECTED   : 'r',
                 RECOVERED  : 'g'}
    
    DERIVATES_1st = {SUSCEPTIBLE : 'd_susceptible',
                     INFECTED    : 'd_infected',
                     RECOVERED   : 'd_recovered'}
    ERROR_2ord    = {SUSCEPTIBLE : 'e_susceptible',
                     INFECTED    : 'e_infected',
                     RECOVERED   : 'e_recovered'}
    
    def getVariableTuple(self, i=-1):
        return tuple([getattr(self, var)[i] for var in self.VARS])
    
    def getResults(self):
        """ Get all the set of results for other uses"""
        return tuple([getattr(self, var) for var in self.VARS])
    
    def __defineDerivates(self):
        
        self._derivates = {}
        
        self._derivates[self.SUSCEPTIBLE]= lambda s,i,r: -self.CONT_RATE*s*i
        self._derivates[self.INFECTED] = lambda s,i,r: (self.CONT_RATE*s*i) - (self.RECO_RATE*i)
        self._derivates[self.RECOVERED]= lambda s,i,r: self.RECO_RATE*i
    
    def __eulerError(self, var_name, deriv):
        _error = 0.5 * (self.t_step**2) * deriv
        _lastErr = getattr(self, self.ERROR_2ord[var_name])[-1]
        getattr(self, self.ERROR_2ord[var_name]).append(_error + _lastErr)
        
    def __euler(self, i):
        new = {}
        for var_name, eq in self._derivates.items():
            
            new[var_name] = getattr(self, var_name)[-1]
            step = (self.t_step) * eq(*self.getVariableTuple())
            
            # Record the the difference
            getattr(self, self.DERIVATES_1st[var_name]).append(step)
#             self.__eulerError(var_name, step)
            
            new[var_name] += step
            new[var_name] = max(min(self.N_population, new[var_name]), 0.01)
            
            getattr(self, var_name).append(new[var_name])
            # Define the maximum
            if ((var_name == self.INFECTED) 
                and (step < 0) and (not self.max_infected)):
                self.max_infected = (round(self.time[i],2), round(new[var_name]))
                if self.STOP_WHEN_MAX_INFECTED_ACHIEVED:
                    break
        
    
    def __call__(self):
        """ run the execution for the object inputs """
        for i in range(self.N_steps-1):
            self.__euler(i)
                        
    
    GRAPH_LABEL = 0
    @classmethod
    def graphLabelIncrement(cls):
        cls.GRAPH_LABEL += 1
    def getDetails(self):
        print(self)
    
    def stopWhenMaxInfectedReached(self, stop=True):
        """ Call this method to set stop the execution when the maximum of 
        infected is reached. Use it before the execution. """
        self.STOP_WHEN_MAX_INFECTED_ACHIEVED = stop
    
    def graph(self, details=True, logY=False):
        """ Graph the results using matplotlib, also prints the object inputs"""
        if details:
            print(self)
        
        import matplotlib.pyplot as plt
        #gn = f"Graph {self.GRAPH_LABEL}"
        self.graphLabelIncrement()
        
        fig, ax = plt.subplots()
        for _var in self.VARS:
            ax.plot(self.time, 
                    getattr(self, _var),
                    self.VAR_COLORS[_var],
                    label=_var)
        
        ax.set_title("Time evolution of Disease\n", loc='center', fontsize=18)
        ax.set_title(f"CR={round(self.CONT_RATE*self.N_population, 4)}[1/pd] 1/RR={round(self.RECO_RATE**(-1), 2)} \
[d] Max={self.max_infected}[(d, p)]", 
                     loc='left', fontsize=13, color='grey')
        ax.set_xlabel("time [days]")
        ax.set_ylabel("People")
        ax.legend()
        
        if logY:
            plt.semilogy()
        plt.show()
        

