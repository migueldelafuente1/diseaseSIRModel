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
    TOP_N = 500000
    PRINT = True
    
    # Initial Values.
    DAY_0 = 0
    INITIALIZERS = {'infected_0'  : 1,
                    'dead_0'      : 0,
                    'recovered_0' : 0
                    }
    
    def __init__(self,
                 t_step  = 0.01, # days
                 days    = 200,  # days
                 N_population    = 200000, # persons
                 contagious_rate = 0.0,    # persons/day
                 recovery_rate   = 0.0,    # person/day (= time to recovery**-1)
                 mortality_rate  = 0.0     # average mortality for the infected
                 ):
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
        
        self.time = None
        self.susceptible = [N_population - sum(self.INITIALIZERS.values())]
        self.infected    = [self.INITIALIZERS['infected_0']]
        self.recovered   = [self.INITIALIZERS['recovered_0']]
        self.dead        = [self.INITIALIZERS['dead_0']]
        self._setCalculationVars()
        
        self.CONT_FACTOR = contagious_rate  # Dimensionless Transmissibility
        self.CONT_RATE = contagious_rate/N_population # Rate for the calculations
        self.RECO_RATE = recovery_rate
        self.MORTALITY = mortality_rate
        self.INFECTEDS_COULD_DIE = mortality_rate > 0.0 
        
        self.R_0 = self.CONT_FACTOR / (self.RECO_RATE + self.MORTALITY)
        self.R_effective = self.R_0 * self.susceptible[0] / self.N_population
        
        self.max_infected = None
        self.__defineDerivates()
        
    def _setCalculationVars(self):
        """ Setting attributes as lists"""
        for _, deriv in self.DERIVATES_1st.items():
            setattr(self, deriv, [])
    
    @classmethod
    def setInitializers(cls, day_0=0, infected_0=1, dead_0=0, recovered_0=0):
        """ modify any of the initial values of the model, also the starting day"""
        cls.DAY_0 = day_0
        cls.INITIALIZERS = {'infected_0'  : infected_0,
                            'dead_0'      : dead_0,
                            'recovered_0' : recovered_0}
    
    def __repr__(self):

        return f"""
  =================================================================
    ***         DISEASE SIMULATION, SIR MODEL: INPUTS       ***
    -----------------------------------------------------------
       time step:\t{self.t_step}\t [days]
       days:     \t{self.days}\t [days]
       N_population:\t{self.N_population}\t [persons]
       
       contagious rate:\t{self.CONT_RATE*self.N_population:6.4f} [1/ persons day]
       recovery rate:  \t{self.RECO_RATE:6.4f} [1/ days]
       
       Considering People Die: {self.INFECTEDS_COULD_DIE} 
           Mortality: {self.MORTALITY*100:6.4f} %
       Considering [{self.GROUPS_OF_RECOVERY}] groups of recovery
       
       R0: {self.R_0:6.4f} \t R_effective: {self.R_effective:6.4f}
       Analytical Max Infections: {int(self.analyticalValueOfMaximumInfected())}
       
  =================================================================
"""
    # =========================================================================
    #   ODE SYSTEM
    # =========================================================================
    SUSCEPTIBLE = 'susceptible'
    INFECTED    = 'infected'
    RECOVERED   = 'recovered'
    DEAD  = 'dead'
    
    VARS = [SUSCEPTIBLE, INFECTED, RECOVERED, DEAD]
    
    VAR_COLORS= {SUSCEPTIBLE: 'b', 
                 INFECTED   : 'r',
                 RECOVERED  : 'g',
                 DEAD       : 'k'}
    
    DERIVATES_1st = {SUSCEPTIBLE : 'd_susceptible',
                     INFECTED    : 'd_infected',
                     RECOVERED   : 'd_recovered',
                     DEAD        : 'd_dead'}
    
    @classmethod
    def setLogsPrint(cls, bool_value):
        assert bool_value in (True, False), "PRINT is boolean, got '{}'".format(bool_value)
        cls.PRINT = bool_value
        
    def _logPrint(self, *msgs):
        if self.PRINT == True:
            print(*msgs)
            
    def getVariableTuple(self, i=-1):
        return tuple([getattr(self, var)[i] for var in self.VARS])
    
    def getResults(self):
        """ Get all the set of results for other uses:
        list of tuples: (day, SUSCEPTIBLE, INFECTED, RECOVERED, DEAD)"""
        if not self._converged:
            self._logPrint("WARNING: DiseaseModel has not converged")
        return zip(*([self.time] + [getattr(self, var) for var in self.VARS]))
    
    def __defineDerivates(self):
        
        self._derivates = {}
        
        self._derivates[self.SUSCEPTIBLE]= lambda s,i,r,d: -self.CONT_RATE*s*i
        self._derivates[self.INFECTED]   = lambda s,i,r,d: (self.CONT_RATE*s*i)\
             - ((self.RECO_RATE + self.MORTALITY)*i)
            
        self._derivates[self.RECOVERED]  = lambda s,i,r,d: self.RECO_RATE*i
        self._derivates[self.DEAD]       = lambda s,i,r,d: self.MORTALITY*i
        
    def __euler(self, i):
        new = {}
        for var_name, eq in self._derivates.items():
            
            new[var_name] = getattr(self, var_name)[-1]
            step = (self.t_step) * eq(*self.getVariableTuple())
            
            # Record the the difference
            getattr(self, self.DERIVATES_1st[var_name]).append(step)
            
            new[var_name] += step
            new[var_name] = max(min(self.N_population, new[var_name]), 0.01)
            
            getattr(self, var_name).append(new[var_name])
            # Define the maximum
            if (var_name == self.INFECTED):
                if (not self.max_infected):
                    if (step < 0) and (self.infected[-1] > self.INITIALIZERS['infected_0']):
                        self.max_infected = (round(self.DAY_0+self.t_step*i, 2),
                                             round(new[var_name]))
                        if self.STOP_WHEN_MAX_INFECTED_ACHIEVED:
                            self._logPrint("STOPPED IN MAX_INFECTED")
                            self._converged = True
                            
                else:
                    # If there is less than a person (after reaching the maximum
                    # of infections), the disease has been eradicated.
                    if (step < 0) and (getattr(self, self.INFECTED)[-1] < 1):
                        self._logPrint("CONVERGENCE ACHIEVED [step {}]".format(i))
                        self._converged = True
                        
    
    def __call__(self):
        """ run the execution for the object inputs """
        self._converged = False
        iterations, ini_step = 1, 0
        while not self._converged:
            self._logPrint(f"Iter {iterations}")
            for i in range(ini_step, (iterations * self.N_steps)-1):
                if self._converged:
                    len_ = len(self.infected)
                    self.time = [self.DAY_0+ j*self.t_step for j in range(len_)]
                    break
                self.__euler(i)
            if i >= self.TOP_N:
                self._logPrint("WARNING: MAX ITERATIONS reached, Convergence NOT ACHIEVED")
                self.time = [self.DAY_0+ j*self.t_step 
                             for j in range(len(self.infected))]
                break
            # Loop again if the process has not reached convergence
            ini_step = self.N_steps * iterations
            iterations += 1
        
    
    GRAPH_LABEL = 0
    @classmethod
    def graphLabelIncrement(cls):
        cls.GRAPH_LABEL += 1
    def getDetails(self):
        print(self)
    
    @classmethod
    def stopWhenMaxInfectedReached(cls, stop=True):
        """ Call this method to set stop the execution when the maximum of 
        infected is reached. Use it before the execution. """
        cls.STOP_WHEN_MAX_INFECTED_ACHIEVED = stop
       
    def analyticalValueOfMaximumInfected(self):
        """ Analytic value for the maximum infected population for the model. """
        _C = self.N_population / self.R_0
        return self.susceptible[0] + self.infected[0]\
            - (_C * (1 + np.log(self.susceptible[0] / _C)))
    
    def graph(self, details=True, logY=False, grid=True):
        """ Graph the results using matplotlib, also prints the object inputs. """
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
        if grid:
            ax.grid()
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
    
            
