# diseaseModel
Model of infection based on SIR model, given by this EDO system:
	
	ds/dt = -CONT_RATE*s*i/N_population
	di/dt = (self.CONT_RATE*s*i/N_population) - (self.RECO_RATE*i)
	dr/dt = self.RECO_RATE*i

Where `s` stands for _susceptibe_, `i` for _infected_ and `r` for _recovered_, `CONT_RATE` or __contagious rate__ is the average number of contacts between people for a person per day times the probability of infection per contact (needs to be normalized to the total population). The __recovery rate__ (`RECO_RATE`) here is the inverse average time to overcome the disease and stop infecting others (in other words, the factor for the number of recoveries per day).

The model could also be extended to take into account different rates of recovery or contagious, death according to mortality or other dependencies.

#Requirements:

python `3.6` or newer, `matplotlib` and `numpy`

	pip3 install matplotlib
	pip3 install numpy

# Parts:
* **fitData** script has a time compilation from some online sources of _corona virus_ infected and deaths for Italy and Spain. The module just have a linear fitting function for base 10 logarithmic increase. The resultant constants are used to plot the perspectives for a certain day (counting from _2020-2-15_) and a plot of this linear trend for both countries.

* **disease** module has the Disease Models. The models has the following parameters:
	| Parameter | Default | Units | Description |
	t_step  = 0.01, # days
    days    = 200,  # days
	N_population    = 200000, # persons
	contagious_rate = 0.0,    # persons/day
	recovery_rate   = 0.0
* **optimizers** has some functions to find the best parameters for the calculations.
* * _stepOptimizer_ find the step necessary for the Euler-Method iteration to be in a relative tolerance range when the time step is split. 
Args:
	h_max (in days), tolerance=0.2 (ratio 1 for 100%), **diseaseKwargs(fixed parameters for the model)

* Run disease model and optimizers rom **main** 