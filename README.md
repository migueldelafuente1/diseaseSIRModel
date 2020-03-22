# disease Model
Model of infection based on SIR model, given by this EDO system:
	
	ds/dt = -CONT_RATE*s*i/N_population
	di/dt = (self.CONT_RATE*s*i/N_population) - (self.RECO_RATE*i)
	dr/dt = self.RECO_RATE*i

$$|\vec{A}|=\sqrt{A_x^2 + A_y^2 + A_z^2}.$$(1)

Where `s` stands for _susceptibe_, `i` for _infected_ and `r` for _recovered_. `CONT_RATE` or __contagious rate__ is the average number of contacts between people for a person per day times the probability of infection per contact (needs to be normalized to the total population). The __recovery rate__ (`RECO_RATE`) here is the inverse average time to overcome the disease and stop infecting others (in other words, the factor for the number of recoveries per day).

The model could also be extended to take into account different rates of recovery or contagious, death according to mortality or other dependencies.

# Requirements:

python `3.6` or newer, `matplotlib` and `numpy`

	pip3 install matplotlib
	pip3 install numpy

# Parts:
* **fitData** script has a time compilation from some online sources of _corona virus_ infected and deaths, for Italy and Spain. The module has a linear fitting function for base 10 logarithmic increase. The resultant constants are used to plot the perspectives for a certain day (counting from _2020-2-15_) and a plot of this linear trend for both countries.

__Arguments__:

~~~~~~~~~~~~~{.py}
    ## INPUTS
    
    date_min = datetime(2020, 3, 5)         # lower bound of data range
    date_max = datetime(2020, 3, 11)        # upper bound of data range
    date_prediction = datetime(2020, 3, 12) # prediction
~~~~~~~~~~~~~
    
The output depends on the date range:

	SPAIN   data from [2020-03-05] to [2020-03-11]
	-----------------------------------------
	Infect Rate:   0.1582   abscI:   0.2600
	Death  Rate:   0.2208   abscD:   0.0002
	Mortality:    0.0201 +/- 0.0052
        2020-03-12 : infected= 3374.0
        2020-03-12 : death= 110.0

	ITALY   data from [2020-03-05] to [2020-03-11]
	-----------------------------------------
	Infect Rate:   0.0882   abscI:  81.8895
	Death  Rate:   0.1273   abscD:   0.5487
	Mortality:    0.0471 +/- 0.0082
        2020-03-12 : infected= 16085.0
        2020-03-12 : death= 1120.0

The prediction depends on the date range selected, notice that my data for February and the beginning of March is not daily.In the other hand, governmental actions affects on the global trend of the progression, which is neither a simple exponential nor a geometric progression of a constant factor. This approximations are reliable for short term data collections (f.e, a week or 5 days) and a couple of days in the future. The approximation won't be valid when the active infections were near to the maximum.

 
        
* **disease** module has the Disease Models. The models has the following parameters:

| Parameter | Default | Units | Description |
| --- | --- | --- | --- |
| t_step |  = 0.01 | days | Step for numeric method, optimize it with _optimizers.stepOptimizer_ tool. |
| days  |   = 200   | days |  |
| N_population | = 200000 | persons |  |
| contagious_rate | = 0.0 | 1/ day * persons |  |
| recovery_rate | = 0.0 | 1/ day |  |

* **optimizers** has some functions to find the best parameters for the calculations.
** _stepOptimizer_ find the step necessary for the Euler-Method iteration to be in a relative tolerance range when the time step is split.

This function iterates to find a value of the step for which the SIR simulation converges under a certain tolerance. The steps are divided by 2 in each iteration (for a maximum of h_max/2^-6).
    | Argument | type | Description |
    | __h_max__ | double | first step value. in days|
    | __tolerance__ | double |(=0.2 by default), tolerance ratio on step (1.0 for 100%).|
    | __diseaseKwargs__ |Dictionary | The parameters for the DiseaseSimulation. |
** _modelOptimizerFromData_
* Run disease model and optimizers from **main** 