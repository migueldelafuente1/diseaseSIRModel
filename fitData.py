'''
Created on 18 mar. 2020

@author: Miguel
'''

# =============================================================================
#   OBTAIN CONTAGIOUS COEFICIENT (data since 12-2-20)
# =============================================================================
from datetime import datetime, timedelta 
import numpy as np
import matplotlib.pyplot as plt

def getContagiousRate(detected_list):
    """ 
    Args:
    :detected_list with tuples:(day, infected, dead)
    
    Returns:
    :tuple 
        with exponential Slope + origin for infected,
        with exponential Slope + origin for deaths
        with the mortality and its stdv"""
    t, infected, dead = zip(*detected_list)
    x = np.vstack([np.array(t), np.ones(len(t))]).T
    infected = list(map(lambda i: np.log10(i), infected))
    dead     = list(map(lambda d: np.log10(d+0.1), dead))

    A_i, B_i = np.linalg.lstsq(x, infected, rcond=None)[0]
    A_d, B_d = np.linalg.lstsq(x, dead, rcond=None)[0]
    
    mortality = list(map(lambda vals: vals[2]/vals[1], detected_list))
    return (
        round(A_i, 4), round(10**(B_i), 4),
        round(A_d, 4), round(10**(B_d), 4), 
        np.mean(mortality), np.std(mortality))



#===============================================================================
## rtve mapa-coronavirus
## t_0 = datetime(2020, 2, 12)
#===============================================================================
#spain_detected = [#(0, 2, 0), 
#                  (16, 31, 0), 
#                  (23, 365, 5), 
#                  (27, 1639, 36), 
#                  (28, 2140, 48), 
#                  (29, 2965, 189), 
#                  (30, 4231, 193), 
#                  (31, 5753, 136),
#                  (32, 7753, 288), 
#                  (33, 9191, 309),
#                  (34, 11178, 491), 
#                  (35, 13716, 598)]
#
#===============================================================================
## from worldmeters.info
## day 0 is 15/2/2020
#===============================================================================
""" TODO: consider inserting a bot for this page:
    click on points: <tabbable-line-cases><div tab-contnet>... inspect to 
         class="highcharts-markers highcharts-series-0 highcharts-line-series
    text and date in: //*[@id="highcharts-9wyla2t-118"]/svg/g[9]/text" 
"""
t_0 = datetime(2020, 2, 15)

spain_detected = [
                  (8, 2, 0),
                  (12, 25, 0),
                  (14, 58, 0),
                  (16, 120, 0),# 2-3-20
                  (17, 165, 1),
                  (19, 282, 3),
                  (21, 525, 10),
                  (22, 674, 17), # 8-3-20
                  (23, 1231, 30),
                  (24, 1695, 36),
                  (25, 2277, 55),
                  (26, 3146, 86),
                  (27, 5232, 133),
                  (28, 6391, 196),
                  (29, 7988, 294), # 15-3-20
                  (30, 9942, 342),
                  (31, 11826, 533),
                  (32, 14769, 638),
                  (33, 18077, 831),
                  (34, 21571, 1093),
                  (35, 25496, 1381),
                  (36, 28768, 1772),
                  (37, 35136, 2311),
                  (38, 42058, 2991),
                  (39, 49515, 3647), # 25-3-20
                  (40, 57786, 4365),
                  (41, 65719, 5138),
                  (42, 73235, 5282)
                  ]

italy_detected = [
                  (0, 3, 0),
                  (6, 21, 1),
                  (9, 229, 7),
                  (11, 470, 12),
                  (13, 889, 21),
                  (15, 1701, 41), # 1-3-20
                  (17, 2502, 79),
                  (18, 3089, 107),
                  (19, 3858, 148),
                  (20, 4636, 197),
                  (21, 5883, 233),
                  (22, 7375, 366),
                  (23, 9172, 463),
                  (24, 10149, 631),
                  (25, 12462, 827),
                  (26, 15113, 1016),
                  (27, 17660, 1266),
                  (28, 21157, 1441),
                  (29, 24747, 1809), # 15-3-20
                  (30, 27980, 2158),
                  (31, 31506, 2503),
                  (32, 35713, 2978),
                  (33, 41035, 3405),
                  (34, 47021, 4032),
                  (35, 53578, 4825),
                  (36, 59138, 5476),
                  (37, 63927, 6077),
                  (38, 69176, 6820),
                  (39, 74386, 7503), # 25-3-20
                  (40, 80589, 8215),
                  (41, 86498, 9134),
                  (42, 92472, 10023),
                  ]


def graphDataAndFit(data_list, A_i, B_i, A_d, B_d, t_min, t_max, name):
    
    t, infected, dead = zip(*data_list)
    dates_t = [(t_0+timedelta(days=t_min + t)).date() for t in range(len(t))]
    fig, ax = plt.subplots()
    ax.plot(dates_t, infected, 'r^-',label=f'{name} infected')
    ax.plot(dates_t, dead, 'kv-',label=f'{name} dead')
    
    x_i, y_i = zip(*[((t_0+timedelta(days=t)).date(), B_i* (10**(A_i*t)))
                     for t in np.arange(t_min, t_max, 1.)])
    x_d, y_d = zip(*[((t_0+timedelta(days=t)).date(), B_d* (10**(A_d*t))) 
                     for t in np.arange(t_min, t_max, 1.)])
    ax.plot(x_i, y_i, 'r--',label=f'{name} infected')
    ax.plot(x_d, y_d, 'k--',label=f'{name} dead')
    plt.xticks(rotation=30)
    ax.legend()
    ax.grid()
    ax.semilogy()
    ax.set_xlabel(f'time (days) from {(t_0+timedelta(days=t_min)).date()}={t_min} to {(t_0+timedelta(days=t_max)).date()}={t_max}')

def getRatesForDateRanges(date_min, date_max, date_prediction=None, graph=True):
    """
    Main Function process here, get data range by dates, call the fitting 
    function and plot and predict for the approximation.
        Args:
    :date_min :date_max :date_prediction(optional) are datetime objects 
    :graph(=True) to graph selected data and trend.
        Return:
    dictionary (for each country) with a tuple:
    the infection and death rate and also the average mortality for the period.
    """
    t_min = (date_min - t_0).days
    t_max = (date_max - t_0).days
    
    dicts_ = {'spain': spain_detected, 
              'italy': italy_detected}
    
    constants = {}
    ## Set data in the range
    for name, detected_list in dicts_.items():
        i_min, i_max = None, 0
        for i in range(len(detected_list)):
            if (i == len(detected_list)-1):
                i_max = i
            if detected_list[i][0] >= t_min and i_min is None:
                i_min = i
            if (detected_list[i][0] >= t_max and not i_max):
                i_max = i
                break
        print(f"{name} -> i_min[{i_min}]  i_max[{i_max}]")
        dicts_[name] = detected_list[max(0, i_min) : min(i_max, len(detected_list))]

    for name, detected_list in dicts_.items():
        print(name.upper()
                +"   data from [{}] to [{}]".format(
                        date_min.strftime("%Y-%m-%d"),
                        date_max.strftime("%Y-%m-%d"))
                +"\n-----------------------------------------")
        vals = getContagiousRate(detected_list)
        i_rate, i_y_abs, d_rate, d_y_abs, mort, mort_std = vals
        print(f"Infect Rate: {i_rate:8.4f}   abscI: {i_y_abs:8.4f}")
        print(f"Death  Rate: {d_rate:8.4f}   abscD: {d_y_abs:8.4f}")
        print(f"Mortality:  {mort:8.4f} +/- {mort_std:6.4f}")
        
        constants[name] = (i_rate, d_rate, mort, mort_std)
        
        if date_prediction:
            t_pred = (date_prediction - t_0).days
            
            _infected = round(i_y_abs * (10**(i_rate*t_pred)))
            _deaths   = round(d_y_abs * (10**(d_rate*t_pred)))
            print(f"\t{(t_0+timedelta(days=t_pred)).date()} : infected= {_infected}")
            print(f"\t{(t_0+timedelta(days=t_pred)).date()} : death= {_deaths}\n")
        if graph:
            graphDataAndFit(detected_list, i_rate, i_y_abs, 
                            d_rate, d_y_abs, t_min, t_max, name)
    
    return constants

def graphRatesEvolution(window_infections, window_deaths, window_mortality,
                        date_min, window_delta_days):
    time = [date_min+timedelta(days=i) for i in range(len(window_infections))]
    #time = [i for i in range(len(window_infections))]
    
    fig1, axs = plt.subplots(2)
    fig1.suptitle(f'Results of fit for infectious and death rate. \n(+- {window_delta_days} days window)')
    axs[0].plot(time, window_infections, 'o-r', label='infectious rate')
    axs[1].plot(time, window_deaths, 'x-k', label='death rate')
    axs[1].set_xlabel(f'time (days) from {date_min.date()}=0')
    axs[0].set_title("Infection Rate")
    axs[1].set_title("Death Rate")
    
    for tick in axs[0].get_xticklabels():
        tick.set_rotation(15)
    for tick in axs[1].get_xticklabels():
        tick.set_rotation(15)
    plt.legend(loc='lower right')
    
    fig2 = plt.figure()
    fig2.suptitle('Mortality averaged.')
    mort  = list(map(lambda x: x[1], window_mortality))
    m_err = list(map(lambda x: x[1], window_mortality))
    plt.errorbar(time, mort, yerr=m_err, label='mortality average')
    plt.legend(loc='lower right')
    plt.xticks(rotation=45)

# =============================================================================
#       MAIN RUN
# =============================================================================
from copy import copy
if __name__ == '__main__':
     
    ## INPUTS
    # =========================================================================
     
    date_min = datetime(2020, 3, 1)        # lower bound of data range
    date_max = datetime(2020, 3, 28)        # upper bound of data range
    date_prediction = datetime(2020, 3, 28) # prediction 
    # =========================================================================
    getRatesForDateRanges(date_min, date_max, date_prediction, graph=True)
    
    
    #===========================================================================
    #     ITERATION OF RANGE WITH delta_days WINDOW
    #===========================================================================
    date_min = datetime(2020, 3, 1)        # lower bound of data range
    date_max = datetime(2020, 3, 28)        # upper bound of data range
     
    delta_days = 1 # window of days for the mean( 2*delta_days + 1) 
    date_min_MIN = date_min + timedelta(days=delta_days)
    date_max_MAX = date_max - timedelta(days=delta_days)
     
    window_mean_infections = []
    window_mean_deaths = []
    window_mean_mortality = []
     
    date_min = copy(date_min_MIN)
    date_max = date_min + timedelta(days= 2*delta_days + 1)
     
    while(date_max <= date_max_MAX):
        results = getRatesForDateRanges(date_min, date_max, graph=False)
         
        window_mean_infections.append(results['spain'][0])
        window_mean_deaths.append(results['spain'][1])
        window_mean_mortality.append(results['spain'][2:4])
         
        date_min = date_min + timedelta(days=1)
        date_max = date_min + timedelta(days= 2*delta_days + 1)
     
    graphRatesEvolution(window_mean_infections, window_mean_deaths, 
                        window_mean_mortality, date_min_MIN, delta_days)
     
     