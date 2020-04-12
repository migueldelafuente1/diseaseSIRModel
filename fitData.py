'''
Created on 18 mar. 2020

@author: Miguel
'''

# =============================================================================
#   OBTAIN CONTAGIOUS COEFICIENT (data since 12-2-20)
# =============================================================================
from datetime import datetime, timedelta 
import numpy as np
from dataWebLoader import DataWebLoader
from copy import copy, deepcopy

def getContagiousRate(detected_list):
    """ 
    Args:
    :detected_list with tuples:(day, infected, dead)
    
    Returns:
    :tuple 
        with exponential Slope + origin for infected,
        with exponential Slope + origin for deaths
        with the mortality and its stdv"""
    tuple_ = list(zip(*detected_list))
    if len(tuple_) == 3:
        t, infected, dead = tuple_
    elif len(tuple_) == 4:
        t, infected, dead, active = tuple_
        
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
                  (42, 73235, 5282),
                  (43, 80110, 6803)
                  ]


spain_detected = [(1, 1990, 81),
                  (4, 4165, 255),
                  (7, 6777, 941),
                  (10, 9702, 1022)]

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
                  (43, 97689, 10779)
                  ]

COUNTRIES_DATA = {'spain': spain_detected, 
                  'italy': italy_detected}

def graphDataAndFit(data_list, A_i, B_i, A_d, B_d, t_min, t_max, name):
    import matplotlib.pyplot as plt
    
    tuple_ = list(zip(*data_list))
    if len(tuple_) == 3:
        t, infected, dead = tuple_
    elif len(tuple_) == 4:
        t, infected, dead, active = tuple_
        
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
    
    dicts_ = deepcopy(COUNTRIES_DATA)
    
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
        if i_min is None:
            return {}
        dicts_[name] = detected_list[max(0, i_min) : min(i_max+1, 
                                                         len(detected_list))]

    for name, detected_list in dicts_.items():
        if detected_list == []:
            continue
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
    import matplotlib.pyplot as plt
    #time = [date_min+timedelta(days=i) for i in range(len(window_infections))]
    time = [i for i in range(len(window_infections))]
    
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

def processWebData(country, data_dict):
    # keys from dataWebLoader.GRAPHS
    if country == DataWebLoader.TIMESTAMP:
        return None
    
    dates, totalInf = data_dict['Total Cases']
    _, totalDeaths  = data_dict['Total Deaths']
    _, totalActive  = data_dict['Active Cases']
    
    days = [(date - t_0.date()).days for date in dates]
    
    COUNTRIES_DATA[country] = list(zip(days, totalInf, totalDeaths, totalActive))

def graphTotalVsActiveCases(t_min, t_max, logXaxis=True):
    """ Graph logarithmically the total cases against the currently active cases
    Requires execute processWebData since theere is no active register locally 
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    
    for country, data_list in COUNTRIES_DATA.items():
        t, infected, dead, active = zip(*data_list)
        _color = np.random.rand(3,)
        ax.plot(infected, active, '.-', c=_color ,label=f'{country}', alpha=0.5)
        ax.plot(infected[-1], active[-1], marker='P',c=_color, mec='k', ms=7)
    
    ax.legend()
    ax.grid()
    plt.suptitle("Total Vs Active Cases Evolution", fontsize=12)
    plt.title(f'data from {t_min.date()} to {t_max.date()}')
    if logXaxis:
        ax.semilogx()
    ax.semilogy()
    ax.set_xlim(10, None)
    ax.set_ylim(10, None)
    ax.set_xlabel('Total Cases')
    ax.set_ylabel('Active Cases')

# =============================================================================
#       MAIN RUN
# =============================================================================

if __name__ == '__main__':
    
    ## INPUTS
    # =========================================================================
     
    date_min = datetime(2020, 3, 1)         # lower bound of data range
    date_max = datetime.now()               # upper bound of data range
    date_prediction = (datetime.now() + timedelta(days=1)) # prediction 

#    date_max = t_0 + timedelta(days=11)              # upper bound of data range
    date_prediction = date_max + timedelta(days=2) # prediction 
    # =========================================================================
#    getRatesForDateRanges(date_min, date_max, date_prediction, graph=True)
    
    #===========================================================================
    #     LOAD DATA FROM WEB
    #===========================================================================
    
    from dataWebLoader import DataWebLoader
    
    countries = ['spain', 'italy', 'germany', 'us', 'uk', 'south-korea']
    
    data = DataWebLoader.getData(countries, save=True)
    
    for country, tables in data.items():
        processWebData(country, tables)
    
    getRatesForDateRanges(date_min, date_max, date_prediction, graph=True)
    
    ## TOTAL CASES VS ACTIVE CASES GRAPH
    #---------------------------------------------------------------------------
    graphTotalVsActiveCases(date_min, date_max)
    
    #===========================================================================
    #     ITERATION OF RANGE WITH delta_days WINDOW (uncomment)
    #===========================================================================
#
#    date_min = datetime(2020, 3, 1)         # lower bound of data range
#    date_max = datetime(2020, 3, 28)        # upper bound of data range
#    
#    countries = ['spain', 'italy']
#    
#    data = DataWebLoader.getData(countries, save=True)
#    for country, tables in data.items():
#        processWebData(country, tables)
#        
#    delta_days = 1 # window of days for the mean( 2*delta_days + 1) 
#    date_min_MIN = date_min + timedelta(days=delta_days)
#    date_max_MAX = date_max - timedelta(days=delta_days)
#    
#    for  country in countries:
#        window_mean_infections = []
#        window_mean_deaths = []
#        window_mean_mortality = []
#         
#        date_min = copy(date_min_MIN)
#        date_max = date_min + timedelta(days= 2*delta_days + 1)
#        
#        while(date_max <= date_max_MAX):
#            results = getRatesForDateRanges(date_min, date_max, graph=False)
#            
#            date_min = date_min + timedelta(days=1)
#            date_max = date_min + timedelta(days= 2*delta_days + 1)
#            
#            if not results:
#                continue
#            window_mean_infections.append(results[country][0])
#            window_mean_deaths.append(results[country][1])
#            window_mean_mortality.append(results[country][2:4])
#            
#        
#        if window_mean_infections and window_mean_deaths and window_mean_mortality:
#            graphRatesEvolution(window_mean_infections, window_mean_deaths, 
#                                window_mean_mortality, date_min_MIN, delta_days)
            
    
