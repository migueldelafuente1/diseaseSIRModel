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

t_0 = datetime(2020, 2, 15)

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
spain_detected = [
#                   (8, 2, 0),
#                   (12, 25, 0),
#                   (14, 58, 0),
#                   (16, 120, 0),
#                  (17, 165, 1),
#                  (19, 282, 3),
#                  (21, 525, 10),
#                  (22, 674, 17),
#                  (23, 1231, 30),
#                  (24, 1695, 36),
#                  (25, 2277, 55),
#                  (26, 3146, 86),
#                  (27, 5232, 133),
                  (28, 6391, 196),
                  (29, 7988, 294),
                  (30, 9942, 342),
                  (31, 11826, 533),
                  (32, 14769, 638),
                  (33, 18077, 831)]

italy_detected = [
#                  (0, 3, 0),
#                  (6, 21, 1),
#                  (9, 229, 7),
#                  (11, 470, 12),
#                  (13, 889, 21),
#                  (15, 1701, 41),
#                  (17, 2502, 79),
#                  (18, 3089, 107),
#                  (19, 3858, 148),
#                  (20, 4636, 197),
#                  (21, 5883, 233),
#                  (22, 7375, 366),
#                  (23, 9172, 463),
#                  (24, 10149, 631),
#                  (25, 12462, 827),
                  (26, 15113, 1016),
                  (27, 17660, 1266),
                  (28, 21157, 1441),
                  (29, 24747, 1809),
                  (30, 27980, 2158),
                  (31, 31506, 2503),
                  (32, 35713, 2978),
                  (33, 41035, 3405)]


def graph(data_list, A_i, B_i, A_d, B_d, t_max, name):
    
    t, infected, dead = zip(*data_list)
    fig, ax = plt.subplots()
    ax.plot(t, infected, 'r^-',label=f'{name} infected')
    ax.plot(t, dead, 'kv-',label=f'{name} dead')
    
    x_i, y_i = zip(*[(t, B_i* (10**(A_i*t)))
                     for t in np.arange(10, t_max, 0.01)])
    x_d, y_d = zip(*[(t, B_d* (10**(A_d*t))) 
                     for t in np.arange(10, t_max, 0.01)])
    ax.plot(x_i, y_i, 'r--',label=f'{name} infected')
    ax.plot(x_d, y_d, 'k--',label=f'{name} dead')
    ax.legend()
    ax.grid()
#    ax.semilogy()
    ax.set_xlabel('time (days)')
    
    

if __name__ == '__main__':
    
    t_max = 34
    dicts_ = {'spain': spain_detected, 
              'italy':italy_detected}
    for name, detected_list in dicts_.items():
        print(name.upper()+"\n-----------------------------------------")
        vals = getContagiousRate(detected_list)
        i_rate, i_y_abs, d_rate, d_y_abs, mort, mort_std = vals
        print(f"Infect Rate: {i_rate:8.4f}   abscI: {i_y_abs:8.4f}")
        print(f"Death  Rate: {d_rate:8.4f}   abscD: {d_y_abs:8.4f}")
        print(f"Mortality:  {mort:8.4f} +/- {mort_std:6.4f}")
        
        _infected = round(i_y_abs * (10**(i_rate*t_max)))
        _deaths   = round(d_y_abs * (10**(d_rate*t_max)))
        print(f"\t{(t_0+timedelta(days=t_max)).date()} : infected= {_infected}")
        print(f"\t{(t_0+timedelta(days=t_max)).date()} : deaad= {_deaths}\n")
        
        graph(detected_list, i_rate, i_y_abs, d_rate, d_y_abs, t_max, name)
