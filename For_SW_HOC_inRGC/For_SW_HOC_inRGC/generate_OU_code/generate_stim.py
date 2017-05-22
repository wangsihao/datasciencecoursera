import numpy as np
import math as m

def generate_stim(t_list, flag, Smean, Sstd, t_refresh):
    #generates stimulus given a flag: 
    #if flag=0 then Gaussian
    #if flag=1 then binary

    if flag == 1:
        Smax = Smean + Sstd
        Smin = Smean - Sstd

    '''takes vector argument t_list in ** SEC **, produces list the same length of
    binary values that switch randomly between Smax and Smin every t_refresh'''

    t_length = len(t_list)
    F = np.zeros(shape = (np.size(t_list)))
    dt = t_list[1] - t_list[0]
    Tmax = t_list[-1] - t_list[0] + dt

    #do a crude loop through t_list, setting stimulus values ...
    stim_index = 1
    stim_start = t_list[0]

    #initial value of stim
    if flag == 1:
        r = np.round(np.random.randint(0, 2))   # Choose 0 or 1, equal prob
        stim = np.dot(r, Smax) + np.dot(1-r, Smin) #takes either value w/ equal prob
    else:
        stim = Smean + np.dot(Sstd, np.random.normal()) #takes either value w/ equal prob
        if stim < -1:
            #resample
            while stim < -1:
                stim  = Smean + np.dot(Sstd, np.random.normal())
            #stim = -1
                
    stim_time = t_list[0]

    for jj in range(1, t_length + 1):
       
        if (t_list[jj-1] - stim_time) > t_refresh:
            if flag == 1:
                r = np.round(np.random.randint(0, 2))   # Choose 0 or 1, equal prob
                stim = np.dot(r, Smax) + np.dot(1-r, Smin) #takes either value w/ equal prob
            else:
                stim = Smean + np.dot(Sstd, np.random.normal()) #takes either value w/ equal prob

            # Keep correct refresh intervals
            stim_time =    stim_time + t_refresh #t_list[jj-1]

        F[jj-1] = stim


    return F
