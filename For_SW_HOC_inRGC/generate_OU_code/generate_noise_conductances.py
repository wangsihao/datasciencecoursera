import numpy as np
import math as m

def generate_noise_conductances(nn, c, total_time, dt):
    #function f = generate_noise_conductances(nn)
    #generates conductances using f(sqrt(c)*n_1 + sqrt(1-c)*n_2)
    #where f is the nonlinearly filtered noise stimuli 
    #nn is 1 for excitatory and 2 for inhibitory stimuli
    #c is the correlation coefficient
    
    #total_time = 10*10**3 #ms
    #dt = 0.01 #ms
    t = np.arange(0, total_time + dt, dt)

    # Note: explicitly cast constants as floating point numbers if applicable! 
    # This prevents division from returning integer division 

    #noise
    noise_mean = np.array([30., -1200.]) #exc and inh
    noise_std = np.array([500., 780.])

    #generate two noise stimuli
    tau_c = np.array([22., 33.]) #correlation timescale for inputs in ms

    #ns = np.zeros(shape = (len(t), 4), dtype = int)
    ns = np.zeros(shape = (len(t), 4))
    #nc = np.zeros(shape = (len(t), 2))

    ex1 = m.exp(-dt/tau_c[nn-1])
    sqex2 = m.sqrt(1-m.exp(-2*dt/tau_c[nn-1]))

    for ii in range(1,len(t)):
        x1 = np.random.normal()
        x2 = np.random.normal()
        x3 = np.random.normal()
        x_c = np.random.normal()

        #OU correlated noise
        ns[ii,] = (ns[ii-1,] * ex1 +
                    noise_mean[nn-1] * (1 - ex1) +
                    np.dot(noise_std[nn-1], (sqex2 * np.array([x1, x2, x3, x_c]))))

        #the following two are the same as the above, which is later combined into stim
        '''stim = nc
        
        nc[ii,0] = (nc[ii-1,0] * m.exp(-dt/tau_c[nn-1]) +
                    noise_mean[nn-1] * (1 - m.exp(-dt/tau_c[nn-1])) +
                    np.dot(noise_std[nn-1],
                    (m.sqrt(1 - m.exp(-2*dt/tau_c[nn-1])) *
                    (np.dot(m.sqrt(c), x_c) + np.dot(m.sqrt(1-c), x1)))))

        nc[ii,1] = (nc[ii-1,1] * m.exp(-dt/tau_c[nn-1]) +
                    noise_mean[nn-1] * (1 - m.exp(-dt/tau_c[nn-1])) +
                    np.dot(noise_std[nn-1],
                    (m.sqrt(1 - m.exp(-2*dt/tau_c[nn-1])) *
                    (np.dot(m.sqrt(c), x_c) + np.dot(m.sqrt(1-c), x2)))))'''
        

    #not temporarily correlated
    #ns = np.random.normal(size = (len(t), 3))

    n_1 = ns[:,0]
    n_2 = ns[:,1]
    n_3 = ns[:,2]
    n_c = ns[:,3]

    #total noise input
    stim1 = np.dot(m.sqrt(c), n_c) + np.dot(m.sqrt(1-c), n_1)
    stim2 = np.dot(m.sqrt(c), n_c) + np.dot(m.sqrt(1-c), n_2)
    stim3 = np.dot(m.sqrt(c), n_c) + np.dot(m.sqrt(1-c), n_3)

    stim = np.array([stim1, stim2, stim3])

    return stim
    #return stim[0,:]
    
    #M = np.array([-0.00005, 0.0001])
    #N = np.array([0.42, 0.37]) #excitatory and inhibitory
    #P = np.array([-57, 250])

    #nonlinear filter
    #f = np.dot(M[nn-1], stim**2) + np.dot(N[nn-1], stim) + P[nn-1]

