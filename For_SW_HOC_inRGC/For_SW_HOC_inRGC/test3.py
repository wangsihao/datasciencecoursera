import math
import numpy as np
import scipy
import matplotlib.pyplot as plt


import sys
sys.path.append('/Users/plesser/NEST/code/trunk/bld_fixes_mpi/install/lib/python2.7/site-packages/')

import nest

def noise_params(V_mean, V_std, dt=1.0, tau_m=10., C_m=250.):
    'Returns mean and std for noise generator for parameters provided; defaults for iaf_psc_alpha.'
    
    return C_m / tau_m * V_mean, math.sqrt(2/(tau_m*dt))*C_m*V_std

def V_asymptotic(mu, sigma, dt=1.0, tau_m=10., C_m=250.):
    'Returns asymptotic mean and std of V_m'
    
    V_mean = mu * tau_m / C_m
    V_std = (sigma * tau_m / C_m) * np.sqrt(( 1 - math.exp(-dt/tau_m) ) / ( 1 + math.exp(-dt/tau_m) ))
                                    
    return V_mean, V_std
  
def V_mean(t, mu, tau_m=10., C_m=250.):
    'Returns predicted voltage for given times and parameters.'
    
    vm, _ = V_asymptotic(mu, sigma, tau_m=tau_m, C_m=C_m)
    return vm * ( 1 - np.exp( - t / tau_m ) )

def V_std(t, sigma, dt=1.0, tau_m=10., C_m=250.):
    'Returns predicted variance for given times and parameters.'
    
    _, vms = V_asymptotic(mu, sigma, dt=dt, tau_m=tau_m, C_m=C_m)
    return vms * np.sqrt(1 - np.exp(-2*t/tau_m))
  



def simulate(mu, sigma, dt=1.0, tau_m=10., C_m=250., N=1000, t_max=50.):

    
    resolution = 0.001

    ng = nest.Create('noise_generator', params={'mean': mu, 'std': sigma, 'dt': dt})
    vm = nest.Create('voltmeter', params={'interval': resolution})
    nrns = nest.Create('iaf_psc_alpha')
    nest.Connect(ng, nrns)
    nest.Connect(vm, nrns)
    
    nest.Simulate(t_max)
    
    ev = nest.GetStatus(vm, keys=['events'])[0][0]
    t, s, v = ev['times'], ev['senders'], ev['V_m']

    
    return t, s, v
  
'''

dt = 0.1
mu, sigma = noise_params(0., 1., dt=dt)
#print "mu = {:.2f}, sigma = {:.2f}".format(mu, sigma)

t, s, v = simulate(mu, sigma, dt=dt)
plt.plot(t,v,'--',color='red',label='dt = 0.1')

plt.legend()
plt.xlabel('time t (ms)')
plt.ylabel('voltage V (mV)')
plt.xlim(0, 50);
'''

dt = 0.01
mu, sigma = noise_params(0., 1., dt=dt)
#print "mu = {:.2f}, sigma = {:.2f}".format(mu, sigma)

t, s, v = simulate(mu, sigma, dt=dt)
plt.plot(t,v,'.',color='blue',label='dt = 0.01')

plt.legend()
plt.xlabel('time t (ms)')
plt.ylabel('voltage V (mV)')
plt.xlim(0, 50);


plt.show()
