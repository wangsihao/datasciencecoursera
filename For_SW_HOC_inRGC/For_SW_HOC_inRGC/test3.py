import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import sys
sys.path.append('/Users/plesser/NEST/code/trunk/bld_fixes_mpi/install/lib/python2.7/site-packages/')

import nest

def noise_params(V_mean, V_std, dt=1.0, tau_m=10., C_m=250.):
    'Returns mean and std for noise generator for parameters provided; defaults for iaf_psc_alpha.'
    
    return C_m / tau_m * V_mean, math.sqrt(2/(tau_m*dt))*C_m*V_std




def simulate(mu, sigma, dt=1.0, t_max=5000.):

    
    resolution = 0.001
    
    nest.ResetKernel()
    nest.SetKernelStatus({'resolution': resolution})
    ng = nest.Create('noise_generator', params={'mean': mu, 'std': sigma, 'dt': dt})
    vm = nest.Create('voltmeter', params={'interval': resolution})
    nrns = nest.Create('iaf_psc_alpha')
    nest.Connect(ng, nrns)
    nest.Connect(vm, nrns)
    
    nest.Simulate(t_max)
    
    ev = nest.GetStatus(vm, keys=['events'])[0][0]
    t, s, v = ev['times'], ev['senders'], ev['V_m']

    
    return t, s, v
  
plt.figure(1)

dt = 0.1
mu, sigma = noise_params(0., 1., dt=dt)
#print "mu = {:.2f}, sigma = {:.2f}".format(mu, sigma)

#t, s, v = simulate(mu, sigma, dt=dt)
t, s, v = simulate(0., 100., dt=dt)
plt.subplot(211)
plt.plot(t,v,'--',color='red',label='dt = 0.1')

n1, bins1, patches1 = plt.hist(v, 50, normed=1, facecolor='green', alpha=0.75)

dt = 0.001
mu, sigma = noise_params(0., 1., dt=dt)
#print "mu = {:.2f}, sigma = {:.2f}".format(mu, sigma)

#t, s, v = simulate(mu, sigma, dt=dt)
t, s, v = simulate(0., 1000., dt=dt)
plt.subplot(212)
plt.plot(t,v,'--',color='blue',label='dt = 0.001')


n2, bins2, patches2 = plt.hist(v, 50, normed=1, facecolor='red', alpha=0.75)


plt.legend()
plt.xlabel('time t (ms)')
plt.ylabel('voltage V (mV)')
plt.xlim(0, 5000);



plt.figure(2)


plt.subplot(221)
plt.plot(bins, label='dt = 0.001')
plt.xlabel('Voltage')
plt.ylabel('Probability')
'''


plt.show()
