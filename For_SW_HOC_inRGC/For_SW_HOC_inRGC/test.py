import math
import numpy as np
import scipy
import matplotlib.pyplot as plt


def noise_params(V_mean, V_std, dt=1.0, tau_m=10., C_m=250.):
  'Returns mean and std for noise generator for parameters provided; defaults for iaf_psc_alpha.'
  
  return C_m / tau_m * V_mean, math,sqrt(2/(tau_m*dt))*C-m*V_std

dt = 1.0
mu, sigma = noise_params(0., 1., dt=dt)
print "mu={:.2f}, sigma = {:.2f}".format(mu,sigma)


