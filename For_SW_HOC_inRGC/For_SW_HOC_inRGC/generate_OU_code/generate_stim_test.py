import numpy as np
import math as m

from generate_stim import *
from generate_noise_conductances import *

# Call generate_noise_conductances

total_time = 10000.
dt = 1.
c  = 0.2
nn = 1


fnc = generate_noise_conductances(nn,c,total_time,dt)

print 'Info about conductance array'
print type(fnc)
print fnc.shape

# Ok: we are ok with shape. Returns a 3 x (total_time/dt)+1 array
# print 'fnc= ', fnc


## First moment (mean)
EX = fnc.mean(axis=1)


### Covariance matrix
CM = np.cov([fnc[0,:],fnc[1,:],fnc[2,:]])

print 'EX = ', EX
print 'CM = ', CM

# Element-wise square root
CMSQ =  np.sqrt(CM)
print 'CMSQ = ', CMSQ

### We are ok with values 






# Call generate_stim

t_list = np.arange(0, total_time+dt, dt)
flag   = 1
Smean  = 1.
Sstd   = 1.5
t_refresh = 2.

F  = generate_stim(t_list,flag,Smean,Sstd,t_refresh)

print 'Info about stim array'
print np.shape(F)
print type(F)

print 'F= ', F

print 'mean = ', F.mean()
print 'var = ', F.var()
