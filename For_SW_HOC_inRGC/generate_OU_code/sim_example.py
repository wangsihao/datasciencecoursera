import nest
import numpy as np
import math as m
import pylab
from generate_noise_conductances_NEST import *
from generate_noise_OU_NEST import *

#generating noise conductances
T = 100.0

tau    = 22.
mu     = 0.
sd     = 1000.
noise_neuron = generate_noise_OU_NEST(tau, mu, sd, T)


#create a neuron of type iaf_neuron
neuron = nest.Create('iaf_neuron',
                     params = {'I_e': 376.0})


noisemeter = nest.Create('multimeter',
                         params = {'withtime': True, 'record_from':['V_m']})


#device used to record membrane voltage of neuron over time
multimeter = nest.Create('multimeter',
                         params = {'withtime': True, 'record_from':['V_m']})


#device that records spiking events produced by neurons
spikedetector = nest.Create('spike_detector',
                            params = {'withgid': True, 'withtime': True})


#connect devices to form network
nest.Connect(noisemeter, noise_neuron)
nest.Connect(noise_neuron, neuron)
nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)

'''The order in which the arguments to Connect are specified
reflects the flow of events: if the neuron spikes,it send an event to the
spike detector. Conversely, the multimeter periodically sends requests
to the neuron to ask for its membrane potential at that point in time.'''


#run simulation in ms
nest.Simulate(T)


#extracting and plotting data from devices
dmm = nest.GetStatus(multimeter)[0]
Vms = dmm['events']['V_m']
ts = dmm['events']['times']
pylab.figure(1)
pylab.plot(ts, Vms)

dmm = nest.GetStatus(noisemeter)[0]
Vms = dmm['events']['V_m']
ts = dmm['events']['times']
pylab.figure(3)
pylab.plot(ts, Vms)

EX = Vms.mean()
CM = m.sqrt(np.var(Vms))

print 'mean = ', EX, 'std = ', CM

autocorr = np.correlate(Vms,Vms,'full')
pylab.figure(4)
taxis = np.arange(len(autocorr))-(len(autocorr))/2
pylab.plot(taxis,autocorr)
pylab.xlim([-200.,200.])
pylab.grid()
print len(autocorr)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(2)
pylab.plot(ts, evs, ".")



pylab.show()
