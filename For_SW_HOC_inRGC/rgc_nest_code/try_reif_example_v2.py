import matplotlib
matplotlib.use('TkAgg')


import matplotlib.pyplot as plt

import nest
import pylab


# Install module
nest.Install("reifmodule") 

nest.Models() #list of all available models within NEST

tmp=nest.GetDefaults('reif_cond_alpha_RK5')
print tmp

tmp=nest.GetDefaults('dc_generator')
print tmp


# Parameters for kernel
nest.SetKernelStatus({'resolution': 0.01})


#create a neuron of type reif_...
neurons = nest.Create('reif_cond_alpha_RK5', 2)

neuronsA = nest.Create('aeif_cond_alpha_RK5', 2)

#device used to record membrane voltage of neuron over time
multimeters = nest.Create('multimeter', 2,
                          params = {'withtime': True, 'record_from':['V_m'], 'interval': 0.01})
multimetersA = nest.Create('multimeter', 2,
                          params = {'withtime': True, 'record_from':['V_m'], 'interval': 0.01})



multimeters_w = nest.Create('multimeter', 2,
                            params = {'withtime': True, 'record_from':['w'], 'interval': 0.1})

#device that records spiking events produced by neurons
spikedetectors = nest.Create('spike_detector', 2,
                            params = {'withgid': True, 'withtime': True})


# WHAT PARAMETERS DOES THIS NEED??
inputnoises  = nest.Create('noise_generator', 2, [{'mean': 600., 'std': 200.}, {'mean': 300., 'std': 1000.}])

inputcurrent  = nest.Create('dc_generator', 2, [{'amplitude': 4000.}, {'amplitude': 0.}])

#connect devices to form network
nest.Connect(multimeters, neurons)
nest.Connect(multimeters_w, neurons)

# From aeif neurons, we will only measure voltage
nest.Connect(multimetersA, neuronsA)


nest.Connect(neurons, spikedetectors)

nest.Connect(inputcurrent, neurons, {'rule': 'one_to_one'})
nest.Connect(inputcurrent, neuronsA, {'rule': 'one_to_one'})

'''The order in which the arguments to Connect are specified
reflects the flow of events: if the neuron spikes,it send an event to the
spike detector. Conversely, the multimeter periodically sends requests
to the neuron to ask for its membrane potential at that point in time.'''




#run simulation in ms
maxT  = 400.0
nest.Simulate(maxT)


#blergh = nest.GetStatus(neurons)
#print blergh


#extracting and plotting data from devices

'''
# Use this to trhy and test X11 forwarding...
Vms  = dmm['events']['V_m']
ts   = dmm['events']['times']

plt.figure(1)
plt.plot(ts, Vms,'b.-')

'''

#### rEIF neurons: Voltage
dmm = nest.GetStatus(multimeters)[0]

Vms  = dmm['events']['V_m']
ts   = dmm['events']['times']
sdrs = dmm['events']['senders']

n1ind = pylab.find(sdrs==1)
n2ind = pylab.find(sdrs==2)

Vms1 = Vms[n1ind]
ts1  = ts[n1ind]

Vms2 = Vms[n2ind]
ts2  = ts[n2ind]

plt.figure(1)
plt.plot(ts1, Vms1,'b-')
plt.plot(ts2, Vms2,'r-')



#### Now for AEIF
dmmA = nest.GetStatus(multimetersA)[0]

VmsA  = dmmA['events']['V_m']
tsA   = dmmA['events']['times']
sdrsA = dmmA['events']['senders']

n1ind = pylab.find(sdrsA==3)
n2ind = pylab.find(sdrsA==4)

VmsA1 = VmsA[n1ind]
tsA1  = tsA[n1ind]

VmsA2 = VmsA[n2ind]
tsA2  = tsA[n2ind]

#print dmmA

plt.plot(tsA1, VmsA1,'g-')
plt.plot(tsA2, VmsA2,'c-')


#### Spike train from rEIF cells
dSD = nest.GetStatus(spikedetectors,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
plt.figure(2)
plt.plot(ts, evs, ".")
plt.axis([0,maxT,0.5,2.5])


#### "w" variable: should track time since last spike
dmw = nest.GetStatus(multimeters_w)[0]

ws  = dmw['events']['w']
tsb = dmw['events']['times']
sdrsb = dmw['events']['senders']
n1ind = pylab.find(sdrsb==1)
n2ind = pylab.find(sdrsb==2)

ws1 = ws[n1ind]
tsb1  = tsb[n1ind]

ws2 = ws[n2ind]
tsb2  = tsb[n2ind]

plt.figure(3)
plt.plot(tsb1, ws1)
plt.plot(tsb2, ws2,'r-')


#
plt.show()
