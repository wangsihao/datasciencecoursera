import nest
import pylab

#nest.Models() #list of all available models within NEST

#create a neuron of type aeif_...
neuron1 = nest.Create('aeif_cond_alpha_RK5', 1)
neuron2 = nest.Create('aeif_cond_alpha_RK5', 1)

#device used to record membrane voltage of neuron over time
multimeter1 = nest.Create('multimeter', 1,
                         params = {'withtime': True, 'record_from':['V_m']})
multimeter2 = nest.Create('multimeter', 1,
                         params = {'withtime': True, 'record_from':['V_m']})

#device that records spiking events produced by neurons
spikedetector1 = nest.Create('spike_detector', 1,
                            params = {'withgid': True, 'withtime': True})
spikedetector2 = nest.Create('spike_detector', 1,
                            params = {'withgid': True, 'withtime': True})


# WHAT PARAMETERS DOES THIS NEED??
inputnoise1  = nest.Create('noise_generator', 1, [{'mean': 600., 'std': 200.}])
inputnoise2  = nest.Create('noise_generator', 1, [{'mean': 300., 'std': 1000.}])

#connect devices to form network
nest.Connect(multimeter1, neuron1)
nest.Connect(multimeter2, neuron2)
nest.Connect(neuron1, spikedetector1)
nest.Connect(neuron2, spikedetector2)


nest.Connect(inputnoise1, neuron1)
nest.Connect(inputnoise2, neuron2)
'''The order in which the arguments to Connect are specified
reflects the flow of events: if the neuron spikes,it send an event to the
spike detector. Conversely, the multimeter periodically sends requests
to the neuron to ask for its membrane potential at that point in time.'''



#run simulation in ms
nest.Simulate(1000.0)


#extracting and plotting data from devices
dmm = nest.GetStatus(multimeter1)[0]

dmm1 = nest.GetStatus(multimeter2)[0]

print dmm
print dmm1

Vms = dmm['events']['V_m']
ts = dmm['events']['times']
pylab.figure(1)
pylab.plot(ts, Vms)

Vms1 = dmm1['events']['V_m']
ts1 = dmm1['events']['times']
pylab.plot(ts1, Vms1,'r-')

dSD = nest.GetStatus(spikedetector1,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(2)
pylab.plot(ts, evs, ".")

dSD1 = nest.GetStatus(spikedetector2,keys='events')[0]
evs1 = dSD1["senders"]
ts1 = dSD1["times"]
pylab.plot(ts1, evs1, "r.")
pylab.axis([0,1000,0.5,2.5])

pylab.show()
