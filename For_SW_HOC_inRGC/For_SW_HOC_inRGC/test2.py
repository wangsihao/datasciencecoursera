import nest
import numpy as np
import pylab as pl

taum = 10.
C_m = 250.
# array of distances between tau_m and tau_ex
epsilon_array = np.hstack(([0.],10.**(np.arange(-6.,1.,1.))))[::-1]
dt = 0.1
fig = pl.figure(1)
NUM_COLORS = len(epsilon_array)
cmap = pl.get_cmap('gist_ncar')
maxVs = []

for i,epsilon in enumerate(epsilon_array):
    nest.ResetKernel() # reset simulation kernel 
    nest.SetKernelStatus({'resolution':dt})

    # Current based alpha neuron 
    neuron = nest.Create('iaf_psc_alpha') 
    nest.SetStatus(neuron,{'C_m':C_m,'tau_m':taum,'t_ref':0.,'V_reset':-70.,'V_th':1e32,
                           'tau_syn_ex':taum+epsilon,'tau_syn_in':taum+epsilon,'I_e':0.})
   
    # create a spike generator
    spikegenerator_ex=nest.Create('spike_generator')
    nest.SetStatus(spikegenerator_ex,{'spike_times': [50.]})
    
    # create a voltmeter
    vm = nest.Create('voltmeter',params={'interval':dt})

    ## connect spike generator and voltmeter to the neuron

    nest.Connect(spikegenerator_ex, neuron,'all_to_all',{'weight':100.})
    nest.Connect(vm, neuron)

    # run simulation for 200ms
    nest.Simulate(200.) 

    # read out recording time and voltage from voltmeter
    times=nest.GetStatus(vm)[0]['events']['times']
    voltage=nest.GetStatus(vm)[0]['events']['V_m']
    
    # store maximum value of voltage trace in array
    maxVs.append(np.max(voltage))

    # plot voltage trace
    if epsilon == 0.:
        pl.plot(times,voltage,'--',color='black',label='singular')
    else:
        pl.plot(times,voltage,color = cmap(1.*i/NUM_COLORS),label=str(epsilon))

pl.legend()
pl.xlabel('time t (ms)')
pl.ylabel('voltage V (mV)')


fig = pl.figure(2)
pl.semilogx(epsilon_array,maxVs,color='red',label='maxV')
#show singular solution as horizontal line
pl.semilogx(epsilon_array,np.ones(len(epsilon_array))*maxVs[-1],color='black',label='singular')
pl.xlabel('epsilon')
pl.ylabel('max(voltage V) (mV)')
pl.legend()


pl.show()
