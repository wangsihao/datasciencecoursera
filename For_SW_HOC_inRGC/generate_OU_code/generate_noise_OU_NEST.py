import nest
import math

def generate_noise_OU_NEST(tau, mu, sd, T):
    #all inputs need to be floats
    #T is the total simulation time
    #mu is the mean for noise gen
    #sd is the standard deviation for noise gen
    #I_e is the constant background current for the neuron
    #V_reset is the reset potential for the neuron
    #V_th is the threshold for the neuron
   
    white_noise_std = math.sqrt(2/tau)*sd

    noise = nest.Create('noise_generator',
                        params = {'mean': 0., 'std':white_noise_std})

    #noise = nest.Create('noise_generator',
    #                    params = {'mean': mu, 'std':sd})

    neuron = nest.Create('iaf_neuron',
                         params = {'I_e': 0.,
                                   'V_reset': 0.,
                                   'E_L': mu,
                                   'V_th': 10.**6,
                                   'tau_m': tau,
                                   'C_m': 1.0})
    
    nest.Connect(noise, neuron)

    #nest.Simulate(T)

    return(neuron)
