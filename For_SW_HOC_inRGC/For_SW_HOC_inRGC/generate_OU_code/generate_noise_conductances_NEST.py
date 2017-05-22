import nest

# Sihao: This function made a neuron whose voltage should 
# have the statistics of an OU process
#
# However, I don't think I ever got it to ACTUALLY function 
#  as an input to another neuron!

def generate_noise_conductances(T, mu, sd, I_e, V_reset, V_th):
    #all inputs need to be floats
    #T is the total simulation time
    #mu is the mean for noise gen
    #sd is the standard deviation for noise gen
    #I_e is the constant background current for the neuron
    #V_reset is the reset potential for the neuron
    #V_th is the threshold for the neuron
    
    white_noise_std = sqrt(2/tau)*sd

    noise = nest.Create('noise_generator',
                        params = {'mean': mu, 'std':white_noise_std})

    neuron = nest.Create('iaf_neuron',
                         params = {'I_e': I_e,
                                   'V_reset': V_reset,
                                   'V_th': V_th})
    
    nest.Connect(noise, neuron)
    nest.Simulate(T)

    return(neuron)
