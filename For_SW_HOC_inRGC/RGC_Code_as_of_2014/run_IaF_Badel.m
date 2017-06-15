function [prob_d,varargout] = run_IaF_Badel(StimParam)
% Old arguments
%(stim_std,flag,t_refresh,stixel_size,varargin)
%
%
% type_stim = 0 for noise only
%            1 for full field
%            2 for stixels
%
%if type_stim = 1, then need
% flag = 1; %1 for binary, 0 for Gaussian inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = StimParam.N;   %number of neurons
%N = 3; 

generate_feature_space_01; %generate feature space of 0 and 1

total_time = 10^5; %ms; up from 10^4 in previous sims 
dt = StimParam.dt;   %ms


%Matlab implementaion of Badel et al fit            
[Spikes, counter_bins, StimParam] = simulate_spikes_Badel(total_time,dt,StimParam);
%counter_bins = total number of extra spikes (more than one spike in a bin)
if (counter_bins > 0)
    sprintf('counter_bins = %d, numspikes = %d',counter_bins,sum(sum(Spikes)))
end

Spikes_one_bin = Spikes;
%put Spikes into -1 and 1, faster computation to get
%probability distribution
Spikes_one_bin(Spikes_one_bin == 0) = -1;
state(state == 0)= -1;

%compute prob_d -- probability of each state occuring in the data
holder2 = zeros(2^N,1);
for i=1:2^N
    holder = Spikes_one_bin * state(i,:)' ;  %for each sample in the distrib 
    %dot with the state in question.
    holder(holder ~= N) = 0;  %if the dot was not N, the sample and state do not match
    holder(holder == N) = 1;  %else count it!
    holder2(i) = sum(holder);
end
prob_d = holder2/sum(holder2);  %probability of each state occuring

nout = max(nargout,1)-1;
if (nout > 0)
    % Pass back StimParam, in case changes have been made here
    varargout(1) = {StimParam};
end
    

%SN = -sum(prob_d.*log2(prob_d));

%print out probabilities
%prob_d
