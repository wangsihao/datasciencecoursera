function [f] = generate_noise_conductances(nn, c, total_time, dt)
%function f = generate_noise_conductances(nn)
%generates conductances using f(sqrt(c)*n_1 + sqrt(1-c)*n_2)
%where f is the nonlinearly filtered noise stimuli 
%nn is 1 for excitatory and 2 for inhibitory stimuli
%c is the correlation coefficient

%total_time = 10*10^3; %ms
%dt = 0.01; %ms
t = [0:dt:total_time];

%noise
noise_mean = [30 -1200]; %exc and inh
noise_std = [500 780];

%generate two noise stimuli
tau_c = [22 33]; %correlation timescale for inputs in ms

ns = zeros(length(t),4);
%nc = zeros(length(t),2);

ex1 = exp(-dt/tau_c(nn));
sqex2 = sqrt(1-exp(-2*dt/tau_c(nn)));

for ii=1:length(t)-1         
    x1 = randn;
    x2 = randn;
    x3 = randn;
    x_c = randn;
     
    
    %OU correlated noise
    ns(ii+1,:) = ns(ii,:).* ex1 + ...
            noise_mean(nn).*(1-ex1) + ...
            noise_std(nn)*sqex2.*[x1 x2 x3 x_c]; 
    
    %the following two are the same as the above, which is later combined into stim
%     stim = nc
%     nc(ii+1,1) = nc(ii,1).* exp(-dt/tau_c(nn)) + ...
%             noise_mean(nn).*(1-exp(-dt/tau_c(nn))) + ...
%             noise_std(nn)*sqrt(1-exp(-2*dt/tau_c(nn))).*(sqrt(c)*x_c + sqrt(1-c)*x1); 
%         
%     nc(ii+1,2) = nc(ii,2).* exp(-dt/tau_c(nn)) + ...
%             noise_mean(nn).*(1-exp(-dt/tau_c(nn))) + ...
%             noise_std(nn)*sqrt(1-exp(-2*dt/tau_c(nn))).*(sqrt(c)*x_c + sqrt(1-c)*x2);    

end

%not temporally correlated
%ns = randn(length(t),3);

n_1 = ns(:,1);
n_2 = ns(:,2);
n_3 = ns(:,3);
n_c = ns(:,4);

%total noise input
stim1 = sqrt(c)*n_c + sqrt(1-c)*n_1;
stim2 = sqrt(c)*n_c + sqrt(1-c)*n_2;
stim3 = sqrt(c)*n_c + sqrt(1-c)*n_3;

stim = [stim1 stim2 stim3];

f = stim;
%M = [-0.00005 0.0001];
%N = [0.42 0.37]; %excitatory and inhibitory
%P = [-57 250];
%
%%nonlinear filter
%f = M(nn)*stim.^2 + N(nn)*stim + P(nn);

