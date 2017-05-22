function [Spikes, counter_bins, varargout] = simulate_spikes_Badel_Recur_Rect(total_time,dt,StimParam)

% uses the integrate-and-fire model with the Badel fit
% to generate probabilities distributions
%

% JG - April 2010 mod Aug 2010
% AKB - mod June 2011 to take general stim structure
%    - mod June 2011 to have corrected Vreset, V_th, treatment of t_ref, etc..
%    - mod June 23 2011 to have RECURRENCE and STEREOTYPED SPIKE WAVEFORM

N = StimParam.N;
%N = 3; %number of neurons

% Get stimulus parameters from attached structure
type_stim = StimParam.type_stim;
marg_flag = StimParam.marg_flag;
t_refresh = StimParam.t_refresh;
stim_std = StimParam.stim_std;

% Parameters for 
tau_abs = 3;    % absolute refractory period base 4-6 (2-3 ms)
tau = 5;        % membrane time constant
%dt = tau/100;                       %integration timestep

if (dt > tau/10)
    disp('dt < tau/10; time step may be too small');
end
t_list = [0:dt:total_time];
t_list_filter = [0:dt:200];          % For filter

T_binning =     5;% 3*tau_abs;              %width of a bin
Tstop = ceil(total_time/T_binning); %number of steps (bins)
Spikes = 0*ones(Tstop,N);

%(each column has num_steps entries for each neuron 1,...,N) 
%1's to be added later if neuron nn (1 to N) spikes at time tt (1:duration)
NumSpikes = zeros(1,N);      %keeps a record of the number of spikes for each neuron
maxdimSp = 1000;
Spike_Times = zeros(N,maxdimSp);

% 7/13/11 - different iterations on model
%
model_number = 2.5;
if (model_number == 2.0)
    % Model 2.0
    %  Vth = -30
    %  Spike Shape from dynamic clamp 
    %  Hold for 2 ms
    varVth_allowed = 0;  % do we allow a variable threshold?
    V_th   = -30;
    HoldT  = 2;
    SSname = 'SpikeShape.mat';
elseif (model_number == 2.5)
    % Model 2.5
    %  Vth VARIABLE
    %  Spike Shape from dynamic clamp
    %  Hold for 0.95 ms
    varVth_allowed = 1;  % do we allow a variable threshold?
    V_th   = -30;
    HoldT  = 0.95;
    SSname = 'SpikeShape.mat';
elseif (model_number == 3.0)
    %
    % Model 3.0
    %  Vth = -30
    %  Spike Shape from CURRENT clamp
    %  Hold for 2 ms
    varVth_allowed = 0;  % do we allow a variable threshold?
    V_th   = -30;
    HoldT  = 2;
    SSname = 'SpikeShape_CurrentClamp.mat';
else
    varVth_allowed = 0;  % do we allow a variable threshold?
    V_th   = -30;
    HoldT  = 2;
    SSname = 'SpikeShape_CurrentClamp.mat';
end
%

%Neuron properties
Vreset = -55; %mV

% For now, threshold is default value
% (If varVth_allowed==0, then should remain here for entire sim)
V_th_temp = V_th;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create stereotyped waveform
% contains Spike waveform for 0:.1:SSTotalT ms.
load(SSname);
SSTotalT      = (length(SpikeShape)-1)*0.1;
sptime        = 0:dt:SSTotalT;
SpikeShapeNew = interp1(0:.1:SSTotalT,SpikeShape,sptime,'spline');

% Extract relevant part: find first point over -30 mV, 
% and 
%   **** HoldT ms ***** after
%
% This number chosen to eliminate hyperpolarization 
%      that might be an artifact...AKB - 7/1/11
% All parameters involved in this SHOULD be set here, so no change
%      needs to be made to integration loop below.
temp = find(SpikeShapeNew > V_th);
temp2 = find(sptime > sptime(temp(1))+HoldT);

SpikeShapeNew = SpikeShapeNew(temp(1):temp2(1));
sptime = sptime(temp(1):temp2(1))-sptime(temp(1));

% will maintain whether or not neuron is spiking
in_spike     = zeros(N,1);
max_in_spike = numel(sptime); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


counter_bins = 0;

V = Vreset*ones(N,1);  % a voltage vector for all the N neurons     
Vold = V;

% Set parameters for reversal potentials
% (*** actually, didn't I get these straight from Fred?
%     Took these from fitting of RGC data using procedure of Badel, et al.
%     ***)
VE = 0;
VI = -80;   % -80 mV

%correlation coefficients for exc and inh conductances
%used in the generation of noise conductances
c_exc = 0.30;
c_inh = 0.15;

%mean conductances i.e. tonic drive
G_exc = 30*0.1611; %nS
G_inh = 40*0.1634; %nS

% Parameters for recurrent gap junction
% 1.1 nS: as specified by Fred 1.1

if (isfield(StimParam,'gapg'))
    gapg = (10/9)*StimParam.gapg;
else
    gapg  = 10/9;   % nS
end

if (N==3)
    % Assume connectivity is symmetric and global
    gapJ  = gapg*[2 -1 -1;
             -1 2 -1;
             -1 -1 2];
end

% Must first generate currents by putting an OU process through 
% Fred's filter
I_exc_noise = generate_noise_conductances(1,c_exc,total_time,dt)'; 
I_inh_noise = generate_noise_conductances(2,c_inh,total_time,dt)'; 

if (type_stim == 0)
    disp('noise only')
    %nonlinear filter
    I_exc = 0;
    I_inh = 0;
    
elseif (type_stim == 1)
    disp('full field stim + noise')
    
    Sstd = StimParam.stim_std;
    Smean = 0;  %NOTE!  Filter set so that the stimuli we use -- as for those Fred used to fit the filter -- should have mean 0
    
    %create a realiz of stim vs. time
    stim = generate_stim(t_list,marg_flag,Smean,Sstd,t_refresh);

%     %parameters of filter f (1 for exc and 2 for inh)
%     A = [-8e4 -1.8e5]/10^3;
%     n = [3.6 3];
%     tau = [12 16]; %ms
%     T = [105 120]; %ms
%     
%     filter_exc = A(1)*(t_list_filter/tau(1)).^n(1).*exp(-(t_list_filter)/tau(1)).*sin(2*pi*t_list_filter/T(1));
%     filter_inh = A(2)*(t_list_filter/tau(2)).^n(2).*exp(-(t_list_filter)/tau(2)).*sin(2*pi*t_list_filter/T(2));
%     

    % P filter
    AP = 4*[-8e4 -8.5e4]/10^3;
    nP = [2 2];
    tauP = [12 13.2]; %ms
    
    TS = [120 132]; %ms
    magC = [0.8 0.8];
    TC = [100 110]; %ms
    filter_exc = AP(1)*(t_list_filter/tauP(1)).^nP(1).*exp(-(t_list_filter)/tauP(1)) ...
        .*(sin(2*pi*t_list_filter/TS(1)) + magC(1) * cos(2*pi*t_list_filter/TC(1)));
    filter_inh = AP(2)*(t_list_filter/tauP(2)).^nP(2).*exp(-(t_list_filter)/tauP(2)) ...
        .*(sin(2*pi*t_list_filter/TS(2)) + magC(2) * cos(2*pi*t_list_filter/TC(2)));
    
    I_exc = conv(filter_exc,stim)*dt;
    I_exc = repmat(I_exc(1:length(t_list)),N,1);
    
    I_inh = conv(filter_inh,stim)*dt;
    I_inh = repmat(I_inh(1:length(t_list)),N,1);
    
elseif (type_stim == 2)
   disp('stixel stim + noise')
   % RANDOMNESS PARAMETER - if 1, then 
   % 1) rotate receptive field randomly and 
   % 2) shift central location of grid
   %
   % IF 2: parameters are passed as variable length argument
    
    %in MICROMETERS
    Stixel_Size = StimParam.stixel_size;
    
    Smean = 0;  %NOTE!  Filter set so that the stimuli we use -- as for those Fred used to fit the filter -- should have mean 0
    Sstd = StimParam.stim_std;
    
    % create a realiz of stim vs. time
    % Here, return new value of StimParam:
    %      this may have changed during call! 
    [stim,StimParam] = generate_stim_stixel(t_list,marg_flag,Smean,Sstd,StimParam);
    
%     %parameters of filter f (1 for exc and 2 for inh)
%     A = [-8e4 -1.8e5]/10^3;
%     n = [3.6 3];
%     tau = [12 16]; %ms
%     T = [105 120]; %ms
%     
%     % Added shorter filter length - 10/21
%     filter_exc = A(1)*(t_list_filter/tau(1)).^n(1).*exp(-(t_list_filter)/tau(1)).*sin(2*pi*t_list_filter/T(1));
%     filter_inh = A(2)*(t_list_filter/tau(2)).^n(2).*exp(-(t_list_filter)/tau(2)).*sin(2*pi*t_list_filter/T(2));
%     
    % P filter
    AP = 4*[-8e4 -8.5e4]/10^3;
    nP = [2 2];
    tauP = [12 13.2]; %ms
    
    TS = [120 132]; %ms
    magC = [0.8 0.8];
    TC = [100 110]; %ms
    filter_exc = AP(1)*(t_list_filter/tauP(1)).^nP(1).*exp(-(t_list_filter)/tauP(1)) ...
        .*(sin(2*pi*t_list_filter/TS(1)) + magC(1) * cos(2*pi*t_list_filter/TC(1)));
    filter_inh = AP(2)*(t_list_filter/tauP(2)).^nP(2).*exp(-(t_list_filter)/tauP(2)) ...
        .*(sin(2*pi*t_list_filter/TS(2)) + magC(2) * cos(2*pi*t_list_filter/TC(2)));
    
    
    I_exc = [];I_inh = [];
    for k=1:N
        I_exc_temp = conv(filter_exc,stim(k,:))*dt;
        I_inh_temp = conv(filter_inh,stim(k,:))*dt;
        I_exc = [I_exc; I_exc_temp(1:length(t_list))]; 
        I_inh = [I_inh; I_inh_temp(1:length(t_list))];
    end

end

I_exc = I_exc + I_exc_noise;
I_inh = I_inh + I_inh_noise;

%parameters of filter g (1 for exc and 2 for inh)
M_param = [-0.00005 0.0001];
N_param = [0.42 0.37]; %excitatory and inhibitory
P_param = [-57 250];

I_exc = M_param(1)*I_exc.^2 + N_param(1)*I_exc + P_param(1);  
I_inh = M_param(2)*I_inh.^2 + N_param(2)*I_inh + P_param(2); 
   
%convert currents to conductances: divide by the driving force
g_exc = I_exc/(-60);
g_inh = I_inh/60;

%offset by tonic drive
g_exc = g_exc + G_exc;
g_inh = g_inh + G_inh;

% Save g_exc, if applicable
if (StimParam.save_gexc_flag && StimParam.random_flag==2)  
    ID_number = StimParam.ID_num;
    if (marg_flag ==1)
        outfname = sprintf('gexc_bin_stixel=%d_ID=%d.mat',Stixel_Size,ID_number);
        save(outfname,'g_exc','StimParam');
    else
        outfname = sprintf('gexc_gauss_stixel=%d_ID=%d.mat',Stixel_Size,ID_number);
        save(outfname,'g_exc','StimParam');
    end
end

%
% Coefficients 
%   (AKB 6/11: NOTE: UNCHANGED SINCE JUNE 2010)
P1 = [0.3719    0.5412   13.2726];
P2 = [-59.4858    5.8966    8.3076  233.1114];
P3 = [20.048672101308   19.055970786380   3.628022920530   2430.376633691003];
P4 = [-44.3323   25.1812    4.7653];

C = 28; %pF

% if (StimParam.ID_num == 1)
%     % Keep some information
%     Vkeep = [];
%     tkeep = [];
%     gexckeep = [];
%     ginhkeep = [];
% end

time_to_quit = total_time+1;

%integrate the membrane potential through LIF and OU for the input
for tt=2:total_time/dt
    
    current_time = tt*dt;   
       
    for nn=1:N %number of neurons
        
        % Figure out time since last spike
        if (NumSpikes(nn) > 0)
            ts = current_time - Spike_Times(nn,NumSpikes(nn));
        else
            ts = current_time;
        end
        f1 = get_function(P1,ts);
        f2 = get_function(P2,ts);
        f3 = get_function(P3,ts);
        f4 = get_function(P4,ts);

        F = f1*(f2 - Vold(nn) + f3*exp((Vold(nn)-f4)/f3));
        
        %forward Euler method integrates membrane potential for the given
        %neuron, unless we are "inside a spike" and using a stereotyped spike waveform
        if (in_spike(nn) == 0 )
            I_ion = g_exc(nn,tt-1)*(V(nn)-VE) + g_inh(nn,tt-1)*(V(nn)-VI);
            
            
            I_gap = gapJ(nn,:)*Vold;    % in same form as I_ion
                                  % gap junction current is linear in V
            
            V(nn) = Vold(nn) + dt*(-I_ion/C - I_gap/C + F);
        
            
           % Here, we try a variable threshold
           % Add an additional 2 mV on, because we should not assume we
           % will go up unless F > 0...otherwise, we may be incorrectly
           % triggering spikes
           if (varVth_allowed)
               if (ts < 6)
                   V_th_temp = 50*exp(-ts/4)-38;
               else
                   V_th_temp = V_th;
               end
           end
           
           % have we exceeded threshold?
           % and are we allowed to detect a spike?
           
           if (V(nn) > V_th_temp && (NumSpikes(nn)==0 || ts > tau_abs) )
               % V
                NumSpikes(nn) = NumSpikes(nn) + 1;
                % Check: must we make Spike_Times larger?
                if (NumSpikes(nn) > maxdimSp)
                    maxdimSp = maxdimSp + 1000;
                    Spike_Times = [Spike_Times zeros(N,1000)];
                end
                Spike_Times(nn,NumSpikes(nn)) = current_time;
                index = ceil(current_time/T_binning);
                %if more than one spike per bin, record how many time it
                %occurs
                if (Spikes(index,nn) > 0)
                    counter_bins = counter_bins+1;
                end
                Spikes(index,nn) = 1;
                
                %Begin using stereotyped voltage
                % Find appropriate place in spike to begin
                %       (If projected voltage is above zero, start near zero)
                
                projind=find(SpikeShapeNew >= min(V(nn),0));projind = projind(1);
                
                in_spike(nn) = projind;
                V(nn) = SpikeShapeNew(in_spike(nn));
                
                % OLD RESET
                %V(nn) = Vreset;             
           end
           
           % Even if we're not allowed to spike, we must 
           % check for an unacceptable voltage level: if present, reject
           % voltage step
           if (abs(V(nn))>1000)
                V(nn) = Vold(nn);
           end

        elseif (in_spike(nn) > 0)
            in_spike(nn) = in_spike(nn)+1;
            V(nn) = SpikeShapeNew(in_spike(nn));
            
            % are we done with the stereotyped spike train? i
            % if so, reset back to 0
            if (in_spike(nn) == max_in_spike)
                in_spike(nn)=0;
            end
            
        else 
            disp('Error, in_spike not legal value!')
           
        end %the if(NumSpikes .... 
        
    end %the nn loop 
 
%     %     Keep some portion of voltage trace
%      if (StimParam.ID_num == 1 && current_time > 1500 && current_time < 2000)
%         Vkeep = [Vkeep V];
%         tkeep = [tkeep current_time];
%         gexckeep = [gexckeep g_exc(:,tt-1)];
%         ginhkeep = [ginhkeep g_inh(:,tt-1)];
%     end
    % Reset Vold
    Vold = V;
end %the dt loop

% truncate spike_times array
nn = max(NumSpikes);
Spike_Times = Spike_Times(:,1:nn);

% fname = sprintf('ST_Recurrent_%d_%4.3g_%d_%d.mat',marg_flag,stim_std,t_refresh,StimParam.ID_num);
% save(fname,'Spike_Times','counter_bins');
% 
% if (StimParam.ID_num == 1)
% % % Save voltage trace, etc..
% fname = sprintf('test_recur_mod2.5_g=%d_sim_%d_%g_%g.mat',StimParam.gapg,marg_flag,stim_std,t_refresh);
% save(fname,'Vkeep','tkeep','gexckeep','ginhkeep','counter_bins','Spike_Times','Spikes');
% end

nout = max(nargout,2)-2;
if (nout > 0)
    % Pass back StimParam, in case changes have been made here
    varargout(1) = {StimParam};
end
    