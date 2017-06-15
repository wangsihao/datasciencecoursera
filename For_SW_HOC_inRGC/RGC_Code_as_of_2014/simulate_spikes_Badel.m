function [Spikes, counter_bins, varargout] = simulate_spikes_Badel(total_time,dt,StimParam)

%  OLD PARAMETERS - NOW PASSED AS STIM_STRUCT
%(total_time,dt,type_stim,marg_flag,t_refresh,stim_std,stixel_size_test,varargin)
%
% uses the integrate-and-fire model with the Badel fit
% to generate probabilities distributions
%
% type_stim can be 
% 0: noise only  
% 1: full field
% 2. stixels
% 
% marg_flag = 1 for binary, 0 for Gaussian
%
% JG - April 2010 mod Aug 2010
% AKB - mod June 2011 to take general stim structure
%    - mod June 2011 to have corrected Vreset, V_th, treatment of t_ref,
%    etc..

N = StimParam.N;
%N = 3; %number of neurons

scaryflag = zeros(N,1);

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

T_binning =   5 %3*tau_abs;              %width of a bin
Tstop = ceil(total_time/T_binning); %number of steps (bins)
Spikes = 0*ones(Tstop,N);

%(each column has num_steps entries for each neuron 1,...,N) 
%1's to be added later if neuron nn (1 to N) spikes at time tt (1:duration)
NumSpikes = zeros(1,N);      %keeps a record of the number of spikes for each neuron
maxdimSp = 1000;
Spike_Times = zeros(N,maxdimSp);

%Neuron properties
Vreset = -55; %mV
V_th   = 0;
counter_bins = 0;

V = Vreset*ones(N,1);  % a voltage vector for all the N neurons     
Vold = V;

% Set parameters for reversal potentials
% (***) actually, didn't I get these straight from Fred?
%     Took these from fitting of RGC data using procedure of Badel, et al.
VE = 0;
VI = -80;   % -80 mV

%correlation coefficients for exc and inh conductances
%used in the generation of noise conductances
c_exc = 0.30;
c_inh = 0.15;

%mean conductances i.e. tonic drive
G_exc = 30*0.1611; %nS
G_inh = 40*0.1634; %nS

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

    %parameters of filter f (1 for exc and 2 for inh)
    A = [-8e4 -1.8e5]/10^3;
    n = [3.6 3];
    tau = [12 16]; %ms
    T = [105 120]; %ms
    
    filter_exc = A(1)*(t_list_filter/tau(1)).^n(1).*exp(-(t_list_filter)/tau(1)).*sin(2*pi*t_list_filter/T(1));
    filter_inh = A(2)*(t_list_filter/tau(2)).^n(2).*exp(-(t_list_filter)/tau(2)).*sin(2*pi*t_list_filter/T(2));
    
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
    
    %parameters of filter f (1 for exc and 2 for inh)
    A = [-8e4 -1.8e5]/10^3;
    n = [3.6 3];
    tau = [12 16]; %ms
    T = [105 120]; %ms
    
    % Added shorter filter length - 10/21
    filter_exc = A(1)*(t_list_filter/tau(1)).^n(1).*exp(-(t_list_filter)/tau(1)).*sin(2*pi*t_list_filter/T(1));
    filter_inh = A(2)*(t_list_filter/tau(2)).^n(2).*exp(-(t_list_filter)/tau(2)).*sin(2*pi*t_list_filter/T(2));
    
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
if (StimParam.save_gexc_flag) % && random_flag==2)  
    if (type_stim == 2)
        Stixel_Size = StimParam.stixel_size;
    else
        Stixel_Size=0;
    end
    ID_number = StimParam.ID_num;
    if (marg_flag ==1)
        outfname = sprintf('gexc_bin_M_stixel=%d_ID=%d.mat',Stixel_Size,ID_number);
        save(outfname,'stim','g_exc','g_inh','StimParam');
    else
        outfname = sprintf('gexc_gauss_M_stixel=%d_ID=%d.mat',Stixel_Size,ID_number);
        save(outfname,'stim','g_exc','g_inh','StimParam');
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

% % % Keep some information
% total_time_keep = 100;
% total_bin_keep = ceil(total_time_keep/dt);
% Vkeep = zeros(N,total_bin_keep);
% %tkeep = [];
% gexckeep = Vkeep;
% ginhkeep = Vkeep;
% tkeep = zeros(1,total_bin_keep);
% 
% next_to_fill = 1;

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
        %neuron, if not refractory
        if (NumSpikes(nn) == 0 || ((current_time - Spike_Times(nn,NumSpikes(nn))) > tau_abs));
            
            I_ion = g_exc(nn,tt-1)*(Vold(nn)-VE) + g_inh(nn,tt-1)*(Vold(nn)-VI);
         
            V(nn) = Vold(nn) + dt*(-I_ion/C + F);
        
           if (V(nn) > V_th)
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
                %reset voltage
                V(nn) = Vreset;
                
            end

        else %same evolution, but don't check for a spike
            I_ion = g_exc(nn,tt-1)*(Vold(nn)-VE) + g_inh(nn,tt-1)*(Vold(nn)-VI);
              
            V(nn) = Vold(nn) + dt*(-I_ion/C + F);
            
            % However, must check for an unacceptable voltage level:
            if (abs(V(nn))>1000)
                V(nn) = Vold(nn);
            end
           
        end %the if(NumSpikes .... 
        
        if (sum(scaryflag)==0 && abs(V(nn))>1000)
            scaryflag(nn)= 1;
            disp(sprintf('Scary flag: voltage %d beyond bounds!: t=%g',nn,current_time));
            current_time
            disp(V);
            
            % Wait for a while, and then quit
            time_to_quit = current_time+50;
        end
        if (current_time > time_to_quit)
            
            % Save what you can
            index = ceil(current_time/T_binning);
            Spikes = Spikes(1:index-10,:);
         
            % truncate spike_times array
            nn = max(NumSpikes);
            Spike_Times = Spike_Times(:,1:nn);
            
%             % % Save voltage trace, etc..
%             fname = sprintf('check_voltage_data_%d_%g_%d.mat',marg_flag,stim_std,StimParam.ID_num);
%             save(fname,'Vkeep','tkeep','gexckeep','ginhkeep','counter_bins','Spike_Times','Spikes');

            nout = max(nargout,2)-2;
            if (nout > 0)
                varargout(1)={StimParam};
            end
            return;
        end
        
    end %the nn loop 
% %      % Keep some portion of voltage trace
%     Vkeep(:,next_to_fill) = V;
%     tkeep(next_to_fill) = current_time;
%     gexckeep(:,next_to_fill) = g_exc(:,tt-1);
%     ginhkeep(:,next_to_fill) = g_inh(:,tt-1);
%     
%     if (next_to_fill == total_bin_keep)
%         next_to_fill = 1;
%     else
%         next_to_fill = next_to_fill+1;
%     end
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

%fname = sprintf('ST_%d_%4.3g_%d_%d.mat',marg_flag,stim_std,t_refresh,StimParam.ID_num);
%save(fname,'Spike_Times','counter_bins');

% % % Save voltage trace, etc..
% fname = sprintf('check_voltage_data_%d_%g_%d.mat',marg_flag,stim_std,StimParam.ID_num);
% save(fname,'Vkeep','tkeep','gexckeep','ginhkeep','counter_bins','Spike_Times','Spikes');


nout = max(nargout,2)-2;
if (nout > 0)
    % Pass back StimParam, in case changes have been made here
    varargout(1) = {StimParam};
end
    