function F = generate_stim(t_list,flag,Smean,Sstd,t_refresh) 
%function F = generate_stim(t_list,flag,Smean,Sstd,t_refresh) 
%
%generates stimulus given a flag: 
%if flag=0 then Gaussian
%if flag=1 then binary
%JG: re-written from Soft_fun.m

if flag == 1
    Smax = Smean + Sstd;
    Smin = Smean - Sstd; 
end


% takes vector argument t_list in ** SEC **, produces list the same length of
% binary values that switch randomly between Smax and Smin every t_refresh
% ** SEC ** 

t_length = length(t_list);
F = zeros(size(t_list));
dt = t_list(2)-t_list(1);
Tmax = t_list(end) - t_list(1) + dt;

%do a crude loop through t_list, setting stimulus values ...
stim_index = 1;
stim_start = t_list(1);

%initial value of stim
if flag == 1
    r = round(rand);
    stim = r*Smax + (1-r)*Smin;  %takes either value w/ equal proba 
else
    stim = Smean + Sstd*randn;  %takes either value w/ equal proba 
    if (stim < -1)
        % resample 
        while (stim < -1)
           stim = Smean + Sstd*randn;
        end
        %stim = -1;
    end
end

stim_time = t_list(1);

for jj=1:t_length

    if t_list(jj)-stim_time > t_refresh
        if flag == 1
            r = round(rand);
            stim = r*Smax + (1-r)*Smin;  %takes either value w/ equal proba 
        else
            stim = Smean + Sstd*randn;  %takes either value w/ equal proba 
        end
        stim_time = t_list(jj);
    end
    
    F(jj) = stim;
    
end

    