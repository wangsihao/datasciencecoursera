function stim_struct = make_stim_struct()

% Some parameters needed for all stim types
%
% type_stim = 0 for noise only
%            1 for full field
%            2 for stixels
%
% if type_stim = 1 or 2, then need
%      marg_flag = 1 for binary, 
%              0 for Gaussian inputs
%      t_refresh = refresh time (in ms)
%      dt = time bin (in ms)
%      stim_std = standard deviation of stimulus <s^2>-<s>^2 
%           (NOTE: <s>=0 by assumption) - given in normalized units (1=4000
%           lumens - MUST CHECK)
%   
% 
stim_struct = struct('type_stim',0,...
    'marg_flag',0,...
    't_refresh',10,'dt',1,...
    'stim_std',1,...
    'stixel_size',60,...
    'random_flag',1,...
    'rot_par',0,'xshift',0,'yshift',0,...
    'N',3,...
    'ID_num',0,...
    'save_gexc_flag',0);
% if type_stim = 2, then need
%      
%    stixel size (in um)
%    random_flag = 1 pick offset, rotation parameters randomly
%             2    use arguments given in this structure
%         rot_par (given as a complex exponential argument)
%         xshift (um), yshift (um)

% N = number of cells
    
% ID_num = optional stim ID number
% save_gexc_flag = save conductances?
    