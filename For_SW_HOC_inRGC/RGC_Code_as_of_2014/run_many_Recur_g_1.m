% Full-field simulations

% For testing counter_bins, for when I'm ready to do stixel sims:
% numsims = 1;
% stim_std = [1/2];
% refresh_size = [8];

% For real sims
numsims = 20;
% 1/3 and decrease by 1/2...
% 1/2 and decrease by 1/2...
stim_std = [1/2,1/3,1/4,1/6,1/8,1/12,1/16];
refresh_size = [8 40 100];


prob_gauss = cell(length(refresh_size),length(stim_std));
prob_bin = cell(length(refresh_size),length(stim_std));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
% Encode (most) of these parameters into a StimStruct
StimParam = make_stim_struct;
StimParam.N = 3;
StimParam.type_stim = 1;
% StimParam.marg_flag = 0;
% StimParam.t_refresh = refresh_size;
StimParam.dt = 0.05;
% StimParam.stim_std = ;
% StimParam.stixel_size = ;
% StimParam.random_flag = 1;
% StimParam.rot_par = 0;
% StimParam.xshift = 0;
% StimParam.yshift = 0;
% Save conductances
StimParam.save_gexc_flag = 0;
StimParam.gapg = 1;

StimParam

for jj=1:length(refresh_size)
        refresh_size_i = refresh_size(jj);
        refresh_size_i
        
    for kk=1:length(stim_std)
        stim_std_i = stim_std(kk);
        stim_std_i

        prob_d_vec_g = [];
        prob_d_vec_b = [];
        
        StimParam.t_refresh = refresh_size_i;
        StimParam.stim_std = stim_std_i;
        
        %run 20 runs with this particular set of parameters, and stim std
        for mm=1:numsims
            StimParam.ID_num=mm;
            mm

            % Gaussian
            StimParam.marg_flag = 0;
            prob_d = run_IaF_Badel_Recur_Rect(StimParam);
            prob_d_vec_g = [prob_d_vec_g prob_d];

            % Binary
            StimParam.marg_flag = 1;
            prob_d = run_IaF_Badel_Recur_Rect(StimParam);
            prob_d_vec_b = [prob_d_vec_b prob_d];
            
            prob_gauss{jj,kk} = prob_d_vec_g;
            prob_bin{jj,kk} = prob_d_vec_b;
            
            save run_FF_P_Recur_g_1.mat refresh_size stim_std prob_gauss prob_bin
            %save FullField_061711.mat refresh_size stim_std prob_gauss prob_bin
        end
        
    end
end
