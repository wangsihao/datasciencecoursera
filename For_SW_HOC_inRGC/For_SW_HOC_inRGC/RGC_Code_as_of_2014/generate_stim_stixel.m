function [F,varargout] = generate_stim_stixel(t_list,marg_flag,Smean,Sstd,StimParam) 
% F: stimulus trains for N cells
% (Optional): varargout{1} = pass back StimParam
%
% function F = generate_stim_stixel(t_list,flag,Smean,Sstd,t_refresh) 
%
%generates stimulus given a flag: 
%if marg_flag=0 then Gaussian
%if marg_flag=1 then binary
%re-written from Soft_fun_stixel.m
carray  = ['b';'m';'g';'y';'r'];

Stixel_Size = StimParam.stixel_size;
t_refresh = StimParam.t_refresh;
random_flag = StimParam.random_flag;

if (random_flag==2)
    rot_par = StimParam.rot_par;
    xshift = StimParam.xshift;
    yshift = StimParam.yshift;
    if (StimParam.ID_num>0)
        ID_number = StimParam.ID_num;
        newID = floor(ID_number/20+0.5);
    end
elseif (random_flag==1)
    % NOTE: THIS WILL BE REDEFINED LATER!!
    dx = 4;num  = floor(Stixel_Size/dx+0.5);

    % Choose parameters randomly here
    rot_par = exp(sqrt(-1)*2*pi*rand);
    xshift = floor(num*rand);
    yshift = floor(num*rand);
    % Save these parameters in Stim Struct, in case we want to remember
    % them later
    StimParam.rot_par = rot_par;
    StimParam.xshift= xshift;
    StimParam.yshift= yshift;
else
    % No shift/rotation
    xshift = 0;yshift=0;rot_par = 1;
end

if marg_flag == 1 %binary
    Smax = Smean + Sstd;
    Smin = Smean - Sstd; 
end

%%%%%%%%%%% STIXEL %%%%%%%%%%%%%%%%%%%%%

% Build spatial grid IN MICROMETERS
dx=4;
dy=4;
xa = -200:dx:200; 
ya = -200:dy:200;
RF1 = zeros(numel(xa),numel(ya));
RF2 = RF1;
RF3 = RF1;

% L  = Distance between centers
% RC = center radius
% RS = surround radius
% RelS  = relative strength of surround
RC=50;
RS=80;  
L = RC*(2/sqrt(3)); %Distance so that RF will "just touch" 
RelS = .3;

% characterizes ellipticity, orientation of RF
Q = [1 0;0 1];

c1x = 0;
c1y=L;

c2x = -sqrt(3)*L/2;
c2y=-L/2; 

c3x = sqrt(3)*L/2;
c3y=-L/2;


% rotate receptive fields by a random parameter
c1_C = c1x + sqrt(-1)*c1y;
c2_C = c2x + sqrt(-1)*c2y;
c3_C = c3x + sqrt(-1)*c3y;


c1_C = c1_C*rot_par;  
c2_C = c2_C*rot_par;  
c3_C = c3_C*rot_par;

c1x = real(c1_C);c1y = imag(c1_C); 
c2x = real(c2_C);c2y = imag(c2_C); 
c3x = real(c3_C);c3y = imag(c3_C);

% Create receptive fields
for j=1:numel(xa)
    for k=1:numel(ya)
        dv = [xa(j)-c1x ya(k)-c1y];
        d2 = dv*Q*dv';RF1(j,k) = exp(-.5*d2/RC^2)- RelS*exp(-.5*d2/RS^2);
        dv = [xa(j)-c2x ya(k)-c2y];
        d2 = dv*Q*dv';RF2(j,k) = exp(-.5*d2/RC^2)- RelS*exp(-.5*d2/RS^2);
        dv = [xa(j)-c3x ya(k)-c3y];
        d2 = dv*Q*dv';RF3(j,k) = exp(-.5*d2/RC^2)- RelS*exp(-.5*d2/RS^2);
    end
end

% Normalize receptive fields
RF1 = RF1/(sum(sum(RF1))*dx*dy);
RF2 = RF2/(sum(sum(RF2))*dx*dy);
RF3 = RF3/(sum(sum(RF3))*dx*dy);

% Get grid
% SIZE OF STIXEL IN MICROMETERS
% STARTS AT DX
stim_grid=zeros(size(RF1));
num  = floor(Stixel_Size/dx+0.5);

% find "first and last" of each x,y block
xbl_first = [1+xshift];
xbl_last = [min(numel(xa),xshift+num)];

curr_ind = 1+xshift;

while (curr_ind+num <=numel(xa))
    curr_ind = curr_ind+num;
    xbl_first = [xbl_first curr_ind];
    xbl_last  = [xbl_last min(numel(xa),curr_ind+num-1)];
end

if (xshift > 0)
    xbl_first = [1 xbl_first];
    xbl_last = [xshift xbl_last];
end
    
ybl_first = [1+yshift];
ybl_last = [min(numel(ya),yshift+num)];
curr_ind = 1+yshift;

while (curr_ind+num <=numel(ya))
    curr_ind = curr_ind+num;
    ybl_first = [ybl_first curr_ind];
    ybl_last  = [ybl_last min(numel(ya),curr_ind+num-1)];
end

if (yshift > 0)
    ybl_first = [1 ybl_first];
    ybl_last = [yshift ybl_last];
end

xdir = numel(xbl_first);
ydir = numel(ybl_first);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% takes vector argument t_list in ** SEC **, produces list the same length of
% binary values that switch randomly between Smax and Smin every t_refresh
% ** SEC ** 

t_length = length(t_list);
F = zeros(3,t_length);
dt = t_list(2)-t_list(1);
Tmax = t_list(end) - t_list(1) + dt;


for jt=1:t_length
    
    if (jt==1 || t_list(jt)-stim_time > t_refresh)
        % Fill up stimulus grid anew; otherwise, just repeat
    
        for j=1:xdir
            for k=1:ydir
                if marg_flag == 1
                    r=round(rand);
                    stim = r*Smax + (1-r)*Smin;  %takes either value w/ equal proba 
                else
                    stim = Smean + Sstd*randn;  %takes either value w/ equal proba 
                end
                
%                 % If I just want checkboards
%                 if (flag==1)
%                     if (mod(j+k,2)==0)
%                         stim = Smax;
%                     else
%                         stim = Smin;
%                     end
%                 else
%                     if (mod(j+k,2)==0)
%                         stim = Smean+Sstd;
%                     else
%                         stim = Smean-Sstd;
%                     end
%                 end
                stim_grid(xbl_first(j):xbl_last(j),ybl_first(k):ybl_last(k))=stim;
            end
        end
             
        % Compute Soft_list1, Soft_list2, Soft_list3, etc..
        stim1 = sum(sum(stim_grid.*RF1))*dx*dy;
        stim2 = sum(sum(stim_grid.*RF2))*dx*dy;
        stim3 = sum(sum(stim_grid.*RF3))*dx*dy;
    
        stim_time = t_list(jt);
    end
    
    F(1,jt) = stim1;
    F(2,jt) = stim2;
    F(3,jt) = stim3;
    
    
    plot_RF_yes = 0;
    if (jt==1 && plot_RF_yes)
        min(min(stim_grid-1))
        max(max(stim_grid-1))
        figure;surf(xa,ya,(stim_grid-1)');axis equal;shading flat;
        caxis([-.75-1 .75-1]);colormap gray;view(0,90);axis off;hold;
        
        % Find SD of gaussian: 
        p1  = max(max(RF1));
        SD1 = exp(-1)*p1;
        SD2 = exp(-2)*p1;
        con_val = [0 SD2 (SD1+SD2)/2 SD1 (SD1+p1)/2];
        %con_val_neg = [-SD2 -SD2/5];
        
        figure;hold;axis equal;view(0,90);axis off;
        contour(xa,ya,RF1',con_val,'c','LineWidth',2);
        %contour(xa,ya,RF1',con_val_neg,'m');
        contour(xa,ya,RF2',con_val,'c','LineWidth',2);
        contour(xa,ya,RF3',con_val,'c','LineWidth',2);
        text(c1x*1.6,c1y*1.6,'1','FontSize',24,'Color',carray(newID,:));
        text(c2x*1.6,c2y*1.6,'2','FontSize',24,'Color',carray(newID,:));
        text(c3x*1.6,c3y*1.6,'3','FontSize',24,'Color',carray(newID,:));
        
        set(gca,'FontSize',16);
        %title(sprintf('RF for stixel size=%d, run=%d  ',Stixel_Size,newID));
    end   
end

nout = max(nargout,1)-1;
if (nout > 0)
    % Pass back StimParam, in case changes have been made here
    varargout(1) = {StimParam};
end
    


    