function [f_samp,t_samp,v_samp] = propagate_pulse(mod_mat,time_mod,vel_mod,L,tau_char)
% Calculate interpolited time-velocity profile at sample position
% 
% mod_mat -- 2D moderator function in units tau(mks) vs 
% time_mod - time axis at moderator  in chop opening
% vel_mod  - velocity 
% mod_Energy - energy at moderator array (in mEv)
% tau -- time at choper to shift function around (in chopper opening time
% units)

% Output:
% Interpolated fime-velocity profile at sample position
% t_samp -- time axis for the profile above (in units of chopper opening
%           time)
% v_samp -- velocity axis for the profile above (in m/s)



tau_shift = (L/tau_char)./vel_mod; %

tau_at_sample = bsxfun(@plus, time_mod, tau_shift'); % t+L/v;

non_zero = mod_mat>0;
t_min = min(min((tau_at_sample(non_zero))));
t_max = max(max((tau_at_sample(non_zero))));

v_at_sample=reshape(repmat(vel_mod,1,numel(time_mod)),numel(vel_mod),numel(time_mod));
v_min = min(min((v_at_sample(non_zero))));
v_max = max(max((v_at_sample(non_zero))));

dt = (t_max-t_min)/(numel(time_mod)-1);
t_samp = t_min:dt:t_max;
dv = (v_max -v_min)/(numel(vel_mod)-1);
v_samp = v_min:dv:v_max;

[xb,yb]=meshgrid(time_mod,vel_mod);
[xi,yi]= meshgrid(t_samp,v_samp);
xi = xi - (L/tau_char)./yi; % t_samp-L/v;


f_samp = interp2(xb,yb,mod_mat,xi,yi,'nearest',0);

