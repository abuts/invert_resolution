function [f_out,t_out,v_out] = propagate_pulse(f_in,time_in,vel_in,L,tau_char)
% Calculate interpolited time-velocity profile at the position L
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



tau_shift = (L/tau_char)./vel_in; %

tau_at_sample = bsxfun(@plus, time_in, tau_shift'); % t+L/v;

non_zero = f_in>0;
t_min = min(min((tau_at_sample(non_zero))));
t_max = max(max((tau_at_sample(non_zero))));

v_at_end=reshape(repmat(vel_in,1,numel(time_in)),numel(vel_in),numel(time_in));
v_min = min(min((v_at_end(non_zero))));
v_max = max(max((v_at_end(non_zero))));

dt = (t_max-t_min)/(numel(time_in)-1);
t_out = t_min:dt:t_max;
dv = (v_max -v_min)/(numel(vel_in)-1);
v_out = v_min:dv:v_max;

[xb,yb]=meshgrid(time_in,vel_in);
[xi,yi]= meshgrid(t_out,v_out);
xi = xi - (L/tau_char)./yi; % t_samp-L/v;


f_out = interp2(xb,yb,f_in,xi,yi,'nearest',0);

