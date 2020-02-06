function [f_out,t_out,v_out] = propagate_pulse(f_in,time_in,vel_in,L,t_exp)
% Calculate interpolited time-velocity profile at the position L
%
% f_in   -- 2D moderator function in units tau(mks) vs
% time_in - time axis at moderator  (sec)
% time_in  - velocity
% mod_Energy - energy at moderator array (in mEv)
% tau -- time at choper to shift function around (in chopper opening time
% units)

% Output:
% Interpolated fime-velocity profile at sample position
% t_samp -- time axis for the profile above (in units of chopper opening
%           time)
% v_samp -- velocity axis for the profile above (in m/s)

if ~exist('t_exp','var')
    t_exp = 0;
end

tau_shift = L./vel_in; %


%tau_final = bsxfun(@plus, time_in, tau_shift'); % t+L/v;


add_time_max = max(tau_shift); % additional time, the slowest particle achieve the detector;
add_time_min = min(tau_shift); % the time max velocity arrive at the detector

t_min = min(time_in)+add_time_min;
t_max = max(time_in)+add_time_max; % the time the slowest particle recorded at the end ot the
% sample time frame arrives at the detector

non_zero = f_in>0;
% t_min = min(min((tau_final(non_zero))));
% t_max = max(max((tau_final(non_zero))))+t_exp;
% 
v_at_end=reshape(repmat(vel_in,1,numel(time_in)),numel(vel_in),numel(time_in));
v_min = min(min((v_at_end(non_zero))));
v_max = max(max((v_at_end(non_zero))));

dt = (t_max-t_min)/(max(numel(time_in),numel(vel_in))-1);
t_out = t_min:dt:t_max;
dv = (v_max -v_min)/(numel(vel_in)-1);
if dv == 0
    v_out = vel_in;
else
    v_out = v_min:dv:v_max;
end

% Regrid on square grid;
[xb,yb]=meshgrid(time_in,vel_in);
[xi,yi]= meshgrid(t_out,v_out);
xi = xi - L./yi; % t_samp-L/v;


f_out = interp2(xb,yb,f_in,xi,yi,'nearest',0);

