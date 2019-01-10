function [f_out,t_out] = propagate_pulse_Int(f_in,time_in,vel_in,L,V_pulse,t_char,v_char)
% Calculate interpolited time-velocity profile at the position L using fft
%
% f_in     -- 2D signal function in units tau(mks) vs
% time_in  -- time axis for signal  (sec)
% vel_in   -- velocity accis for signal (m/sec)
% L        -- distance to the target
% tau -- time at choper to shift function around (in chopper opening time
%
% Output:
%  Interpolated fime-velocity profile at position L
% f_out  --
% t_out -- time axis for the profile above (sec)
% v_out -- velocity axis for the profile above (in m/s)
v_min = min(vel_in);
v_max = max(vel_in);

tp_min = min(time_in);
tp_max = L/v_min+max(time_in);
dt = time_in(2) - time_in(1);


t_out = tp_min:dt:tp_max;
Nt = numel(t_out);
f_out = zeros(1,Nt);

for i=1:Nt
    fv = @(v)f_tv(f_in,time_in,vel_in,L,t_out(i),v);
    f_out(i) = integral(fv,v_min,v_max);
    if rem(i,100) == 0
        fprintf('step %d#%d\n',i,Nt);
    end
end



function f= f_tv(f_in,time_in,vel_in,L,t0,vel)

persistent xb;
persistent yb;
if isempty(xb)
    [xb,yb] = meshgrid(time_in,vel_in);
end
ti = t0-L./vel;

f = interp2(xb,yb,f_in,ti,vel,'linear',0);

