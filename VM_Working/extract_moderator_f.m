function [f_mod,t_chop,v_chop] = extract_moderator_f(in_mat,time_mod,mod_Energy,tau,L,tau_char,nip)
% in_mat -- 2D moderator function in units tau(mks) vs energy (mEv)
% time_mod - time axis at moderator  in mks
% lambda_mod - wavelength axis at moderator (in A)
% mod_Energy - energy at moderator array (in mEv)
% tau -- time at choper to shift function around (in chopper opening time
% units)
% Output:
% Interpolated fime-velocity profile at choper position
% t_chop -- time axis for the profile above (in units of chopper opening
%           time)
% v_chop -- velocity axis for the profile above (in m/s)

dt = 2/(nip-1);

t_chop = tau-1:dt:tau+1; % chopper opening time in units of characteristic time

En_transf = 437.384; % conversion from root of energy (in mEv) to neutron velocity m/sec
%v_transf = 3.962e+3; % m/s (transform wavelength (A) into velocity)
V_mod = En_transf*sqrt(mod_Energy); % velocity at moderator (m/sec);
%V_mod = v_transf/lambda_mod; % velocity at moderator (m/sec);

tau_mod = (L/tau_char)./V_mod; %

tau_at_chop = bsxfun(@minus, t_chop, tau_mod'); % t-L/v;
time_mod    = time_mod*(1.e-6/tau_char); % time at moderator 
tau_mod_min = min(time_mod);
tau_mod_max = max(time_mod);

non_zero = tau_at_chop>=tau_mod_min & tau_at_chop<=tau_mod_max;

nL = numel(V_mod);
V_mod_w   = reshape(repmat(V_mod,1,nip),nL,nip); % moderator velocity axis
v_chop    = V_mod_w(non_zero);
v_ch_min  = min(v_chop);
v_ch_max  = max(v_chop);
if v_ch_min == v_ch_max
    v_ch_min = v_ch_min-1;
    v_ch_max = v_ch_max+1;    
    dV = 2/(0.5*nip-1);
else
    dV = (v_ch_max-v_ch_min)/(0.5*nip-1);
end
v_chop = v_ch_min:dV:v_ch_max;

[xb,yb]=meshgrid(time_mod,V_mod);
[xi,yi]= meshgrid(t_chop,v_chop);
xi = xi - (L/tau_char)./yi;

f_mod = interp2(xb,yb,in_mat',xi,yi,'nearest',0);

% convert velocity m/s to mEv
%vSq2mEv = 5.227e-6;  % s^2/m^2
%v_chop = v_chop.^2*vSq2mEv;






