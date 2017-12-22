function [f_out,v_out] = fft_invert_propagation(ModPulse,f_det_vs_t,t_det,L)
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
vel_in = L/t_det;
v_min  = min(vel_in);
v_max  = max(vel_in);
tp_min = min(time_in);
tp_max = max(time_in);

DV0 = v_max - v_min;
dt = time_in(2)-time_in(1);
%Nt = size(f_in,2);
%dt = (tp_max-tp_min)/(Nt-1);
% if dt<2e-6
%     dt = 2e-6;
% end
t_out = tp_min:dt:tp_max;
Nt = numel(t_out);



v_min_norm = v_min/DV0;
v_max_norm = v_max/DV0;

v_index = fft_ind(Nv);
t_index = fft_ind(Nt);

Imk = @(n,k)I_mk(n,k,v_min_norm,v_max_norm,v_index,t_index,Nv,Nt);

I1 = Imkl(0,0);
%


% fft frequencies do not correspond to indexes as index m in fft array 
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft(f_det_vs_t,Nt);
% indexes of frequncies corresponging to V-frequencies in f_in_sp;

t_start=tic;
f_out_sp = zeros(Nv,Nt);
for k=1:Nv
    fprintf('k=%d#%d\n',k,Nv);
    for l=1:Nt
        f_out_sp(k,l) = Imk(k,l);
    end
end
f_out_sp = f_in_sp(k,n).*f_out_sp;

time_c = toc(t_start)/60 % convert in minutes
f_out = ifft2(f_out_sp);

function int = I_mk(n,k,v_min,v_max,v_index,t_index,Nv,Nt)
persistent cash;
if isempty(cash)
    cash = ones(Nv,Nt)*NaN;
end
int = cash(n,k);
if isnan(int)
    ki = v_index(k);
    ni = t_index(k);    
    fun = @(v)(exp(2i*pi*(ki*v-v_min*ni./v)));
    int = integral(fun,v_min,v_max);
    cash(n,k) = int;
end



