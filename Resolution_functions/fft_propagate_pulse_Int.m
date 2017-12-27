function [f_out,t_out] = fft_propagate_pulse_Int(f_in,time_in,vel_in,L)
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
DV0 = v_max - v_min;
tp_min = L/v_max+min(time_in);
tp_max = L/v_min+max(time_in);
dt = time_in(2) - time_in(1);
% if dt<2e-6 % debugging -- not necessary in reality
%     dt = 2e-6;
% end
t_out = tp_min-dt:dt:tp_max;
Nt = numel(t_out);
Nv = numel(vel_in);


v_min_norm = v_min/DV0;
v_max_norm = v_max/DV0;

v_index = fft_ind(Nv);
t_index = fft_ind(Nt);

Imk = @(nv,kt)I_mk(nv,kt,v_min_norm,v_max_norm,v_index,t_index);


I1 = Imk(1,2);
%
[xi,yi]= meshgrid(t_out,vel_in);
xi = xi-L./yi;
[xb,yb] = meshgrid(time_in,vel_in);
f_in = interp2(xb,yb,f_in,xi,yi,'linear',0);


% fft frequencies do not correspond to indexes as index m in fft array
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
persistent f_con_mat;
if isempty(f_con_mat)
    cn = sprintf('CashNv%dNt%d',Nv,Nt);
    if exist([cn,'.mat'],'file')
        cns = load(cn);
        f_con_mat = cns.f_con_mat;
    else
        f_con_mat = ones(Nv,Nt)*NaN;
    end
end

t_start=tic;
call_count = 0;
for k=1:Nv
    undef = isnan(f_con_mat(k,:));
    if ~any(undef)
        continue;
    end
    for n=1:Nt
        if undef(n)
            f_con_mat(k,n) = Imk(k,n);
        end
    end
    if call_count > 10
        fprintf('k=%d#%d\n',k,Nv);
        cn = sprintf('CashNv%dNt%d',Nv,Nt);
        save(cn,'f_con_mat');
        call_count = 0;
    end
    call_count = call_count+1;
    
end
f_out_sp = sum(f_in_sp.*f_con_mat,1);
time_c = toc(t_start)/60 % convert in minutes
f_out = ifft(f_out_sp);



