function [f_out,dv_out] = fft_invert_propagation(f_samp,t_samp,v_samp,f_det_vs_t,t_det,L_det,v_max,t_char,v_char)
% retrieve samle gain/loss velocity distribution given signal on the detector and
% and the smple incident beam time/velocity distribution.

% shift sample time frame to 0
%t0 = min(t_samp);
%t_det = t_det-t0;
%t_samp_lf = t_samp-t0;
t_samp_lf = t_samp;

dt = t_det(2)-t_det(1);
if dt<2e-6 % debugging -- not needed in reality
    dt = 2e-6;
end

tp_min = min(t_samp_lf)+L_det/v_max;
tp_max = max(t_det);
[t_det_r,dt,Nt] = adjust_step(tp_min,tp_max,dt);
DT0 = max(t_det_r)-min(t_det_r);


f_det = interp1(t_det,f_det_vs_t,t_det_r,'linear',0);
if ~exist('v_max','var')
    i1 = find(f_det,1);
    v_max = L_det/t_in(i1);
end
v_min = L_det/(tp_max-min(t_samp_lf));
 % expand sample pulse onto the whole detectir integration range
dv = 4*(v_samp(2)-v_samp(1));
%
%dv = (v_max-v_min)/(Nt-1);
[v_pulse,dv,Nv] = adjust_step(v_min,v_max,dv);

[xb,yb] = meshgrid(t_samp,v_samp);
[xi,yi] = meshgrid(t_det_r-L_det/v_max,v_pulse);
f_in = interp2(xb,yb,f_samp,xi,yi,'linear',0);



% %
% % % [xi,yi]= meshgrid(t_out,vel_in);
% % [xb,yb] = meshgrid(t_int,vel_in);
% surf(xi/t_char,yi/v_char,f_in,'Edgecolor','None');
% view(0,90);

% f_in_int = interp2(xb,yb,f_in,xi,yi,'linear',0);
%
%
% fft frequencies do not correspond to indexes as index m in fft array
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
% v_phase =  exp(-1i*pi*v_index*DNv/(Nv-1)); % appears to shift zeros added at the end to zeros, added to the beginning
% %
% f_in_sp = bsxfun(@times,f_in_sp,v_phase');
%f_in_sp = bsxfun(@times,f_in_sp,t_phase);
% f_in_t  = ifft2(f_in_sp);
% figure(21)
% [xi,yi] = meshgrid(t_in,dv_out);
% surf(xi,yi,abs(f_in_t),'Edgecolor','None');
% view(0,90);
f_in_t  = ifft2(f_in_sp);
figure(21)
[xi,yi] = meshgrid(t_det_r/t_char,v_pulse/v_char);
surf(xi,yi,abs(f_in_t),'Edgecolor','None');
view(0,90);



v_index = fft_ind(Nv);
t_index = fft_ind(Nt);
DV0 = max(v_pulse)-min(v_pulse);


v_min_norm = v_min/(DV0);
v_max_norm = v_max/(DV0);
betta0 = (L_det/DT0)/DV0;

%fun = @(m,n)I_mkv4(m,n);
f_con_mat  = build_dif_mat(betta0,v_min_norm,v_max_norm,'IntW',Nv,Nt);


f_in_sp = f_in_sp.*f_con_mat;
v_phase =  exp(1i*pi*v_index);
f_in_sp = bsxfun(@times,f_in_sp,v_phase');


f_det_sp = fft(f_det,Nt);
%
f_vel_sp  = linsolve((f_in_sp)',(f_det_sp)');


%t_phase = (DV0)*exp(-1i*pi*t_index*((tp_min+tp_max)/(2*tp_max))); %
%t_phase = (DV0)*exp(1i*pi*t_index); %
%f_vel_sp   = f_vel_sp .*v_phase;
%
f_out = ifft(f_vel_sp);
