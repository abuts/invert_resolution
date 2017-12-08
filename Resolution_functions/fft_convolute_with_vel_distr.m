function [f_s_out,t_samp,v_s_out] = fft_convolute_with_vel_distr(f_samp,t_samp,v_samp,V_char)
% Calculate sample velocity distribution for modelling and recovering
%
%
%   Detailed explanation goes here




e_max = 25;

e_transf_const = 5.22725e-6; % sec^2/m^2
dV_max = sqrt(e_max/e_transf_const);

%v_samp_norm = v_samp/V_char;
dV = v_samp(2) - v_samp(1);
V_min = min(v_samp)-dV_max;
V_max = max(v_samp)+dV_max;


v_s_out = V_min:dV:V_max;
Nv = numel(v_s_out);
Nt = size(f_samp,2);

ft = fft2(f_samp,Nv,Nt);

%ftm = abs(fftshift(ft));
%surf(ftm,'EdgeColor','none');


% caclculate velocity transfer distrubution in its own range but with
% joing accuracy.
[v_t,f_d] = vel_distribution(-dV_max:2*dV_max/(Nv-1):dV_max);
figure(111)
acolor 'b';
plot(v_t/V_char,f_d);
ax = gca;
ax.XLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
ffv = fft(f_d);



fm = bsxfun(@times,ft',ffv);
%fm = ft.*repmat(ffv',1,numel(t_samp));
%fm = fftshift(fm);
%ftm = abs(fftshift(fm));
%surf(ftm,'EdgeColor','none');


%fm = repmat(fm(1,:),size(fm,1),1);
f_s_out = ifft2(fm)';





function fc = f1_kernel(dv,v_samp,f_samp,v_s_out,V_char)

Mf = @(v)(interp1(v_samp,f_samp,v,'linear',0));

[~,fc]= vel_distribution(dv,V_char);
fc = fc.*Mf(v_s_out-dv);



