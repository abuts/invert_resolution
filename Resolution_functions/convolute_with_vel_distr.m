function [f_s_out,t_samp,v_s_out] = convolute_with_vel_distr(f_samp,t_samp,v_samp,tau_char,V_char)
% Calculate sample velocity distribution for modelling and recovering
%
%
%   Detailed explanation goes here

%Np = 1000;
Np = numel(v_samp);
e_max = 25;

e_transf_const = 5.22725e-6; % sec^2/m^2
dV_max = sqrt(e_max/e_transf_const)/V_char;

v_samp_norm = v_samp/V_char;
V_min = min(v_samp_norm)-dV_max;
V_max = max(v_samp_norm)+dV_max;

dV = (V_max-V_min)/(Np-1);
figure;
v_s_out = V_min:dV:V_max;
[v_t,f_d,v_peaks] = vel_distribution(-dV_max:2*dV_max/(Np-1):dV_max,V_char);
plot(v_t,f_d);

% is_pos_peak = v_peaks>1.e-6;
% v_peaks_neg = -(v_peaks(is_pos_peak));
% v_peaks = sort([v_peaks_neg,v_peaks]);
% vI_min = min(v_peaks)-5*dV_max;
% vI_max = max(v_peaks)+5*dV_max;
% vI_range = [vI_min,v_peaks,vI_max];

f_s_out = f_samp;
for i=1:numel(t_samp)
    conv_profile = zeros(size(v_s_out));
    for j=1:numel(v_samp)
        v_int = 0;
        if sum(f_samp(:,i)) ~= 0
%             for k=1:numel(vI_range)-1
%                 v_int = v_int+integral(@(dv)f1_kernel(dv,v_samp_norm,f_samp(:,i) ,v_s_out(j),V_char),vI_range(k),vI_range(k+1));
%             end
            v_int = integral(@(dv)f1_kernel(dv,v_samp_norm,f_samp(:,i) ,v_s_out(j),V_char),V_min-4*dV_max,V_max+4*dV_max);
        end
        conv_profile(j) = v_int;
    end
    f_s_out(:,i) = conv_profile;
    plot(v_s_out.*v_s_out*(V_char*V_char*e_transf_const),conv_profile);
end




function fc = f1_kernel(dv,v_samp,f_samp,v_s_out,V_char)

Mf = @(v)(interp1(v_samp,f_samp,v,'linear',0));

[~,fc]= vel_distribution(dv,V_char);
fc = fc.*Mf(v_s_out-dv);



