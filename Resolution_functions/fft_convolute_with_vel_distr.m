function [f_s_out,t_samp,v_s_out] = fft_convolute_with_vel_distr(f_samp,t_samp,v_samp,tau_char,V_char)
% Calculate sample velocity distribution for modelling and recovering
%
%
%   Detailed explanation goes here




e_max = 25;

e_transf_const = 5.22725e-6; % sec^2/m^2
dV_max = sqrt(e_max/e_transf_const)/V_char;

v_samp_norm = v_samp/V_char;
dV = v_samp_norm (2) - v_samp_norm(1);
V_min = min(v_samp_norm)-dV_max;
V_max = max(v_samp_norm)+dV_max;


v_s_out = V_min:dV:V_max;
Nv = numel(v_s_out);
Nt = size(f_samp,2);

ft = fft2(f_samp,Nv,Nt);

ftm = abs(fftshift(ft));
surf(ftm,'EdgeColor','none');


% caclculate velocity transfer distrubution in its own range but with
% joing accuracy.
[v_t,f_d] = vel_distribution(-dV_max:2*dV_max/(Nv-1):dV_max,V_char);
figure
plot(v_t,f_d);
ffv = fft(f_d);



fm = bsxfun(@times,ft',ffv);
%fm = ft.*repmat(ffv',1,numel(t_samp));
%fm = fftshift(fm);
ftm = abs(fftshift(fm));
surf(ftm,'EdgeColor','none');


%fm = repmat(fm(1,:),size(fm,1),1);
f_s_out = ifft2(fm);

% f_s_out = f_samp;
% figure;
% for i=1:numel(t_samp)
%     conv_profile = zeros(size(v_s_out));
%     for j=1:numel(v_samp)
%         v_int = 0;
%         if sum(f_samp(:,i)) ~= 0
%             for k=1:numel(vI_range)-1
%                 v_int = v_int+integral(@(dv)f1_kernel(dv,v_samp_norm,f_samp(:,i) ,v_s_out(j),V_char),vI_range(k),vI_range(k+1));
%             end
%             %            v_int = integral(@(dv)f1_kernel(dv,v_samp_norm,f_samp(:,i) ,v_s_out(j),V_char),V_min-4*dV_max,V_max+4*dV_max);
%         end
%         conv_profile(j) = v_int;
%     end
%     f_s_out(:,i) = conv_profile;
%     plot(v_s_out.*v_s_out*(V_char*V_char*e_transf_const),conv_profile);
% end




function fc = f1_kernel(dv,v_samp,f_samp,v_s_out,V_char)

Mf = @(v)(interp1(v_samp,f_samp,v,'linear',0));

[~,fc]= vel_distribution(dv,V_char);
fc = fc.*Mf(v_s_out-dv);



