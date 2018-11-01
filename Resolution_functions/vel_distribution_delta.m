function [dv,f_d,v_peaks,dV_scat_max] = vel_distribution_delta(dvs)
% Calculate sample probabilityh distribution for velocity transfer 
% for modelling and recovering
% 
% The range of the energy to calculate the delta-function:
de_range = 0.4; %25;
% energy peaks
%e_exc = [0,0.3,2,10,15];
e_exch = 0.3;
e_min = min(e_exch);
e_max = max(e_exch);


e_transf_const = 5.22725e-6; % sec^2/m^2
dV_scat_max = sqrt((e_max+de_range)/e_transf_const);
%dV_scat_min = sqrt((e_min-de_range)/e_transf_const);
if numel(dvs) ==  1 % step provided    
    [dv,dV,NpDV] = adjust_step(-dV_scat_max,dV_scat_max,dvs);
else   % velocity points are provided
    dv = dvs;
    dV = dvs(2) - dvs(1);
end
dv_min = min(dv);

e_transf_const = 5.22725e-6; % sec^2/m^2
v_peaks  = sqrt(e_exch/e_transf_const);

f_d = zeros(size(dv));
for i=1:numel(e_exch)
    peak_ind = floor((v_peaks(i)-dv_min)/dV);
    f_d(peak_ind) = 1/dV; 
end

%plot(dv.*dv*(e_transf_const),f_d);

