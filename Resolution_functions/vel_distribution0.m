function [dv,f_d,v_peaks,dV_scat_max] = vel_distribution0(dvs)
% Calculate sample probabilityh distribution for velocity transfer 
% for modelling and recovering
% 
%
e_max = 25;

e_transf_const = 5.22725e-6; % sec^2/m^2
dV_scat_max = sqrt(e_max/e_transf_const);
if numel(dvs) ==  1 % step provided
    [dv,dV,NpDV] = adjust_step(-dV_scat_max,dV_scat_max,dvs);
else   % velocity points are provided
    dv = dvs;
end

% energy peaks
e_exc = [0,0.3,2,10,15];
sigma = [10,60,5,80,160];
%sigma = [50,300,50,400,160];
Ampl  = [5,9,1,5,5];

e_transf_const = 5.22725e-6; % sec^2/m^2
v_peaks  = sqrt(e_exc/e_transf_const);

sig22 = 2*sigma.*sigma;
par = [Ampl',v_peaks',sig22'];

f_d = zeros(size(dv));
for i=1:numel(e_exc)
    f_d = f_d +par(i,1).*(exp(-(dv-par(i,2)).^2./par(i,3))+exp(-(dv+par(i,2)).^2./par(i,3)))./sqrt(pi*par(i,3));
end

%plot(v.*v*(V_char*V_char*e_transf_const),f_d);

