function [v,f_d,v_peaks] = vel_distribution(v)
% Calculate sample velocity distribution for modelling and recovering
% 
%
%   Detailed explanation goes here

% energy peaks
e_exc = [0,0.3,2,10,15];
sigma = [10,60,10,80,160];
Ampl  = [5,9,2,10,10];

e_transf_const = 5.22725e-6; % sec^2/m^2
v_peaks  = sqrt(e_exc/e_transf_const);

sig22 = 2*sigma.*sigma;
par = [Ampl',v_peaks',sig22'];

f_d = zeros(size(v));
for i=1:numel(e_exc)
    f_d = f_d +par(i,1).*(exp(-(v-par(i,2)).^2./par(i,3))+exp(-(v+par(i,2)).^2./par(i,3)))./sqrt(pi*par(i,3));
end



%plot(v.*v*(V_char*V_char*e_transf_const),f_d);







