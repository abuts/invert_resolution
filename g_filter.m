function s_int = g_filter(omega,s_int,sigma)
% function applies Gaussian filter to the input spectrum
%
omega_max = max(omega);
g_f = exp(-(omega/(omega_max*sigma)).^2);
Norm = sum(g_f);
g_f  = g_f/Norm;
s_int = s_int.*g_f;

