function [omega,sf] = sft(t,f,ind)
% slow Fourier transformation.
% Inputs:
% x -- vector of x-axis of the signal
% f -- vector of function values
% Outputs:
% omega -- signal frequencies
% sf    -- amplitudes of Fourier harmonics
%
if any(size(t) ~= size(f))
    error('SFT:invalid_argument',' number of x values and function values have to be equal');
end
if size(t,1) ~= 1
    t = t';
    f = f'; % f shoule be row.
end

Np = numel(t);
if ~exist('ind','var')
    Np2 = floor(Np/2);
    if  rem(Np,2) >0
        ind = -Np2:Np2;
    else
        ind = -Np2:Np2-1;
    end
end
t_min = t(1);
t_max = t(end);
dT = t_max-t_min;
omega = (2*pi/dT)*ind;
e_w_matrix = exp(-1i*omega.*t'); % omega changes along rows.

sf = sum(e_w_matrix.*repmat(f',1,Np),1)/(Np-1);
%omega = omega;