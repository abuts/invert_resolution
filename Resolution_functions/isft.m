function [t,f] = isft(omega,sf,t)
% slow inverse Fourier transformation.
% Inputs:
% omega -- signal frequencies. Should be caclulated by sft.
% sf    -- amplitudes of Fourier harmonics (from sft)
%
% Outputs:
% x -- vector of x-axis of the signal
% f -- vector of function values

if any(size(omega) ~= size(sf))
    error('SFT:invalid_argument',' number of omega values and number of harmonics coefficients have to be equal');
end
if size(omega,1) ~= 1
    omega = omega';
    sf = sf'; % f shoule be a row.
end

Np = numel(omega);
if ~exist('t','var')
    Np2 = floor(Np/2);
    omega_max = max(abs(omega));
    T = 2*pi*Np2/omega_max;
    dt = T/(Np-1);
    t = dt*(0:Np-1);
end

e_w_matrix = exp(1i*omega.*t'); % omega changes along rows.

f = sum(e_w_matrix.*repmat(sf,Np,1),2);
