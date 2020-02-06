function [omega,sf] = sfft(t,f)
% slow Fourier transformation conserving phase
% Inputs:
% x -- vector of x-axis of the signal
% f -- vector of function values or matxix 
% Outputs:
% omega -- signal frequencies
% sf    -- amplitudes of Fourier harmonics
%
if size(t,2) ~= 1
    t = t';
end
if size(f,1) ~= numel(t)
    f = f'; %  change in f shoule be located in a column to work with fft
end

if (size(t,1) ~= size(f,1))
    error('SFT:invalid_argument',' number of x values and function values have to be equal');
end

Np = size(t,1);

[t_edges,t_bins] = build_bins(t);
dT = max(t_edges)-min(t_edges);
dt = dT/Np;
ind = fft_ind(Np);
omega = (2*pi/dT)*ind; % Important! period is wider than t_max - t_min by one step

t_min_s = min(t);
shift = exp(-1i*t_min_s*omega);

t_bins = conj(t_bins.*shift)'/dt;
%
sf = fft(f).*t_bins;

