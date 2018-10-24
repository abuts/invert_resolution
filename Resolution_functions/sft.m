function [omega,sf] = sft(t,f,ind)
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
    f = f'; % t change in f shoule be located in a column.
end

if (size(t,1) ~= size(f,1))
    error('SFT:invalid_argument',' number of x values and function values have to be equal');
end

Np = size(t,1);
Nmat = size(f,2);
if ~exist('ind','var')
    ind = fft_ind(Np);
end

[t_edges,t_bins] = build_bins(t);
dT = max(t_edges)-min(t_edges);
%
omega = (2*pi/dT)*ind; % Important! period is wider than t_max - t_min
e_w_matrix = exp(-1i*omega.*t).*repmat(t_bins',1,Np); % omega changes along rows, t-along columns.

sf = zeros(Nmat,Np);
for i=1:Nmat
    sf(i,:) = sum(e_w_matrix.*repmat(f(:,i),1,Np),1)/dT;
end
%omega = omega';