function [omega,sf,e_w_matrix] = sft(t,f,ind)
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
    f = f'; % t change in f shoule be located in a column.
end

if (size(t,1) ~= size(f,1))
    error('SFT:invalid_argument',' number of x values and function values have to be equal');
end

Np = size(t,1);
Nmat = size(f,2);
if ~exist('ind','var')
    if  rem(Np,2) >0
        Np2 = floor((Np-1)/2);
        ind = [0:Np2,-Np2:1:-1];
    else
        Np2 = floor(Np/2);
        ind = [0:Np2,-(Np2-1):1:-1];
    end
end
% only equal step works
dt = max(t(2:end)-t(1:end-1));
%
omega = (2*pi/(Np*dt))*ind; % Important! period is wider than t_max - t_min
e_w_matrix = exp(-1i*omega.*t); % omega changes along rows.

sf = zeros(Nmat,Np);
for i=1:Nmat
    sf(i,:) = sum(e_w_matrix.*repmat(f(:,i),1,Np),1)/(Np);
end
%omega = omega';