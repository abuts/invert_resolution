function [t,f] = isft(omega,sf,t)
% slow inverse Fourier transformation.
% Inputs:
% omega -- signal frequencies. Should be caclulated by sft.
% sf    -- amplitudes of Fourier harmonics (from sft)
%
% Outputs:
% x -- vector of x-axis of the signal
% f -- vector of function values

if size(omega,1) ~= 1 % let's arrange omege into rows.
    omega = omega';
end
Np = numel(omega);
if ~any(size(sf)==Np)
    error('ISFT:invalid_argument',' number of omega values and number of harmonic coefficients have to be equal');
end

if size(sf,1) == Np %let's arrange sf in rows
    sf = conj(sf');
end
Nspec = size(sf,1);

if ~exist('t','var') || any(size(t) ~=size(omega))
    if ~exist('t','var') 
        t_min = 0;
    else
        t_min = min(t);        
    end
    Np2 = floor(Np/2);
    omega_max = max(abs(omega));
    T = 2*pi*Np2/omega_max;
    dt = T/Np;
    t = t_min+dt*(0:Np-1);
else
    if size(t,1)~=1
        t=t'; % let's arrange t into rows
    end
    if any(size(t) ~= size(omega))
        error('ISFT:invalid_argument',...
            ' sizes of frequency and time arrays (if one is provided) have to be equal');
    end
    
end

e_w_matrix = exp(1i*omega.*t'); % omega changes along rows.

f = zeros(Np,Nspec);
for i=1:Nspec
    f(:,i) = sum(e_w_matrix.*repmat(sf(i,:),Np,1),2);
end
f = f';
