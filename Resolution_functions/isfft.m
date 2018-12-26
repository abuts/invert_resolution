function [t,f] = isfft(omega,sf,t)
% slow inverse Fourier transformation.
% Inputs:
% omega -- signal frequencies. Should be caclulated by sft.
% sf    -- amplitudes of Fourier harmonics (from sft)
%
% Outputs:
% x -- vector of x-axis of the signal
% f -- vector of function values

if size(omega,2) ~= 1 % let's arrange omega into columns.
    omega = omega';
end
Np = numel(omega);
if ~any(size(sf)==Np)
    error('ISFFT:invalid_argument',' number of omega values and number of harmonic coefficients have to be equal');
end

if size(sf,2) == Np %let's arrange sf in columns;
    sf = conj(sf');
end
Nspec = size(sf,2);
if ~exist('t','var')
    t=build_t(omega,0);
    t_bins = ones(1,Np);
else
    if size(t,2) ~= 1
        t=t'; % let's arrange t into columns
    end
    if any(size(t,1) ~= size(omega,1))
        t=build_t(omega,min(t));
    end
    [t_edges,t_bins] = build_bins(t);
    dT = max(t_edges)-min(t_edges);
    omega_1 = omega(2);
    T1 = 2*pi/omega_1;
    if(abs(dT-T1)>1.e-6)
        error('ISFFT:invalid_argument',...
            ' inconsistent time and omega axis');
    end
    t_min_s = min(t);
    shift = exp(1i*t_min_s*omega');
    t_bins = (t_bins.*shift)'*(Np/dT);
end

f = ifft(sf./t_bins);

function t=build_t(omega,t_min)
Np = numel(omega);
omega_1 = omega(2);
T = 2*pi/omega_1;
dt = T/Np;
t = (t_min+dt*(0:Np-1))';

