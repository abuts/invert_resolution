function [omega_t,omega_v,sf] = sft2(t,v,f,dT,dV)
% slow 2D Fourier transformation conserving phase
% Inputs:
% t -- vector of x-axis of the bin centers of the signal
% v -- vector of y-axus of the bin centers of the signal
% f -- 2D matrix of function values
% Outputs:
% omega_t -- signal frequencies along axis x;
% omega_v -- signal frequencies along axis y;
% sf      -- matxis of the Fourier harmonics amplitudes.
%
if size(t,2) ~= 1
    t = t';
end
if size(v,1) ~= 1
    v = v';
end

if size(f,1) ~= size(t,1) || size(f,2) ~= size(v,2)
    error('SFT2:invalid_argument','size of input matrix and input axis are inconsistent');
end
Npt = size(t,1);
Npv = size(v,2);
ind_t = fft_ind(Npt);
ind_v = fft_ind(Npv);
[t_edges,t_bins] = build_bins(t);
if ~exist('dT','var') || isempty(dT)
    dT = max(t_edges) - min(t_edges);
end
omega_t = (2*pi/dT)*ind_t;

[v_edges,v_bins] = build_bins(v);
if ~exist('dV','var') || isempty(dV)
    dV = max(v_edges) - min(v_edges);
end
omega_v = (2*pi/dV)*ind_v;

dv_dt= t_bins'.*v_bins; % omega changes along rows, t -- along columns

sf = zeros(Npt,Npv);
parfor m=1:Npv
    fprintf(' step N%d#%d\n',m,Npv);
    for n=1:Npt       
        omgf = exp(-1i*(omega_t(n)*t+omega_v(m)*v)).*dv_dt;
        sf(n,m) = sum(reshape(f.*omgf,1,Npt*Npv))/(dV*dT);
    end
end
omega_t = omega_t';
