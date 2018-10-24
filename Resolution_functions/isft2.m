function [t,v,f] = isft2(omega_t,omega_v,sf,t,v)
% slow inverse 2D Fourier transformation conserving phase
% Inputs:
% omega_t -- time frequencies. Should be caclulated by sft2.
% omega_v -- velocity frequencies. Should be caclulated by sft2.
% sf    -- amplitudes of Fourier harmonics (from sft2) where time harmonics
%          are arranged along columns
%
% Outputs:
% t -- vector of x-axis of the signal
% v -- vector of y-axis of the signal
% f -- 2D matrix of function values

if size(omega_t,2) ~= 1 % let's arrange omega_t into columns.
    omega_t = omega_t';
end
if exist('t','var') && size(t,2) ~=1
    t = t';
end

if size(omega_v,1) ~= 1 % let's arrange omega_v into rows.
    omega_v = omega_v';
end
if exist('v','var') && size(v,1) ~=1
    v = v';
end

Nt = size(omega_t,1);
Nv = size(omega_v,2);

if size(sf,1)~= Nt || size(sf,2)~= Nv
    error('ISFT2:invalid_argument',' inconsistent number of omega axis and harmonics matrix');
end


if ~exist('t','var') || any(size(t) ~=size(omega_t))
    if ~exist('t','var')
        t_min = 0;
    else
        t_min = min(t);
    end
    t = build_inverse_axis(omega_t,t_min);
end
if ~exist('v','var') || any(size(v) ~=size(omega_v))
    if ~exist('v','var')
        v_min = 0;
    else
        v_min = min(v);
    end
    v = build_inverse_axis(omega_v,v_min);
end


f = zeros(Nt,Nv);
for j=1:Nv
    for i=1:Nt
        ex_mat = exp(1i*(omega_t*t(i)+omega_v*v(j)));
        f(i,j) = sum(reshape(sf.*ex_mat,1,Nv*Nt));
    end
end


function  t = build_inverse_axis(omega,t_min)
Np = numel(omega);
Np2 = floor(Np/2);
omega_max = max(abs(omega));
T = 2*pi*Np2/omega_max;
dt = T/Np;
t = t_min+dt*(0:Np-1);
