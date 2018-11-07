function [t,v,f] = isfft2(omega_t,omega_v,sf,t,v)
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


if exist('t','var')
    t_min_s = min(t);
    shift_t = exp(-1i*t_min_s*omega_t);
    
    if any(size(t) ~=size(omega_t))
        t = build_inverse_axis(omega_t,t_min_s);
    end
    [t_edges,t_bins] = build_bins(t);
    
    dT = max(t_edges)-min(t_edges);
    t_bins = (t_bins'.*shift_t)*(Nt/dT);    
    
else
    t = build_inverse_axis(omega_t,0);    
    t_bins =ones(Nt,1);
end
if exist('v','var')
    v_min_s = min(v);    
    shift_v = exp(-1i*v_min_s *omega_v);    
    if size(v) ~=size(omega_v)
        t = build_inverse_axis(omega_v,v_min_s);
    end
    [v_edges,v_bins] = build_bins(v);
    dV = max(v_edges)-min(v_edges);
    v_bins = (v_bins.*shift_v)*(Nv/dV);    
else
    t = build_inverse_axis(omega_t,0);    
    v_bins = ones(1,Nv);
end

bins = t_bins.*v_bins;
sf = sf./bins;

f = ifft2(sf);

% f = zeros(Nt,Nv);
% for j=1:Nv
%     fprintf(' step: N%d#%d\n',j,Nv);    
%     for i=1:Nt
%         ex_mat = exp(1i*(omega_t*t(i)+omega_v*v(j)));
%         f(i,j) = sum(reshape(sf.*ex_mat,1,Nv*Nt));
%     end
% end
% 

function  t = build_inverse_axis(omega,t_min)
Np = numel(omega);
omega_1 = omega(2);
T = 2*pi/omega_1;
dt = T/Np;
t = t_min+dt*(0:Np-1);
