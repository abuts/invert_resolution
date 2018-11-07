function [omega_t,omega_v,sf] = sfft2(t,v,f)
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

dT = max(t_edges) - min(t_edges);

omega_t = (2*pi/dT)*ind_t;
t_min_s = min(t); % shift from the area used by fft2 to the input area
shift_t = exp(1i*t_min_s*omega_t);
t_bins = (t_bins.*shift_t)'*(Npt/dT);    


[v_edges,v_bins] = build_bins(v);
dV = max(v_edges) - min(v_edges);
omega_v = (2*pi/dV)*ind_v;
v_min_s = min(v);
shift_v = exp(1i*v_min_s *omega_v);
v_bins = v_bins.*shift_v*(Npv/dV);    




dv_dt= t_bins.*v_bins; % omega changes along rows, t -- along columns

sf = fft2(f).*dv_dt;

% sf = zeros(Npt,Npv);
% parfor m=1:Npv
%     fprintf(' step N%d#%d\n',m,Npv);
%     for n=1:Npt                                      % let's try equal step for the time being.
%         omgf = exp(-1i*(omega_t(n)*t+omega_v(m)*v)); %.*dv_dt;
%         sf(n,m) = sum(reshape(f.*omgf,1,Npt*Npv))/(Npt*Npv);
%     end
% end
omega_t = omega_t';
