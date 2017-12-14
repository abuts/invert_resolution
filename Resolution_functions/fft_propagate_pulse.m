function [f_out,t_out,v_out] = fft_propagate_pulse(f_in,time_in,vel_in,L)
% Calculate interpolited time-velocity profile at the position L using fft
%
% f_in     -- 2D signal function in units tau(mks) vs
% time_in  -- time axis for signal  (sec)
% vel_in   -- velocity accis for signal (m/sec)
% L        -- distance to the target
% tau -- time at choper to shift function around (in chopper opening time
%
% Output:
%  Interpolated fime-velocity profile at position L
% f_out  --
% t_out -- time axis for the profile above (sec)
% v_out -- velocity axis for the profile above (in m/s)
v_min = min(vel_in);
v_max = max(vel_in);
tp_min = L/v_max+min(time_in);
tp_max = L/v_min+max(time_in);
dt = time_in(2) - time_in(1);
if dt<2e-6
    dt = 2e-6;
end
t_out = tp_min-dt:dt:tp_max;
Nt = numel(t_out);
dN = 1;
v_out = vel_in;

Del = dN/Nt;
j=0:(Nt-dN);
bin_bnd = Del/(1-Del)./(1+j/Nt);
dV = bin_bnd(1:end-1)-bin_bnd(2:end);
j = (0:Nt-dN-1)/Nt;
Imkl = @(m,k,l)(sum(dV.*exp(2i*pi*((m-k)./((1-Del)*(Del+j))-l*(Del+j)))));
%I1 = Imkl(1,1,1);
%



Nv = size(f_in,1);
%Nt = size(f_in,2);
% fft frequencies do not correspond to indexes as index m in fft array 
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
% indexes of frequncies corresponging to V-frequencies in f_in_sp;
Nv2 = floor(Nv/2);
nW = [0:Nv2-1,Nv2:-1:1];

t_start=tic;
f_out_sp = zeros(Nv,Nt);
parfor k=1:Nv
    fprintf('k=%d#%d\n',k,Nv);
    for l=1:Nt
        Im = zeros(1,numel(nW));        
        for mj=1:Nv
            Im(mj) = Imkl(nW(mj),nW(k),nW(l));
        end
        f_out_sp(k,l) = sum(f_in_sp(:,l).*Im');
    end
end
time_c = toc(t_start)/60 % convert in minutes
f_out = ifft2(f_out_sp);



