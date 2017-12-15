function [f_out,t_out,v_out] = fftv_propagate_pulse(f_in,time_in,vel_in,L)
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
DT0 = max(time_in)-min(time_in);
DV0 = v_max - v_min;
dt = time_in(2)-time_in(1);
%Nt = size(f_in,2);
%dt = (tp_max-tp_min)/(Nt-1);
% if dt<2e-6
%     dt = 2e-6;
% end
t_out = tp_min:dt:tp_max;
Nt = numel(t_out);

v_out = vel_in;
dV = vel_in(2)-vel_in(1);
dN = v_min/dV;
Nv = numel(vel_in);


Del = dN/Nv;
DelV = L/(DT0*DV0);
j = (0:Nv-1)/Nv;
Imkl = @(m,k,l)(sum(exp(2i*pi*((m-k)*(Del+j)-l*DelV./(Del+j)))))/Nv;
%I1 = Imkl(1,1,1);
%



% fft frequencies do not correspond to indexes as index m in fft array 
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
% indexes of frequncies corresponging to V-frequencies in f_in_sp;
Nv2 = floor(Nv/2);
if  rem(Nv,2) >0
    nW = [0:Nv2,-Nv2:1:-1];
else
    nW = [0:Nv2-1,-Nv2:1:-1];
end
Nt2 = floor(Nt/2);
if  rem(Nt,2) >0
    nT = [0:Nt2,-Nt2:1:-1];
else
    nT = [0:Nt2-1,-Nt2:1:-1];
end


t_start=tic;
f_out_sp = zeros(Nv,Nt);
parfor k=1:Nv
    fprintf('k=%d#%d\n',k,Nv);
    for l=1:Nt
        Im = zeros(1,numel(nW));        
        for mj=1:Nv
            Im(mj) = Imkl(nW(mj),nW(k),nT(l));
        end
        f_out_sp(k,l) = sum(f_in_sp(:,l).*Im');
    end
end
time_c = toc(t_start)/60 % convert in minutes
f_out = ifft2(f_out_sp);



