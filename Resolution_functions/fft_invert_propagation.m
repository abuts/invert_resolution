function [f_out,dv_out] = fft_invert_propagation(f_samp,t_samp,v_samp,f_det_vs_t,t_det,L_det)
% retrieve samle gain/loss velocity distribution given signal on the detector and 
% and the smple incident beam time/velocity distribution. 

fv_samp = sum(f_samp,2);
ft_samp  = sum(f_samp,1);
fv_norm = sum(fv_samp);
ft_norm = sum(ft_samp);
v_av = sum(fv_samp.*v_samp')/fv_norm;
ts_av = sum(ft_samp.*t_samp)/ft_norm;

[xi,yi] = meshgrid(t_det,v_samp);
[xb,yb] = meshgrid(t_samp,v_samp);
f_samp_in = interp2(xb,yb,f_samp,xi,yi,'nearest',0);

t_det = t_det - ts_av;
vel_in = L_det./(t_det);
v_min  = min(vel_in);
v_max  = max(vel_in);
tp_min = min(t_det);
tp_max = max(t_det);

DV0 = v_max - v_min;
dt = t_det(2)-t_det(1);
%Nt = size(f_in,2);
%dt = (tp_max-tp_min)/(Nt-1);
% if dt<2e-6
%     dt = 2e-6;
% end
Nt= numel(t_samp);
Nv= floor(0.5*Nt);


v_min_norm = v_min/DV0;
v_max_norm = v_max/DV0;

v_index = fft_ind(Nv);
t_index = fft_ind(Nt);

Imk = @(n,k)I_mk(n,k,v_min_norm,v_max_norm,v_index,t_index,Nv,Nt);

I1 = Imk(1,2);
%


% fft frequencies do not correspond to indexes as index m in fft array 
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft(f_det_vs_t,Nt);
% indexes of frequncies corresponging to V-frequencies in f_in_sp;

t_start=tic;
f_out_sp = zeros(Nv,Nt);
for k=1:Nv
    fprintf('k=%d#%d\n',k,Nv);
    for l=1:Nt
        f_out_sp(k,l) = Imk(k,l);
    end
end
f_out_sp = f_in_sp(k,n).*f_out_sp;

time_c = toc(t_start)/60 % convert in minutes
f_out = ifft2(f_out_sp);

function int = I_mk(n,k,v_min,v_max,v_index,t_index,Nv,Nt)
persistent cash;
if isempty(cash)
    cash = ones(Nv,Nt)*NaN;
end
int = cash(n,k);
if isnan(int)
    ki = v_index(k);
    ni = t_index(n);    
    fun = @(v)(exp(2i*pi*(ki*v-v_min*ni./v)));
    int = integral(fun,v_min,v_max);
    cash(n,k) = int;
end



