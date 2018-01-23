function [f_out,t_out,v_max] = fft_propagate_pulse_IntWC(f_samp,t_samp,v_samp,L,V_pulse,t_char,v_char)
% Calculate interpolited time-velocity profile at the position L using fft
% W - version integrates from t_min = min(time_in)+L/v_max
%                        to   t_max = max(time_in)+L/v_min;
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

[f_in,time_in,vel_in,Norma,ft_sample,ft_vt] = fft_convolute_with_vel_distr(f_samp,t_samp,v_samp,v_char,'noplot');



v_min = min(vel_in);
v_max = max(vel_in);
%DV0 = v_max - v_min;
% detector time range:
tp_min = min(time_in)+L/v_max;
tp_max = max(time_in)+L/v_min;
DT0 = tp_max-tp_min;
dt = time_in(2) - time_in(1);
if dt<2e-6 % debugging -- not needed in reality
    dt = 2e-6;
end
% detector time scale
[t_out,dt,Nt] = adjust_step(tp_min,tp_max,dt);
% the sample time to reflect onto the detector time
t_orig_min = min(time_in);
t_orig_max = tp_max-L/v_max;
dt0 = (t_orig_max-t_orig_min)/(Nt-1);
[t_in_expanded,dt0,Nt] = adjust_step(t_orig_min,t_orig_max,dt0);

dv = (vel_in(2)-vel_in(1));
%dv = 3*(vel_in(2)-vel_in(1));
dv  = (v_max-v_min)/(Nt-1);
[v_out,dv,Nv] = adjust_step(v_min,v_max,dv);

[xb,yb] = meshgrid(time_in,vel_in);
[xi,yi]= meshgrid(t_in_expanded,v_out);
f_in = interp2(xb,yb,f_in,xi,yi,'linear',0);

v_index = fft_ind(Nv);
t_index = fft_ind(Nt);
DV0 = max(v_out)-min(v_out);

%
%
% fft frequencies do not correspond to indexes as index m in fft array
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = (fft2(f_in,Nv,Nt));


% DNv = Nv-numel(vel_in);
%v_phase =  exp(1i*pi*v_index*((v_min_norm+v_max_norm))); % appears to shift zeros added at the end to zeros, added to the beginning

%
f_in_t  = ifft2(f_in_sp);
figure(20)
[xfi,yfi] = meshgrid(t_out/t_char,v_out/v_char);
surf(xfi,yfi,abs(f_in_t),'Edgecolor','None');
view(0,90);

v_steps_norm = v_out/DV0;
v_min_norm = v_min/(DV0);
v_max_norm = v_max/(DV0);
betta0 = (L/DT0)/DV0;

%fun = @(m,n)I_mk4(m,n,betta0,v_min_norm,v_max_norm,v_index,t_index);
f_con_mat  = build_dif_mat(betta0,v_min_norm,v_max_norm,'IntW',Nv,Nt);

% for n=1:Nt
%     nt = t_index(n);
%     if nt ~= 0
%         f_con_mat(:,n) = phase_fitter(f_con_mat(:,n),v_steps_norm,betta0,nt,v_min_norm,v_max_norm);
%     end
% end
v_phase =  exp(1i*pi*v_index);
f_con_mat = bsxfun(@times,f_con_mat,v_phase');

f_in_sp = f_in_sp.*f_con_mat;

% % %f_in_sp = bsxfun(@times,f_in_sp,v_phase');
% % %v_phase =  exp(1i*pi*v_index);


%t_phase = exp(1i*pi*t_index*(L/(v_min*DT0))); %
t_phase = exp(1i*pi*t_index); %
f_in_sp = bsxfun(@times,f_in_sp,t_phase );
%

f_t = ifft2(f_in_sp);
[xfi,yfi] = meshgrid(t_out/t_char,v_out/v_char);
surf(xfi,yfi,abs(f_t),'Edgecolor','None');
view(0,90);

f_out_sp = sum(f_in_sp,1);
%t_phase =  exp(1i*pi*t_index*(1+2*v_min_norm-(tp_min+tp_max)/(DT0))); %
%t_phase =  exp(-1i*pi*t_index*((tp_min+tp_max)/(DT0))); %
% t_phase =  exp(-2i*pi*t_index*v_min_norm); %
%
f_out = ifft(f_out_sp);
%--------------------------------------------------------------------------
% Testing inverse problem
[xb,yb]=meshgrid(t_samp,v_samp);
f_in_expanded  =  interp2(xb,yb,f_samp,xi,yi,'linear',0);
ft_f_in_exp = ifft2(f_in_expanded);

f_in_sp_test = ft_f_in_exp.*f_con_mat;

if numel(ft_vt) ~=size(f_in_sp_test,1)
    [vel_transf,f_d] = vel_distribution0(dv);
    
    if numel(ft_vt) ~=size(f_in_sp_test,1)
        % extend velocity transfer function onto the full scattering scale
        v_sh = 0.5*(v_min+v_max);
        f_d = interp1(vel_transf,f_d,v_out-v_sh,'linear',0);
    end
    
    ft_vt = fft(f_d);
    ind = fft_ind(numel(ft_vt));
    phase = exp(1i*pi*ind); % correct for symmetric integration range only
    ft_vt = ft_vt.*phase;
end


ft_signal_test = bsxfun(@times,f_in_sp_test,ft_vt');
ft_signal_test = bsxfun(@times,ft_signal_test,t_phase );
signal_test2D  = ifft2(ft_signal_test);


figure(22);
surf(xfi,yfi,abs(signal_test2D),'Edgecolor','None');
view(0,90)


ft_stignal1D = sum(ft_signal_test,1);
signal_test1D = ifft(ft_stignal1D);
figure(23)
hold off
plot(t_out/t_char,signal_test1D);
hold on
plot(t_out/t_char,f_out);
hold off
% v_phase =  exp(1i*pi*v_index);
% f_in_sp = bsxfun(@times,f_in_sp,v_phase');


%f_det_sp = fft(f_det,Nt);
%
t_phase = exp(-1i*pi*t_index); %
ft_stignal1D = ft_stignal1D.*t_phase ;
np= 155;
while isnumeric(np) &&np>0
    fprintf(' filtering using %d harmonics\n',np);
    [f_vel_sp,vel_ind] = filter_and_test(f_in_sp_test,ft_stignal1D,np);
    f_vel_sp= f_vel_sp.*exp(1i*pi*vel_ind)' ;
    vel_dist = ifft(f_vel_sp);
    
    %max_h_num = max(abs(vel_ind));
    dV_period = DV0/2;
    dv_pres = 2*dV_period/(numel(f_vel_sp)-1);
    dv_steps = -dV_period:dv_pres:dV_period;
    figure(25);
    plot(dv_steps/v_char,real(vel_dist));
    figure(26);
    plot(dv_steps/v_char,imag(vel_dist));
    
    np = input('Enter number of harmonics to use in the filter or 0 to exit:>> ');
end

persistent stor;
if isempty(stor)
    stor = struct();
end
pulse_name = ['pls',num2str(round(V_pulse))];
stor.(pulse_name) = [vel_ind',f_vel_sp,dv_steps',vel_dist];

fnms = fieldnames(stor);
[dv,f_d] = vel_distribution0(dv);
f_d = f_d/sum(f_d);
p_con= IX_dataset_1d(dv/v_char,f_d);
p_con.x_axis = sprintf('Velocity Transfer/(%3.2g sec)',v_char);
p_con.s_axis = 'probability';
acolor('b');
%dl(p_con)
colors = {'k','r','g','b','m'};
for i=1:numel(fnms)
    dat = stor.(fnms{i});
    p_con.x = real(dat(:,3)/v_char);
    sig = dat(:,4)/sum(real(dat(:,4)));
    p_con.signal = real(sig);
    p_con.error = imag(sig);
    
    acolor(colors{i});
    pl(p_con);
end
keep_figure;


function [vel_spectra,Nv_left] = filter_and_test(dirf_matrix,signal_spectra,Nv_left,Nt_left)

%dirf_matrix = fftshift(dirf_matrix);
%signal_spectra = fftshift(signal_spectra);
[Nv,Nt] = size(dirf_matrix);
if ~exist('Nt_left','var')
    Nt_left = Nt;
end

if Nv_left>Nv
    Nv_left = Nv;
end
Nv_left = floor(Nv_left/2);
Nt_left = floor(Nt_left/2);
dirf_matrix = dirf_matrix([1:Nv_left+1,Nv-Nv_left:Nv],[1:Nt_left+1,Nt-Nt_left:Nt]);
signal_spectra = signal_spectra([1:Nt_left+1,Nt-Nt_left:Nt]);


vel_spectra  = linsolve((dirf_matrix)',(signal_spectra)');
Nv_ind = fft_ind(Nv);
Nv_left = Nv_ind([1:Nv_left+1,Nv-Nv_left:Nv]);
Nv_left = fftshift(Nv_left);
vel_left = fftshift(vel_spectra);

figure(24);
plot(Nv_left,abs(vel_left));
Nv_left = fftshift(Nv_left);


