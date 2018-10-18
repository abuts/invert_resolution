function [v_distr,vel_steps] =  InvertPulse2(f_samp,t_samp,v_samp,t_det,f_det_vs_t,L_det,V_pulse,tau_char,V_char,conv_pl_h)

% velocities of the incoming pulse
v0_min  = min(v_samp);
v0_max  = max(v_samp);
v0_avrg = 0.5*(v0_max+v0_min);
Nv0     = numel(v_samp);
dV0     = v0_max-v0_min;
dv      = dV0/(Nv0-1);

% velocities at the arrival
dT_samp  = max(t_samp)-min(t_samp);
T_min = min(t_det)-min(t_samp);
T_max = max(t_det)-min(t_samp)+dT_samp;

V_max = L_det/T_min;
V_min = L_det/T_max;
dV = (V_max -V_min);
V_avrg = 0.5*(V_max + V_min);

% maximal accessible interval of velocity change due to the sample scattering
Dv = (dV-dV0)/2;


v2_range = V_min+0.5*dv:dv:V_max;
Nv = numel(v2_range);
ind = fft_ind(Nv);

[xb,yb] = meshgrid(t_samp,v_samp);
[xi,yi] = meshgrid(t_samp,v2_range);

f_sampI = interp2(xb,yb,f_samp,xi,yi,'linear',0);


ti = L_det./(v2_range);

cache_file_name = pulse_name(V_pulse,'resolution_matrix');
if exist([cache_file_name,'.mat'],'file')
    load([cache_file_name,'.mat'],'res_matrix','omega_v','ti');
else
    
    [omega_v,rm_t] = sft(v2_range,f_sampI,ind);
    
    res_matrix = zeros(Nv,Nv);
    t_rel = t_samp-min(t_samp);
    for i=1:numel(ti)
        t_start = ti(i);
        ti_shift   = t_start-ti;
        for m=1:Nv
            rm_ij      = interp1(t_rel,rm_t(:,m),ti_shift,'linear',0);
            res_matrix(m,i) = rm_ij*exp(1i*omega_v(m)*v2_range')/Nv;
        end
    end
    save(cache_file_name,'res_matrix','omega_v','ti');
end
check_propagation(res_matrix,ti+min(t_samp),v2_range-V_avrg,tau_char,conv_pl_h)



function check_propagation(res_matrix,ti,vel_transf,tau_char,conv_pl_h)

[vel_transf,f_d] = vel_distribution0(vel_transf);
[omega_dv,sv] = sft(vel_transf,f_d);

f_t = sum(res_matrix.*sv',1);
[tis,ind] = sort(ti);
f_t = f_t(ind);
pn = IX_dataset_1d(tis/tau_char,real(f_t));
pn.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
pn.s_axis = 'Signal';
acolor('r');
if ~isempty(conv_pl_h)
    make_current(conv_pl_h);
end

pl(pn);
pimg = IX_dataset_1d(tis/tau_char,imag(f_t));
pl(pimg)


