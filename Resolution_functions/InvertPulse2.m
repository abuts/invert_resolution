function [v_distr,vel_steps] =  InvertPulse2(f_samp,t_samp,v_samp,t_det,f_det_vs_t,L_det,V_pulse,tau_char,V_char,conv_pl_h)

if ~exist('conv_pl_h','var')
    conv_pl_h = [];
end

% velocities of the incoming pulse
v0_min  = min(v_samp);
v0_max  = max(v_samp);
v0_avrg = 0.5*(v0_max+v0_min);
Nv0     = numel(v_samp);
dV0     = v0_max-v0_min;
dv      = dV0/(Nv0-1);

% velocities at the arrival
dT_samp  = max(t_samp)-min(t_samp);
if isempty(t_det) % event mode, f_det_vs_t is sequence of events
    event_mode = true;
    td_min = min(f_det_vs_t);    
    td_max = max(f_det_vs_t);            
else
    if all(size(t_det) ==[1,2])
        event_mode = true;            
        td_min = t_det(1);    
        td_max = t_det(2);       
        
    else
        event_mode = false;    
        td_min = min(t_det);    
        td_max = max(t_det);        
    end
end
T_min = td_min -min(t_samp);
T_max = td_max -min(t_samp)+dT_samp;


V_max = L_det/T_min;
V_min = L_det/T_max;
dV = (V_max -V_min);
V_avrg = 0.5*(V_max + V_min);

% maximal accessible interval of velocity change due to the sample scattering
%Dv = (dV-dV0)/2;


v2_range = V_min+0.5*dv:dv:V_max;

Nv = numel(v2_range);
ind = fft_ind(Nv);
% if  rem(Nv,2) >0
%     Nv2 = floor((Nv-1)/2);
%     ind = -Nv2:Nv2;
% else
%     Nv2 = floor(Nv/2);
%     ind = -Nv2:Nv2-1;
% end



[xb,yb] = meshgrid(t_samp,v_samp);
[xi,yi] = meshgrid(t_samp,v2_range);

f_sampI = interp2(xb,yb,f_samp,xi,yi,'linear',0);


ti = L_det./(v2_range);
ti = sort(ti);
if event_mode
   [tbin_edges,t_bins] = build_bins(ti+min(t_samp));    
   fd = histcounts(f_det_vs_t,tbin_edges);
   f_det_vs_t = fd./t_bins;
else
end

cache_file_name = pulse_name(V_pulse,'resolution_matrix');
if exist([cache_file_name,'.mat'],'file')
    load([cache_file_name,'.mat'],'res_matrix','omega_v','ti');
else
    t_rel = t_samp-min(t_samp);    
    [omega_v,rm_t] = sft(v2_range,f_sampI,ind);
    
    res_matrix = zeros(Nv,Nv);
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
vel_steps = v2_range-V_avrg;

check_propagation(res_matrix,ti+min(t_samp),vel_steps,tau_char,conv_pl_h)


% invert propagation:
if event_mode
    intensity = f_det_vs_t;
else
    intensity =  interp1(t_det-min(t_samp) ,f_det_vs_t,ti,'linear',0);
end

in  = input('Enter number of harmonics to keep or "q" to finish: ','s');

while true
    n_harm_left = textscan(in,'%d');
    n_harm_left = n_harm_left{1};
    fprintf(' processing %d harmonics\n',n_harm_left);
    [rm,omega_vr] = p_filter(res_matrix,omega_v,n_harm_left);
    
    Sm = linsolve(conj(rm'),intensity');
    
    [vel_steps,v_distr] = isft(omega_vr,Sm,vel_steps);
    fn = sprintf('Recoverted velocity transfer distribuion');
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    plot(vel_steps/V_char,real(v_distr),vel_steps/V_char,imag(v_distr))
    in = input('Enter number of harmonics to keep or "q" to finish: ','s');
    if strncmpi(in,'q',1)
        break;
    end
end
%


function check_propagation(res_matrix,ti,vel_transf,tau_char,conv_pl_h)

[vel_transf,f_d] = vel_distribution0(vel_transf);
[omega_dv,sv] = sft(vel_transf,f_d);

f_t = sum(res_matrix.*sv',1);

[tis,ind] = sort(ti);
f_t = f_t(ind);

[~,dt]=build_bins(tis);
% bin_centers = 0.5*(tis(1:end-1)+tis(2:end));
% dt_0 = tis(2)-tis(1);
% dt_e = tis(end)-tis(end-1);
% bin_centers = [bin_centers(1)-dt_0,bin_centers,bin_centers(end)+dt_e];
% dt = bin_centers(2:end)-bin_centers(1:end-1);
Norm =tau_char*(real(f_t)*dt');
f_t = f_t/Norm;

pn = IX_dataset_1d(tis/tau_char,real(f_t));
pn.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
pn.s_axis = 'Signal/Per unit time';
acolor('r');
if ~isempty(conv_pl_h)
    make_current(conv_pl_h);
end

pl(pn);
pimg = IX_dataset_1d(tis/tau_char,imag(f_t));
pl(pimg)


