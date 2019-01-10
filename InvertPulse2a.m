function [v_distr,vel_steps] =  InvertPulse2a(f_samp,t_samp,v_samp,t_det,f_det_vs_t,L_det,V_pulse,tau_char,V_char,conv_pl_h,vel_distr_fun)
%                                             fsample,tsample,vsample,t_det,f_det_vs_t,L_det,V_pulseI,tau_char,V_char,conv_pl_h,vel_distr_fun
% Ragged Fourier transformation
% The method does not give correct direct problem solution
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

% time scale at the arrival to sample
T_samp_min = min(t_samp); % Zero point, all time should start from this
dT_samp  = max(t_samp)-T_samp_min;
dt_samp = t_samp(2)-t_samp(1);
t_samp  = t_samp-T_samp_min;

if isempty(t_det) % event mode, f_det_vs_t is sequence of events.
    % Check the events starging time!
    f_det_vs_t = f_det_vs_t-T_samp_min;
    event_mode = true;
    td_min = min(f_det_vs_t);
    td_max = max(f_det_vs_t);
else
    if all(size(t_det) ==[1,2])
        event_mode = true;
        td_min = t_det(1)-T_samp_min;
        td_max = t_det(2)-T_samp_min;
        
    else
        event_mode = false;
        t_det = t_det-T_samp_min;
        td_min = min(t_det);
        td_max = max(t_det);
    end
end
% Detector's signal time interval
T_min = td_min;
T_max = td_max;



V_max = L_det/T_min;
V_min = L_det/T_max;
dV = (V_max -V_min);
V_avrg = 0.5*(V_max + V_min);
v2_range = V_min:dv:V_max;

% maximal accessible interval of velocity change due to the sample scattering
%Dv = (dV-dV0)/2;



%Nt =

ti = L_det./(v2_range);
ti = sort(ti);
dt_v = min(ti(2:end)-ti(1:end-1));
if dt_v > dt_samp
    dv = L*dt_samp/(T_min*(T_min-dt_samp));
    v2_range = V_min:dv:V_max;
    ti = L_det./(v2_range);
    ti = sort(ti);
end
Nv = numel(v2_range);
Nt = numel(ti);


[xb,yb] = meshgrid(t_samp+T_min,v_samp);
[xi,yi] = meshgrid(ti,v2_range);

f_sampI = interp2(xb,yb,f_samp,xi,yi,'linear',0);


if event_mode
    [tbin_edges,t_bins] = build_bins(ti+min(t_samp));
    fd = histcounts(f_det_vs_t,tbin_edges);
    f_det_vs_t = fd./t_bins;
else
end

name = vel_distr_fun('name');
cache_file_name = pulse_name(V_pulse,[name,'_2a_resolution_matrix'],Nt,Nv);
if exist([cache_file_name,'.mat'],'file')
    
    load([cache_file_name,'.mat'],'difr_matrix','res_spectrum','omega_v','omega_t');
    if isempty(difr_matrix)
        difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,L_det);
        save(cache_file_name,'difr_matrix','res_spectrum','omega_v','omega_t');
    end
else    
    [omega_t,omega_v,res_spectrum] = sfft2(ti,v2_range,f_sampI);
    difr_matrix = [];
    save(cache_file_name,'difr_matrix','res_spectrum','omega_v','omega_t');
    
    difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,L_det);
    save(cache_file_name,'difr_matrix','res_spectrum','omega_v','omega_t');
    
end
Err = check_difraction_matrix(difr_matrix,v2_range,omega_v,omega_t,L_det);
fprintf(' Total error from the diffraction matrix: (%g,%g)\n',real(Err),imag(Err));

phase_shift = exp(-1i*omega_t*T_min);
res_spectrum = res_spectrum.*phase_shift;

res_matrix = res_spectrum.*difr_matrix;


vel_steps = v2_range-V_avrg;
check_propagation(res_matrix,ti,omega_t,vel_steps,tau_char,conv_pl_h,vel_distr_fun);



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
    [res_spectrum,omega_vr] = p_filter(res_matrix,omega_v,n_harm_left);
    
    Sm = linsolve(conj(res_spectrum'),intensity');
    
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

function [f_t,t_range]=check_propagation(res_matrix,t_range,omega_t,vel_transf,tau_char,conv_pl_h,vel_distr,nt_range,nv_range)

if ~exist('nv_range','var')
    filter = false;
else
    filter = true;    
end

v_min = min(vel_transf);
v_max = max(vel_transf);
V_av = 0.5*(v_min+v_max);
[vel_transf,f_d] = vel_distr(vel_transf-V_av);
[~,sv] = sfft(vel_transf,f_d);
if filter
    svv = zeros(size(sv));
    svv(nv_range) = sv(nv_range);
    svv = conj(svv');
    rm = res_matrix(nt_range,:);
else
    svv = conj(sv');
    rm = res_matrix;
end

fr_nm = sum(rm.*svv,2);
if filter
    f_nm(nt_range) = fr_nm;
else
    f_nm = fr_nm;
end
[t_range,f_t] = isfft(omega_t,f_nm,t_range);
%[t_range,ind] = sort(ti);
%f_t = f_t(ind);

[~,dt]=build_bins(t_range);
Norm =tau_char*sum(real(f_t).*dt');
f_t = f_t/Norm;

pn = IX_dataset_1d(t_range/tau_char,real(f_t));
pn.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
pn.s_axis = 'Signal/Per unit time';
acolor('b');
if ~isempty(conv_pl_h)
    make_current(conv_pl_h);
end

pl(pn);
pimg = IX_dataset_1d(t_range/tau_char,imag(f_t));
acolor('r');
pl(pimg)
fprintf('max real f_t = %f max imag f_t = %f\n',max(real(f_t)),max(imag(f_t)));

function Err = check_difraction_matrix(difr_matrix,v_range,omega_v,omega_t,L_det)

Nt = numel(omega_t);
Nv = numel(omega_v);

v_peak = v_range(end-1);

exps = exp(-1i*omega_v*v_peak);
Err_row = 1i*zeros(1,Nt);
for n=1:Nt
    expi  = exps*exp(1i*omega_t(n)*L_det/v_peak);
    
    
    LH = sum(expi.*difr_matrix(n,:))/Nv;
    difr  = 1-LH;
    Err_row(n) = difr;
    %     difr_ph = atan2(imag(difr),real(difr))*180/pi;
    %     fprintf('n: %d Difr: (%f, %f) mod: %f, phase: %f\n',n,real(difr),imag(difr),abs(difr),difr_ph);
end
Err = sum(Err_row);
