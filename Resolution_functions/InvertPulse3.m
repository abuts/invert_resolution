function [v_distr,vel_steps] =  InvertPulse3(f_samp,t_samp,v_samp,t_det,f_det_vs_t,L_det,V_pulse,tau_char,V_char,conv_pl_h,varargin)

if ~exist('conv_pl_h','var')
    conv_pl_h = [];
end
if nargin > 5
    vel_distr = varargin{1};
else
    vel_distr = @vel_distribution0;
end


% velocities of the incoming pulse
v0_min  = min(v_samp);
v0_max  = max(v_samp);
v0_avrg = 0.5*(v0_max+v0_min);
Nv0     = numel(v_samp);
dV0     = v0_max-v0_min;
dv      = dV0/(Nv0-1);

% velocities at the arrival
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
        td_min = t_det(1);
        td_max = t_det(2);
        
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
V_min = L_det/T_max; %
dV = (V_max -V_min);
if dV <= 0
    error('INVERT_PULSE3:invalid_argument',' wrong detector"s time range (smaller than sample time range)')
end
V_avrg = 0.5*(V_max + V_min);

% maximal accessible interval of velocity change due to the sample scattering
%Dv = (dV-dV0)/2;


v2_range = V_min+0.5*dv:dv:V_max;

Nv = numel(v2_range);


ti = L_det./(v2_range);
ti = sort(ti);
if event_mode
    [tbin_edges,t_bins] = build_bins(ti);
    fd = histcounts(f_det_vs_t,tbin_edges);
    f_det_vs_t = fd./t_bins;
end

cache_file_name = pulse_name(V_pulse,'resolution_delta_matrix');
if exist([cache_file_name,'.mat'],'file')
    load([cache_file_name,'.mat'],'rm','difr_matrix','omega_v','omega_t','ti');
    t_range = T_min:dt_samp:T_max;    
    if difr_matrix == 0
        difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,ti);
        save(cache_file_name,'difr_matrix','rm','omega_v','omega_t','ti');        
    end
else
    [tb,vb] = meshgrid(t_samp+T_min,v_samp); % time is shifted -- phase shift ignored ?
    t_range = T_min:dt_samp:T_max;
    [tpi,vpi] = meshgrid(t_range,v2_range);
    f_samp_extended = interp2(tb,vb,f_samp,tpi,vpi,'linear',0);
    
    [omega_t,omega_v,rm] = sft2(t_range,v2_range,f_samp_extended');
    difr_matrix = 0;
    save(cache_file_name,'difr_matrix','rm','omega_v','omega_t','ti');
    
    difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,ti);
    save(cache_file_name,'difr_matrix','rm','omega_v','omega_t','ti');
end
res_matrix = rm.*difr_matrix;
vel_steps = v2_range-V_avrg;

check_propagation(res_matrix,t_range,omega_t,vel_steps,tau_char,conv_pl_h,vel_distr)


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


function check_propagation(res_matrix,t_range,omega_t,vel_transf,tau_char,conv_pl_h,vel_distr)

[vel_transf,f_d] = vel_distr(vel_transf);
[omega_dv,sv] = sft(vel_transf,f_d);

f_nm = sum(res_matrix.*sv,2);

[t_range,f_t] = isft(omega_t,f_nm,t_range);



[~,dt]=build_bins(t_range);
Norm =tau_char*(real(f_t)*dt');
f_t = f_t/Norm;

pn = IX_dataset_1d(t_range/tau_char,real(f_t));
pn.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
pn.s_axis = 'Signal/Per unit time';
acolor('r');
if ~isempty(conv_pl_h)
    make_current(conv_pl_h);
end

pl(pn);
pimg = IX_dataset_1d(t_range/tau_char,imag(f_t));
pl(pimg)

function difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,ti)
Nv = numel(omega_v);
Nt = numel(omega_t);
difr_matrix = zeros(Nt,Nv);
for m=1:Nv
    for n=1:Nt
        difr_matrix(n,m) =sum(exp(1i*(omega_v(m)*v2_range-omega_t(n)*ti)));
    end
end
