function [v_distr,vel_steps] =  InvertPulse3(in_data,conv_pl_h)

[t_samp,v_samp,f_samp,  ...
    V_pulse,tau_chop,...
    t_det,f_det_vs_t,...
    L_det,L_samp,tau_char,V_char,...
    vel_distr] = in_data.get_data();
rebin_to_velocity = in_data.use_velocity_scale;
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
T_samp_min = min(t_samp); % Zero point, all time should start from this
dT_samp  = max(t_samp)-T_samp_min;
dt_samp = t_samp(2)-t_samp(1);
t_samp  = t_samp-T_samp_min;

[td_min,td_max,event_mode,t_det,f_det_vs_t] = select_det_scale(t_det,f_det_vs_t,T_samp_min);
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


v2_range = V_min:dv:V_max;

Nv = numel(v2_range);
if rebin_to_velocity
    dt_samp = L_det/V_max - L_det/(V_max+dv);
end
t_range = T_min:dt_samp:T_max;
Nt = numel(t_range);

%ti = L_det./(v2_range);
%ti = sort(ti);
if event_mode
    [tbin_edges,t_bins] = build_bins(t_range);
    fd = histcounts(f_det_vs_t-T_samp_min,tbin_edges);
    f_det_vs_t = fd./t_bins;
end
fprintf('Defining diffraction equation with [%d,%d] elements\n',Nt,Nv);

% name = vel_distr('name');
% cache_file_name = pulse_name(V_pulse,[name,'_FFTresolution_matrix'],Nt,Nv);
% if exist([cache_file_name,'.mat'],'file')
%     load([cache_file_name,'.mat'],'rm','difr_matrix','omega_v','omega_t');
%     if difr_matrix == 0
%         difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,L_det);
%         save(cache_file_name,'difr_matrix','rm','omega_v','omega_t');
%     end
% else
% original signal at sample projected to detector position.
[tb,vb] = meshgrid(t_samp+T_min,v_samp); % time is shifted -- phase shift ignored ?
[tpi,vpi] = meshgrid(t_range,v2_range);
f_samp_extended = interp2(tb,vb,f_samp,tpi,vpi,'linear',0);

[omega_t,omega_v,rm] = sfft2(t_range,v2_range,f_samp_extended');
%    difr_matrix = 0;
%    save(cache_file_name,'difr_matrix','rm','omega_v','omega_t');

difr_matrix=calc_difr_matrix(omega_v,omega_t,v2_range,L_det);
%    save(cache_file_name,'difr_matrix','rm','omega_v','omega_t');
%end
%[difr_matrix,Err]=calc_difr_matrix(omega_v,omega_t,v2_range,L_det);
%fprintf(' Total error of the diffraction matrix: %g\n',Err);
Err = check_difraction_matrix(difr_matrix,v2_range,omega_v,omega_t,L_det);
fprintf(' Total error from the diffraction matrix: (%g,%g)\n',real(Err),imag(Err));

phase_shift = exp(-1i*omega_t*T_min);
rm = rm.*phase_shift;
%
%[n_max,m_max]  = find_max_ind(rm,1.e-6);
res_matrix = rm.*difr_matrix;
vel_steps = v2_range;


[fte,t_steps]=check_propagation(res_matrix,t_range,omega_t,vel_steps,tau_char,conv_pl_h,vel_distr);
%check_propagation(res_matrix,t_range,omega_t,vel_steps,tau_char,conv_pl_h,vel_distr);


% invert propagation:
if event_mode
    intensity = f_det_vs_t;
else
    %intensity = real(fte);
    intensity =  interp1(t_steps,real(fte),t_range,'linear',0);
    %intensity =  interp1(t_det,f_det_vs_t,t_range,'linear',0);
    in_data.ds.f_det_vs_t = interp1(t_steps,real(fte),t_det,'linear',0);
    in_data.save_data();
    if  ~isempty(conv_pl_h)
        make_current(conv_pl_h);
        [~,dt]=build_bins(t_range);
        Norm =tau_char*(real(intensity)*dt');
        intensity_v = intensity/Norm;
        
        pn = IX_dataset_1d(t_range/tau_char,real(intensity_v));
        acolor('g');
        pl(pn);
    end
    %     pulse_data_file_name = pulse_name(V_pulse,[name,'_input_data']);
    %     load(pulse_data_file_name,'tsample','fsample','vsample','V_pulseI','t_det','f_det_vs_t','L_det','L_samp','t_chop','tau_char','V_char');
    %     t_det = t_range;
    %     f_det_vs_t = intensity_v;
    %     save(pulse_data_file_name,'tsample','fsample','vsample','V_pulseI','t_det','f_det_vs_t','L_det','L_samp','t_chop','tau_char','V_char');
    
end
[omega_t,s_int] = sft(t_range,intensity);

sv = svds(res_matrix,1);
%in  = input('Enter accuracy and sigma avrg if requested or "q" to finish: ','s');
in = '1000 1.e-4';
eps = 1.e-4;
while true
    in_val = textscan(in,'%f %f');
    sigma = in_val{1};
    if ~isempty(in_val{2})
        eps = in_val{2};
    end
    if eps < 0
        break;
    end
    fprintf('To select conditionality %f and keeping singular values larger than: %f and applying Gaussian filder with sigma=%f \n ',eps,eps*sv,sigma);
    int_r = g_filter(omega_t,s_int,sigma);
    Sm = lsqminnorm(res_matrix, conj(int_r'), eps*sv);
    %     neglect = abs(S)<=1.e-4*max(abs(diag(S)));
    %     S(neglect) = 0;
    %     Sm = s_int*U*S*V;
    %[rm,int_r,omega_vt] = p_filter3(res_matrix,s_int,omega_v,omega_t,n_harm_left);
    %Sm = pinv(res_matrix,1.e-6)*conj(int_r');% linsolve(res_matrix,conj(int_r'));
    %Sm = linsolve(rm,conj(int_r'));
    [Sm,omega_vt] = symmeterize_spectrum(Sm,omega_v);
    
    
    [vel_steps,v_distr] = isft(omega_vt,Sm,min(v2_range)-V_avrg);
    fn = sprintf('Recoverted velocity transfer distribuion');
    fh = findobj('type','figure', 'Name', fn);
    
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    plot(vel_steps/V_char,abs(v_distr),vel_steps/V_char,imag(v_distr))
    in  = input('Enter accuracy and sigma avrg if requested or "q" to finish: ','s');
    if strncmpi(in,'q',1)
        break;
    end
end
close(fh);
%


function [f_t,t_range]=check_propagation(res_matrix,t_range,omega_t,vel_transf,tau_char,conv_pl_h,vel_distr)

v_min = min(vel_transf);
v_max = max(vel_transf);
V_av = 0.5*(v_min+v_max);
[vel_transf,f_d] = vel_distr(vel_transf-V_av);
[omega_dv,sv] = sft(vel_transf,f_d);

f_nm = sum(res_matrix.*sv,2);

[t_range,f_t] = isft(omega_t,f_nm,t_range);
%[t_range,ind] = sort(ti);
%f_t = f_t(ind);

[~,dt]=build_bins(t_range);
Norm =tau_char*(real(f_t)*dt');
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
%--------------------------------------------------------------------------
function difr_matrix=calc_difr_matrix(omega_v,omega_t,v_range,L_det)
%Nv = numel(omega_v);
%Nt = numel(omega_t);
ti = L_det./v_range;

% if nargout>1
%     calc_error = true;
% else
%     calc_error = false;
% end
% 
% v_peak = v_range(20);
% dv_j = (v_range-v_peak);
% exp2 = exp(1i*omega_v.*dv_j');
% SM = sum(exp2,2)/Nv;
% ST = exp(1i*omega_t(1)*(L_det/v_peak-ti'));
% Int = sum(SM.*ST);


v_mat = exp(1i*(omega_v.*v_range'));
t_mat = exp(-1i*(omega_t.*ti));
difr_matrix = t_mat*v_mat;
% if calc_error
%     test_row = zeros(1,Nt);
%     for n=1:Nt
%         exp2 = exp(-1i*(omega_v*v_peak-omega_t(n)*L_det/v_peak));
%         test_row(n) = sum(exp2.*difr_matrix(n,:))/Nv;
%     end
%     Err = sum(test_row)/Nt-1;
% end
%difr_matrix = zeros(Nt,Nv);
% for n=1:Nt
%     t_vec = t_mat(n,:);
%     for m=1:Nv
%         %difr_matrix(n,m) =sum(exp(1i*(omega_v(m)*v_range-omega_t(n)*ti)));
%         difr_matrix(n,m) =sum(v_mat(m,:).*t_vec);
%     end
%     if calc_error
%         exp2 = exp(-1i*(omega_v*v_peak-omega_t(n)*L_det/v_peak));
%         test_row(n) = sum(exp2.*difr_matrix(n,:))/Nv;
%     end
% end
%--------------------------------------------------------------------------
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

function [td_min,td_max,event_mode,t_det,f_det_vs_t] = select_det_scale(t_det,f_det_vs_t,T_samp_min)
if isempty(t_det) % event mode, f_det_vs_t is sequence of events.
    % Check the events starting time!
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
