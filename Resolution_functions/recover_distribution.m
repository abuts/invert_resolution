function recover_distribution(data_file_name)

if ~exist('data_file_name','var')
    data_file_name= 'Pulse_V520588_delta_input_data.mat';
end
load(data_file_name,'tsample','fsample','vsample','V_pulseI','t_det','f_det_vs_t','L_det','L_samp','t_chop','tau_char','V_char');


f_det_counts = build_distribution(t_det,f_det_vs_t,10000);
t_min = min(t_det);
t_max = max(t_det);
%                         InvertPulse3(f_samp,t_samp, v_samp,   t_det,       f_det_vs_t,  L_det,V_pulse, tau_char,V_char,conv_pl_h,varargin)
[f_det_conv,v_det_conv] = InvertPulse3(fsample,tsample,vsample,[t_min,t_max],f_det_counts,L_det,V_pulseI,tau_char,V_char,[],@vel_distribution_delta);

[~,dv_four] = build_bins(v_det_conv);
Norm0  = abs(f_det_conv*dv_four');

%---------------------------------------------------

acolor('b');
p1 = IX_dataset_1d(v_det_conv/V_char,abs(f_det_conv));
p1.x_axis = sprintf('Velocity transfer/(%3.2g m/sec)',V_char);
p1.s_axis = 'probability density ';

dl(p1);
dvs = vsample(2)-vsample(1);
[vel_transf_source,f_d_source] = vel_distribution0(dvs);
acolor('g');
[~,dv_four] = build_bins(vel_transf_source);
NormI = f_d_source*dv_four';
p2 = IX_dataset_1d(vel_transf_source/V_char,f_d_source*(Norm0/NormI));
pl(p2);
acolor('r');
p1.signal = imag(f_det_conv);
p1.s_axis = 'Img error';
pl(p1);


v_transf = convert2v_transf(t_det,L_det,L_samp,t_chop,V_pulseI);
% HACK!!!!

v0 = 0;
%[~,im] = max(f_det_vs_t);
%v0 = v_transf(im);
v_transf = v_transf-v0; % fixing elastic line
[v_transf,ind] = sort(v_transf);
f_det_vs_t_co = f_det_vs_t(ind);
[~,dv_trans] = build_bins(v_transf);

Norm1 = f_det_vs_t_co*dv_trans';
f_det_vs_t_co = f_det_vs_t_co*(Norm0/Norm1);
%pn = IX_dataset_1d(v_transf/V_char,f_det_vs_t./dv);
pn = IX_dataset_1d(v_transf/V_char,f_det_vs_t_co);
acolor('k')
pl(pn);
