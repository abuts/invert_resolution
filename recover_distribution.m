function recover_distribution(data_file_name)

if ~exist('data_file_name','var')
    data_file_name= 'Pulse_V520588_5peaks_input_data.mat';
end
in_data = data_saver();
in_data = in_data.load_data(data_file_name);
t_det = in_data.ds.t_det;

f_det_vs_t  = in_data.ds.f_det_vs_t;


f_det_counts = build_distribution(t_det,f_det_vs_t,100000000);
t_min = min(t_det);
t_max = max(t_det);
hh=histogram(f_det_counts,100);

in_data.ds.t_det = [t_min,t_max];
in_data.ds.f_det_vs_t = f_det_counts;

[f_det_conv,v_det_conv] = InvertPulse3(in_data,hh);

[~,dv_four] = build_bins(v_det_conv);
Norm0  = abs(f_det_conv*dv_four');
V_char = in_data.ds.V_char;
%---------------------------------------------------

acolor('b');
p1 = IX_dataset_1d(v_det_conv/V_char,abs(f_det_conv));
p1.x_axis = sprintf('Velocity transfer/(%3.2g m/sec)',V_char);
p1.s_axis = 'probability density ';

dl(p1);
vsample = in_data.ds.vsample;
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

L_det  = in_data.ds.L_det;
L_samp = in_data.ds.L_samp;
t_chop = in_data.ds.t_chop;
V_pulseI=in_data.ds.V_pulseI;
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
