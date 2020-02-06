function propagate_deltaf
persistent conv_pl_h;

V_char = 1;
tau_char = 1;
V_pulse = 50; % delta function at given V;
t_afs = 0:0.1:100;
dv = 0.1;
v_afs = 1:dv:100;
t_ind = 50;
ind = floor(((V_pulse - min(v_afs)))/dv)+1;
f_afs = zeros(numel(v_afs),numel(t_afs));
f_afs(ind,t_ind) = 1;


L_det = 50;
t_range = max(t_afs)-min(t_afs);
[f_det,t_det,v_det] = propagate_pulse(f_afs,t_afs,v_afs,L_det,t_range);
ft = sum(f_det,1);
[omg1,fts] = sfft(t_det,ft);
[tr,ftr] = isfft(omg1,fts,t_det);
plot(t_det,ft);
hold on
plot(tr,ftr);



f_det_vs_t = sum(f_det,1);
Norm0 = sum(f_det_vs_t)*tau_char*(t_det(2)-t_det(1));
f_det_vs_t = f_det_vs_t/Norm0;

fn = sprintf('Detector time/velocity profile N 1');
fh = findobj('type','figure', 'Name', fn);
if  isempty(fh)
    figure('Name',fn);
else
    figure(fh);
end
[xi,yi]=meshgrid((t_det-min(t_afs))/tau_char,v_det/V_char);
surf(xi,yi,abs(f_det),'EdgeColor','none');
ax = gca;
ax.XLabel.String = sprintf('(Det time -min(t_{samp}))/(%3.2g sec)',tau_char);
ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
view(0,90);
t_sampl_min = min(t_afs);
pn = IX_dataset_1d((t_det-t_sampl_min)/tau_char,f_det_vs_t);
pn.x_axis = sprintf('Time -min(t_{samp})/(%3.2g sec)',tau_char);
pn.s_axis = 'Signal/per unit time';
acolor('k');
if ~isempty(conv_pl_h)
    make_current(conv_pl_h);
end
conv_pl_h=pl(pn);

L_samp=10;
t0_chop = 1;

vel_distr_fun = @delta_v;
in_data = data_saver(vel_distr_fun,t_afs,v_afs,f_afs,...
    V_pulse,t0_chop,...
    t_det,f_det_vs_t,L_det,L_samp,tau_char,V_char);
in_data.use_velocity_scale=false;
[f_det_dec,v_det_dec] = InvertPulse3(in_data,conv_pl_h);
end

function [vel,fr] = delta_v(vel)
delta_v = vel(2:end)-vel(1:end-1);
del_v = max(abs(delta_v));
nz = vel>-0.5*del_v & vel < 0.5*del_v;
fr = zeros(size(vel));
fr(nz) = 1/del_v;
end
