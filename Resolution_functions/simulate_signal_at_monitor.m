%function calc_

[f_samp,t_samp,v_samp,tau_char,V_char] = propagate_pulse_to_sample(5);

num_pulses = numel(t_samp);
for i=1:num_pulses
    [f_as,t_as,v_as] = convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [xi,yi]=meshgrid(t_as,v_as);
    figure('Name',sprintf('Sample time/velocity profile N %d',i));
    surf(xi,yi,f_as,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
end