%function calc_

[f_samp,t_samp,v_samp,tau_char,V_char,V_pulse,L_samp] = propagate_pulse_to_sample(1);
% [V_char]  = m/sec
% [L_samp] = m;
% [tau_char] = sec
num_pulses = numel(t_samp);
L_det = 2.5;
for i=1:num_pulses
    %[f_as,t_as,v_as] = convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [f_as,t_as,v_as] = fft_convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [xi,yi]=meshgrid(t_as,v_as);
    figure('Name',sprintf('Sample time/velocity profile N %d',i));
    surf(xi,yi,f_as,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    
    [f_det,t_det,v_det] = propagate_pulse(f_as,t_as,v_as*V_char,L_det,tau_char);

    [xi,yi]=meshgrid(t_det,v_det/V_char);
    figure('Name',sprintf('Detector time/velocity profile N %d',i));
    surf(xi,yi,f_det,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    f_det_vs_t = sum(f_det,1)/size(f_det,1);
    
    
    t0 = (L_samp+L_det)/(V_char*V_pulse(i))/tau_char;
    figure;
    plot(t_det-t0,f_det_vs_t);
    ax = gca;    
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Signal');
    
end