%function calc_

[f_samp,t_samp,v_samp,tau_char,V_char,V_pulse,L_samp,t0_chop,norm] = propagate_pulse_to_sample(3);
% [V_char]  = m/sec
% [L_samp] = m;
% [tau_char] = sec
num_pulses = numel(t_samp);
%t0 = [1]; % the equvalent time of moderator*chopper pulse arrival maximum, recalculated to
% moderator position (like in SNS reduction)
L_det = 2.5;
colors = {'r','g','b','k','m'};
f_max_1f = [];
for i=1:num_pulses
    %[f_as,t_as,v_as] = convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [f_as,t_as,v_as] = fft_convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},V_char);
    [xi,yi]=meshgrid(t_as/tau_char,v_as/V_char);
    %
    fn = sprintf('Sample time/velocity profile N %d',i);
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    surf(xi,yi,f_as,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
    
    [f_det,t_det,v_det] = propagate_pulse(f_as,t_as,v_as,L_det);
    %[f_det,t_det,v_det] = fft_propagate_pulse(f_as,t_as,v_as,L_det);
    
    [xi,yi]=meshgrid(t_det/tau_char,v_det/V_char);
    
    fn = sprintf('Detector time/velocity profile N %d',i);
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    surf(xi,yi,abs(f_det),'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
    f_det_vs_t = sum(f_det,1)*norm(i)/size(f_det,1);
    
    
    
    t0 = (L_samp+L_det)/(V_pulse(i));
    %     figure(111);
    %     hold on;
    %     acolor 'r';
    %     plot((t_det-t0)/tau_char,f_det_vs_t);
    %     ax = gca;
    %     ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    %     ax.YLabel.String = sprintf('Signal');
    
    v_transf = convert2v_transf(t_det,L_det,L_samp,t0_chop(i),V_pulse(i));
    % HACK!!!!
    [f_max,im] = max(f_det_vs_t);
    v0 = v_transf(im);
    v_transf = v_transf-v0; % fixing elastic line
    if isempty(f_max_1f)
        f_max_1f = f_max;
    else
        mult = f_max_1f/f_max;
        f_det_vs_t = f_det_vs_t*mult;
    end
    figure(112)
    plot(v_transf/V_char,f_det_vs_t);
    hold on
    ax = gca;
    ax.XLabel.String = sprintf('Velocity transfer/(%3.2g m/s)',V_char);
    ax.YLabel.String = sprintf('Signal');
    
end


