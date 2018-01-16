%function calc_

[f_samp,t_samp,v_samp,tau_char,V_char,V_pulse,L_samp,t0_chop,norm] = propagate_pulse_to_sample(1);% [V_char]  = m/sec
% [L_samp] = m;
% [tau_char] = sec
num_pulses = numel(t_samp);
%t0 = [1]; % the equvalent time of moderator*chopper pulse arrival maximum, recalculated to
% moderator position (like in SNS reduction)
L_det = 2.5;
colors = {'k','r','g','b','m'};
f_max_1f = [];
for i=1:num_pulses
    %[f_as,t_as,v_as] = convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [f_afs,t_afs,v_afs,Norm] = fft_convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},V_char);
    [xi,yi]=meshgrid(t_afs/tau_char,v_afs/V_char);
    %
    fn = sprintf('Sample time/velocity profile N %d',i);
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    surf(xi,yi,f_afs,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
    
    [f_det,t_det,v_det] = propagate_pulse(f_afs,t_afs,v_afs,L_det);
    %[f_det,t_det,v_det] = fftv_propagate_pulse(f_as,t_as,v_as,L_det);
    
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
    %---------------------------------------------------
    % signal at detector
    f_det_vs_t = sum(f_det,1)*norm(i)/size(f_det,1);
    %
    pn = IX_dataset_1d(t_det/tau_char,f_det_vs_t);
    pn.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
    pn.s_axis = 'Signal';
    acolor(colors(i));
    dl(pn);
    keep_figure
    [f_det_conv,t_det_conv] = fft_propagate_pulse_IntC(f_afs,t_afs,v_afs,L_det,V_pulse(i),tau_char,V_char);
    %[f_det_conv,t_det_conv] = propagate_pulse_Int(f_afs,t_afs,v_afs,L_det,V_pulse(i),tau_char,V_char);
    p1 = IX_dataset_1d(t_det_conv/tau_char,abs(f_det_conv));
    p1.x_axis = sprintf('Time/(%3.2g sec)',tau_char);
    p1.s_axis = 'Signal';
    dl(p1);
    keep_figure
    p1.signal = imag(f_det_conv);
    p1.s_axis = 'Img error';
    dl(p1);
 
    
    
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
    dv = v_transf(1:end-1)-v_transf(2:end);
    dv = [dv,dv(end)];
    v0 = 0;
    %[~,im] = max(f_det_vs_t);
    %v0 = v_transf(im);
    v_transf = v_transf-v0; % fixing elastic line
    Norm1 = sum(f_det_vs_t.*dv);
    mult = Norm/Norm1;
    f_det_vs_t = f_det_vs_t*mult;
    [v_transf,ind] = sort(v_transf);
    f_det_vs_t = f_det_vs_t(ind);
    
    %pn = IX_dataset_1d(v_transf/V_char,f_det_vs_t./dv);
    pn = IX_dataset_1d(v_transf/V_char,f_det_vs_t);
    acolor(colors(i));
    if i==1
        [reduced_fh,ax]=dl(pn);
    else
        [reduced_fh,ax]=pl(pn,'name',reduced_fh.Name);
    end
    
    %ax = gca;
    ax.XLabel.String = sprintf('Velocity transfer/(%3.2g m/s)',V_char);
    ax.YLabel.String = sprintf('Signal');
    [f_out,v_out] = fft_invert_propagation(f_samp{i},t_samp{i},v_samp{i},f_det_vs_t,t_det,L_det);
    p_con = IX_dataset_1d(v_out/V_char,abs(f_out));
    p_con.x_axis = sprintf('Velocity Transfer/(%3.2g sec)',V_char);
    p_con.s_axis = 'probability';
    dl(p_con);
    keep_figure
    p_con.signal = imag(f_out);
    p_con.s_axis = 'Img error';
    dl(p1);
    
end

