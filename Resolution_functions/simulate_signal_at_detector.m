function simulate_signal_at_detector()
persistent conv_pl_h;

[f_samp,t_samp,v_samp,tau_char,V_char,V_pulse,L_samp,t0_chop,norm] = propagate_pulse_to_sample(3);% [V_char]  = m/sec
% [L_samp] = m;
% [tau_char] = sec
num_pulses = numel(t_samp);
%t0 = [1]; % the equvalent time of moderator*chopper pulse arrival maximum, recalculated to
% moderator position (like in SNS reduction)
L_det = 2.5;
colors = {'k','r','g','b','m'};
f_max_1f = [];
recovered_dirst_h = [];
vel_distr_fun = @vel_distribution0;
for i=1:num_pulses
    %[f_as,t_as,v_as] = convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},tau_char,V_char);
    [f_afs,t_afs,v_afs,Norm0] = fft_convolute_with_vel_distr(f_samp{i},t_samp{i},v_samp{i},V_char,true,vel_distr_fun);
    [xi,yi]=meshgrid(t_afs/tau_char,v_afs/V_char);
    %
    fn = sprintf('After-Sample time/velocity profile N %d',i);
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    surf(xi,yi,abs(f_afs),'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
    
    t_range = max(t_afs)-min(t_afs);
    [f_det,t_det,v_det] = propagate_pulse(f_afs,t_afs,v_afs,L_det,t_range);
    %[f_det,t_det,v_det] = fftv_propagate_pulse(f_as,t_as,v_as,L_det);
    
    [xi,yi]=meshgrid((t_det-min(t_afs))/tau_char,v_det/V_char);
    
    fn = sprintf('Detector time/velocity profile N %d',i);
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    surf(xi,yi,abs(f_det),'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('(Det time -min(t_{samp}))/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
    %---------------------------------------------------
    % signal at detector
    
    f_det_vs_t = sum(f_det,1);
    Norm0 = sum(f_det_vs_t)*tau_char*(t_det(2)-t_det(1));
    f_det_vs_t = f_det_vs_t/Norm0;
    
    
    %Norm1 = sum(sum(f_afs,1))*tau_char*(t_afs(2)-t_afs(1));
    f_samp_vs_t = sum(f_afs,1)/Norm0;
    
    
    %
    t_sampl_min = min(t_samp{i});
    pn = IX_dataset_1d((t_det-t_sampl_min)/tau_char,f_det_vs_t);
    pn.x_axis = sprintf('Time -min(t_{samp})/(%3.2g sec)',tau_char);
    pn.s_axis = 'Signal/per unit time';
    acolor(colors(i));
    if ~isempty(conv_pl_h)
        make_current(conv_pl_h);
    end
    acolor('k');
    conv_pl_h=pl(pn);
%     acolor('g');
%     pn = IX_dataset_1d((t_afs-t_sampl_min)/tau_char,f_samp_vs_t);
%     pl(pn);
    %
    %---------------------------------------------------
    ds = data_saver(vel_distr_fun,t_samp{i},v_samp{i},f_samp{i},...
        V_pulse(i),t0_chop(i),...
        t_det,f_det_vs_t,L_det,L_samp,tau_char,V_char);
    [f_det_dec,v_det_dec] = InvertPulse3(ds,conv_pl_h);
    stop;
    [~,dv_four] = build_bins(v_det_dec);
    Norm0  = abs(f_det_dec*dv_four');
    
    %---------------------------------------------------
    
    acolor('b');
    p1 = IX_dataset_1d(v_det_dec/V_char,abs(f_det_dec));
    p1.x_axis = sprintf('Velocity transfer/(%3.2g m/sec)',V_char);
    p1.s_axis = 'probability density ';
    if isempty(recovered_dirst_h)
        recovered_dirst_h= dl(p1);
        [vel_transf_source,f_d_source] = vel_distr_fun(v_det_dec);
        acolor('g');
        [~,dv_four] = build_bins(vel_transf_source);
        NormI = f_d_source*dv_four';
        p2 = IX_dataset_1d(vel_transf_source/V_char,f_d_source*(Norm0/NormI));
        pl(p2);
    else
        make_current(recovered_dirst_h);
        recovered_dirst_h= pl(p1);
    end
    acolor('r');
    p1.signal = imag(f_det_dec);
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
    keep_figure;
    
end

