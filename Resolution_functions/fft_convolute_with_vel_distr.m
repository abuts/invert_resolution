function [f_s_out,t_samp,v_s_out,Norma,ft_sample,ft_vt] = fft_convolute_with_vel_distr(f_samp,t_samp,v_samp,V_char,no_plot,varargin)
% Calculate time-velocity distribution for beam after propagating though sample given
% time/velicity distribution of incident beam and model sample velocity
% gain/loss probability.
%
%
%   Detailed explanation goes here
if exist('no_plot','var')
    noplot = no_plot;
else
    noplot = false;
end
if nargin > 5
    vel_distr = varargin{1};
else
    vel_distr = @vel_distribution0;
end



dvs = v_samp(2)-v_samp(1);
[~,~,~,dV_scat_max] = vel_distr(dvs);

%dV_scat_max = max(vel_transf);



%v_samp_norm = v_samp/V_char;
V_min_fin = min(v_samp)-dV_scat_max;
if V_min_fin<0
    V_min_fin = 0;
end
V_max_fin = max(v_samp)+dV_scat_max;



dV = dvs;
nSource_pts = floor((max(v_samp)-min(v_samp))/dV);
if nSource_pts < 16
    dV = (max(v_samp)-min(v_samp))/16;
    %nSource_pts = 16;
end
vi =V_min_fin:dV:V_max_fin;
[~,v_bins] = build_bins(vi);
[~,t_bins] =  build_bins(t_samp);
[xi,yi] = meshgrid(t_samp,vi);
[xb,yb] = meshgrid(t_samp,v_samp);

f_samp = interp2(xb,yb,f_samp,xi,yi,'linear',0);


Nv = numel(vi);
Nt = size(f_samp,2);

bin_mat = v_bins'.*t_bins;
Norma = sum(reshape(f_samp.*bin_mat,1,Nt*Nv));
f_samp = f_samp/Norma;

v_sh = 0.5*(V_min_fin+V_max_fin);
% caclculate velocity transfer distrubution in its own range but with
% joing accuracy.
% extend velocity transfer function onto the full scattering scale
[vel_transf,f_d] = vel_distr(vi-v_sh);



Norma = f_d*v_bins';
f_d = f_d/Norma;

if ~noplot
    fh = findobj('type','figure', 'Name', 'sample distribution');
    if  isempty(fh)
        figure('Name','sample distribution');
    else
        figure(fh);
    end
    surf(xi,yi/V_char,f_samp,'EdgeColor','none');
    ax = gca;
    %ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    view(0,90);
end
ft_sample = fft2(f_samp);
ind = fft_ind(Nv);
phase = exp(1i*pi*ind); % correct for symmetric integration range only

%ftm = abs(fftshift(ft));
%surf(ftm,'EdgeColor','none');
ft_vt = fft(f_d);

if ~noplot
    figure(111)
    acolor 'b';
    plot(vel_transf/V_char,f_d);
    ax = gca;
    ax.XLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    figure(112);
    ind = fft_ind(numel(ft_vt));
    freq = fftshift((1/(dV*Nv))*ind);
    plot(freq,fftshift(abs(ft_vt)));
end

ft_vt = ft_vt.*phase;
%fm = bsxfun(@times,ft_sample,ft_vt');
fm = ft_sample.*ft_vt';

%fm = ft.*repmat(ffv',1,numel(t_samp));
%fm = fftshift(fm);
%ftm = abs(fftshift(fm));
%surf(ftm,'EdgeColor','none');


%fm = repmat(fm(1,:),size(fm,1),1);
f_s_out = ifft2(fm);

if ~noplot
    fh = findobj('type','figure', 'Name', 'convoluted distribution');
    if  isempty(fh)
        figure('Name','convoluted distribution');
    else
        figure(fh);
    end
    %[xi,yi]= meshgrid();
    surf(xi,yi/V_char,abs(f_s_out),'EdgeColor','none');
    view(0,90);
    
end
% ax = gca;
% %ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
% ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
if max(max(imag(f_s_out)))>1.e-9
    warning('substantial imaginary part of the propagated distribution')
end
f_s_out = abs(f_s_out);
v_s_out = vi;





function fc = f1_kernel(dv,v_samp,f_samp,v_s_out,V_char)

Mf = @(v)(interp1(v_samp,f_samp,v,'linear',0));

[~,fc]= vel_distribution(dv,V_char);
fc = fc.*Mf(v_s_out-dv);



