function [f_s_out,t_samp,v_s_out,Norma] = fft_convolute_with_vel_distr(f_samp,t_samp,v_samp,V_char)
% Calculate time-velocity distribution for beam after propagating though sample given
% time/velicity distribution of incident beam and model sample velocity
% gain/loss probability.
%
%
%   Detailed explanation goes here



Nvp= 1024;
e_max = 25;

e_transf_const = 5.22725e-6; % sec^2/m^2
dV_scat_max = sqrt(e_max/e_transf_const);

%v_samp_norm = v_samp/V_char;
V_min_fin = min(v_samp)-dV_scat_max;
if V_min_fin<0
    V_min_fin = 0;
end
V_max = max(v_samp)+dV_scat_max;

%dV = (V_max-V_min)/(Nvp-1);
dV = 2*V_max/(2*Nvp-1);
nSource_pts = floor((max(v_samp)-min(v_samp))/dV);
if nSource_pts < 16
    dV = (max(v_samp)-min(v_samp))/16;
    %nSource_pts = 16;
end
vi = -V_max:dV:V_max;
[xi,yi] = meshgrid(t_samp,vi);
[xb,yb] = meshgrid(t_samp,v_samp);

f_samp = interp2(xb,yb,f_samp,xi,yi,'nearest',0);


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


Nv = numel(vi);
Nt = size(f_samp,2);

ft = fft2(f_samp);
ind = fft_ind(Nv);
phase = exp(1i*pi*ind); % correct for symmetric integration range only

%ftm = abs(fftshift(ft));
%surf(ftm,'EdgeColor','none');


% caclculate velocity transfer distrubution in its own range but with
% joing accuracy.
[v_t,f_d] = vel_distribution(vi);
Norma = sum(f_d)*(max(v_t)-min(v_t))/(Nv-1);
f_d = f_d/Norma;
figure(111)
acolor 'b';
plot(v_t/V_char,f_d);
ax = gca;
ax.XLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
ffv = fft(f_d);
ffv = ffv.*phase;



fm = bsxfun(@times,ft,ffv');
%fm = ft.*repmat(ffv',1,numel(t_samp));
%fm = fftshift(fm);
%ftm = abs(fftshift(fm));
%surf(ftm,'EdgeColor','none');


%fm = repmat(fm(1,:),size(fm,1),1);
f_s_out = ifft2(fm);

fh = findobj('type','figure', 'Name', 'convoluted distribution');
if  isempty(fh)
    figure('Name','convoluted distribution');
else
    figure(fh);
end
surf(xi,yi/V_char,f_s_out,'EdgeColor','none');
view(0,90);
% ax = gca;
% %ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
% ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);


v_s_out = V_min_fin:dV:V_max;
[xi,yi] = meshgrid(t_samp,v_s_out);
[xb,yb] = meshgrid(t_samp,vi);
f_s_out = interp2(xb,yb,f_s_out,xi,yi,'cubic',0);







function fc = f1_kernel(dv,v_samp,f_samp,v_s_out,V_char)

Mf = @(v)(interp1(v_samp,f_samp,v,'linear',0));

[~,fc]= vel_distribution(dv,V_char);
fc = fc.*Mf(v_s_out-dv);



