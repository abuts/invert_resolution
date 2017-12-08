function [f_samp,t_samp,v_samp,tau_char,V_char,V_pulse,L_sample,t0] = propagate_pulse_to_sample(n_en)
% function simulates the shape of tabulated neutronic moderator pulse at
% sample after beeng cut by chopper.
%
% simulation done for MERLIN, the run N22346

filename = 'mod_e_table.csv';
persistent mat;

if isempty(mat)
    mat = xlsread(filename);
end
% reshape input data in the format x,y,e into matrix with axis
% time, energy and matrix, containing moderator pulse as function of these
% variables, stored in square matrix. Need to know the real time axis
% lenght
[in_mat,time_mod,mod_Energy]= get_sq_matr(mat,400);

e_i = [150,65.21,35.79,22.52,15.48];
if exist('n_en','var')
    n_en = min(n_en,numel(e_i));
    e_i = e_i(1:n_en); %
end
En_transf = 437.384; % conversion from root of energy (in mEv) to neutron velocity m/sec
%v_transf = 3.962e+3; % m/s (transform wavelength (A) into velocity)
V_mod = En_transf*sqrt(mod_Energy); % velocity at moderator (m/sec);
%V_mod = v_transf/lambda_mod; % velocity at moderator (m/sec);

L_chop = 10.0; % m
L_samp = 1.8; % m

[tau,tau_char,R_chop] = t_vs_e(e_i,L_chop);
%L_det = 2.5; % m
V_char = R_chop/tau_char; % m/sec

f_chop =  cell(1,numel(tau));
t_chop = cell(1,numel(tau));
%e_chop = cell(1,numel(tau));
f_samp =  cell(1,numel(tau));
t_samp = cell(1,numel(tau));
v_samp = cell(1,numel(tau));
V_pulse = zeros(size(tau));
t0     = zeros(size(tau));
norm   = zeros(size(tau));

fh = findobj('type','figure', 'Name', 'Moderator pulse');
if  isempty(fh)
    figure('Name','Moderator pulse');
else
    figure(fh);
end
%[xi,yi] = my_meshgrid(time_mod/tau_char,V_mod/V_char);
%surf(xi,yi,in_mat,'EdgeColor','none');

[xi,yi] = meshgrid(time_mod/tau_char,V_mod/V_char);
surf(xi,yi,in_mat','EdgeColor','none');

ax = gca;
ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);




for i=1:numel(tau)
    V_pulse(i) = L_chop/tau(i); %
    % calculate time-velocity profile at chopper position for given chopper
    % opening time.
    [fmod,tchop,vchop] = extract_moderator_f(in_mat,time_mod,V_mod,tau(i),L_chop,tau_char,V_char,200);
    % convert velocity m/s to mEv
    %vSq2mEv = 5.227e-6;  % s^2/m^2
    %Echop = vchop.^2*vSq2mEv;
    
    t_chop{i} = tchop;
    %e_chop{i} = Echop;
    fn = sprintf('chop pulse for Ei=%3.1f',e_i(i));
    fh = findobj('type','figure', 'Name', fn);
    if  isempty(fh)
        figure('Name',fn);
    else
        figure(fh);
    end
    
    [xi,yi]=meshgrid(tchop/tau_char,vchop/V_char);
    surf(xi,yi,fmod,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    
    
    %v_chop_h  = vchop*tau_char; % Convert velocity into characteristic velocity used by chop pulse
    % Calculate chopper pulse shape as function of time and neutrons velocity
    % at chopper position
    fchop = chop_pulse(vchop*tau_char,tchop/tau_char,tau(i)/tau_char,R_chop);
    %fchop = chop_pulse(vchop,tchop,tau(i),R_chop);
    surf(xi,yi,fchop,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    
    %
    % velocity profile cut by choper opening:
    f_chop{i} = fmod.*fchop;
    surf(xi,yi,f_chop{i},'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    %
    norm(i) = sum(reshape(f_chop{i},1,numel(f_chop{i})));
    
    fchs =  sum(f_chop{i},1);
    [~,im] = max(fchs);
    t0(i) = tchop(im)-tau(i);
    
    
    % Calculate profile propagated to a sample.
    [fsamp,tsamp,vsamp] = propagate_pulse(f_chop{i},tchop,vchop,L_samp);
    [xi,yi]=meshgrid(tsamp/tau_char,vsamp/V_char);
    surf(xi,yi,fsamp,'EdgeColor','none');
    ax = gca;
    ax.XLabel.String = sprintf('Time/(%3.2g sec)',tau_char);
    ax.YLabel.String = sprintf('Velocity/(%3.2g m/s)',V_char);
    
    f_samp{i} = fsamp;
    t_samp{i} = tsamp;
    %e_samp{i} = vchop.^2*vSq2mEv;
    v_samp{i} = vchop;
    
    %     % time spread at detector due to energy spread and velocity spread
    %     tau_v  = L_det./vchop/tau_char; % time spread due to velocity distribution variations
    %     tau_av = 0.5*(min(tau_v)+max(tau_v));
    %     tau_v  = tau_v-tau_av;
    %     t_chop = t_chop-tau(i);
    %     [xi,yi]=meshgrid(t_chop,tau_v );
    %     surf(xi,yi,fmod{i},'EdgeColor','none');
    %
end
%surf(min_mat);
L_sample = L_samp+L_chop;
max_pulse = max(norm);
norm = norm/max_pulse;
for i=1:numel(tau)
    f_samp{i} =  f_samp{i}/norm(i);
    t_samp{i}  = t_samp{i};
    v_samp{i}  = v_samp{i};
end




