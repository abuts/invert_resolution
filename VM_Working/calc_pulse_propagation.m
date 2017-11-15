filename = 'mod_e_table.csv';

if ~exist('mat','var')
    mat = xlsread(filename);
end

[in_mat,time_mod,Mod_energy]= get_sq_matr(mat,400);

e_i = [150,65.21,35.79,22.52,15.48];

L = 10.0; % m
L_samp = 1.8; % m

[tau,tau_char,R_chop] = t_vs_e(e_i,L);
L_det = 4.3; % m

f_chop =  cell(1,numel(tau));
t_chop = cell(1,numel(tau));
e_chop = cell(1,numel(tau));
f_samp =  cell(1,numel(tau));
t_samp = cell(1,numel(tau));
e_samp = cell(1,numel(tau));

for i=1:numel(tau)
    [fmod,tchop,vchop] = extract_moderator_f(in_mat,time_mod,Mod_energy,tau(i),L,tau_char,200);
    % convert velocity m/s to mEv
    vSq2mEv = 5.227e-6;  % s^2/m^2
    Echop = vchop.^2*vSq2mEv;
    
    t_chop{i} = tchop;
    e_chop{i} = Echop;
    figure('Name',sprintf('chop pulse for Ei=%3.1f',e_i(i)));
    [xi,yi]=meshgrid(tchop,Echop);
    surf(xi,yi,fmod,'EdgeColor','none');
    v_chop_h  = vchop*tau_char; % Convert velocity into characteristic velocity used by chop pulse
    fchop = chop_pulse(v_chop_h,tchop,tau(i),R_chop);
    surf(xi,yi,fchop,'EdgeColor','none');
    f_chop{i} = fmod.*fchop;
    surf(xi,yi,f_chop{i},'EdgeColor','none');
    
    [fsamp,tsamp,vsamp] = propagate_pulse(f_chop{i},tchop,vchop,L_samp,tau_char);
    [xi,yi]=meshgrid(tsamp,vsamp);
    surf(xi,yi,fsamp,'EdgeColor','none');
    f_samp{i} = fsamp;
    t_samp{i} = t_samp;
    e_samp{i} = vchop.^2*vSq2mEv;
    
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




