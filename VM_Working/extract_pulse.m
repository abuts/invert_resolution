filename = 'mod_e_table.csv';

if ~exist('mat','var')
    mat = xlsread(filename);
end

[in_mat,time_mod,Mod_energy]= get_sq_matr(mat,400);

e_i = [150,65.21,35.79,22.52,15.48];

L = 10.0; % m

[tau,tau_char,R_chop] = t_vs_e(e_i,L);
L_det = 4.3; % m

fmod =  cell(1,numel(tau));
t_mod = cell(1,numel(tau));
e_mod = cell(1,numel(tau));
for i=1:numel(tau)
    [f_mod,t_chop,v_chop] = extract_moderator_f(in_mat,time_mod,Mod_energy,tau(i),L,tau_char,200);
    % convert velocity m/s to mEv
    vSq2mEv = 5.227e-6;  % s^2/m^2
    E_chop = v_chop.^2*vSq2mEv;
    tau_v  = L_det./v_chop/tau_char; % time spread due to velocity distribution variations
    
    t_mod{i} = t_chop;
    e_mod{i} = E_chop;
    figure('Name',sprintf('chop pulse for Ei=%3.1f',e_i(i)));
    [xi,yi]=meshgrid(t_chop,E_chop);
    surf(xi,yi,f_mod,'EdgeColor','none');
    v_chop  = v_chop*tau_char; % ?onvert velocity into characteristic velocity used by chop pulse
    f_chop = chop_pulse(v_chop,t_chop,tau(i),R_chop);
    surf(xi,yi,f_chop,'EdgeColor','none');
    fmod{i} = f_mod.*f_chop;
    surf(xi,yi,fmod{i},'EdgeColor','none');
    
    % time spread at detector due to energy spread and velocity spread
    tau_av = 0.5*(min(tau_v)+max(tau_v));
    tau_v  = tau_v-tau_av;
    t_chop = t_chop-tau(i);    
    [xi,yi]=meshgrid(t_chop,tau_v );
    surf(xi,yi,fmod{i},'EdgeColor','none');
    
end

%surf(min_mat);




