function [tau,tau_char,R_chop] = t_vs_e(en,L)
% Get time of flight to the distance L in the units of the chopper
% opening for given energy
%
% L in meters
% energy in mEv;
%
% Chopper parameters:
R_ch = 5; %mm
W_ch = 2*pi*500; % Hz
H_ch = 0.2; % mm
Gamma = W_ch*R_ch*R_ch/H_ch;
tau_char = R_ch/Gamma; % sec

e_transf_const = 5.22725e-6; % sec^2/m^2

%tau1 = L/tau_char*sqrt(e_transf_const./en);
%
tau0 = L*sqrt(e_transf_const./max(en));
t_step = pi/W_ch; %sec;
t_steps =tau0:t_step:tau0+t_step*(numel(en)-1);

tau = t_steps/tau_char;
R_chop = R_ch/1000; % convert in meters
