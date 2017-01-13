function mat=get_sq_matr(filename,n_tau,n_lambda)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mat = xlsread(filename);

% wavelength, A
lambda = mat(:,1);
% pulse time (mks)
time   = mat(:,2);
%
signal = mat(:,2);
%
Lambda_min= min(lambda);
Lambda_max= max(lambda);
Tay_min= min(time);
Tay_max= max(time);
dTau = (Tau_max-Tay_min)/(n_tau-1);

lam0 = lambda(1);
block0 = lambda==lam0;
Nblock = sum(block0);







