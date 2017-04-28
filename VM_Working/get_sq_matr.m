function [sig_mat,time_int,En_out]=get_sq_matr(mat,npoints)
% return Summary of this function goes here
%
%   Detailed explanation goes here


% incident energy (mEv) %wavelength, A
energy = mat(:,1);
% pulse time (mks)
time   = mat(:,2);
%
signal = mat(:,3);
error  = mat(:,4);
%
%e_min= min(energy );
%e_max= max(energy );
Tau_min= min(time);
Tau_max= max(time);
dTau = (Tau_max-Tau_min)/(npoints-1);
time_int = Tau_min:dTau:Tau_max;

E0 = energy(1);
block0 = energy== E0;
block_size = sum(block0);
n_blocks = numel(time)/block_size;

time  = reshape(time,block_size,n_blocks);
energy = reshape(energy,block_size,n_blocks);
signal = reshape(signal,block_size,n_blocks);
error = reshape(error,block_size,n_blocks);

En_out = energy(1,:);

sig_mat = zeros(npoints,n_blocks);
for i=1:n_blocks
    tf = time(:,i);
    sf = signal(:,i);
    %ef = error(:,i);
    
    sig_mat(:,i) = interp1(tf,sf,time_int,'linear');
    t_min = min(tf);
    t_max = max(tf);
    out = time_int< t_min | time_int > t_max;
    sig_mat(out,i) = 0;
    
    %lambda_i = lambda(1,i);
    %[tau1,tau2,t0,R,a0] = evaluate_tau(tf,sf,lambda_i);
    %
    %
    %pin = [a0,R,tau1,tau2,t0];
    %[yfit,fitpar]=fit(tf,sf,ef,@ikeda,pin,[1,0,1,1,0]);
    %pout = fitpar.p;
    %in =  time_int>t0;
    %sig_mat(in,i) = ikeda(time_int(in),pout);
end



