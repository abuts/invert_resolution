function [out_sp,omega_out] = symmeterize_spectrum(in_sp,omega_in)
% make Fourier spectrum symmetric and conjugeted to get
% real solution from inverse Fourier transformation
nin = numel(in_sp);
nin2 = floor(nin/2);
if rem(nin,2) == 0
    half1 = 2:(nin2+1);
    half2 = nin:-1:(nin2+1); 
    omg_m = max(abs(omega_in));
    omg_h = [omega_in(2:nin2),omg_m];
    omega_out = [omega_in(1),omg_h,-fliplr(omg_h)];
else
    half1 = 2:(nin2+1);
    half2 = nin:-1:(nin2+2);
    
    omega_out = omega_in;    
end


ampl = 0.5*(abs(in_sp(half1))+abs(in_sp(half2)));
phase = 0.5*(in_sp(half1)+in_sp(half2))./ampl;
zer_point = (ampl == 0);
phase(zer_point) = 0;

half_sp = ampl.*phase;
out_sp  = [in_sp(1);half_sp ;conj(half_sp(end:-1:1))];
