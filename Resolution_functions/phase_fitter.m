function fint = phase_fitter(integrals,v_steps,nt,v_min,v_max)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nv = numel(v_steps);
betta = v_min*v_max*nt;

k=1:Nv;
CM  = zeros(Nv);
for i=1:Nv
    phase_l = exp(2i*pi*(betta/v_steps(i)-k*v_steps(i)));
    CM(i,:) = phase_l.*integrals';
end
delta = CM\eye(Nv);
fint = integrals.*delta;
