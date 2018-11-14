function [difr_matrix,Err]=calc_difr_matrix(omega_v,omega_t,v_range,L_det)
% Calculates regularized diffraction matrix, binding time and velocity
% scales
Nv = numel(omega_v);
Nt = numel(omega_t);
ti = L_det./v_range;

if nargout>1
    calc_error = true;
else
    calc_error = false;
end

v_peak = v_range(20);
% dv_j = (v_range-v_peak);
% exp2 = exp(1i*omega_v.*dv_j');
% SM = sum(exp2,2)/Nv;
% ST = exp(1i*omega_t(1)*(L_det/v_peak-ti'));
% Int = sum(SM.*ST);
if calc_error
    test_row = zeros(1,Nt);
end
difr_matrix = zeros(Nt,Nv);
for n=1:Nt
    for m=1:Nv
        difr_matrix(n,m) =sum(exp(1i*(omega_v(m)*v_range-omega_t(n)*ti)));
    end
    if calc_error
        exp2 = exp(-1i*(omega_v*v_peak-omega_t(n)*L_det/v_peak));
        test_row(n) = sum(exp2.*difr_matrix(n,:))/Nv;
    end
end
if calc_error
    Err = sum(test_row)/Nt-1;
end
