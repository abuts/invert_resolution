function [rm,ft_reduced,omega_vt,omega_tt] = p_filter3b(input_spectra,ft_signal,omega_v,omega_t,n_harm_left)
% P-filter to work together with  InvertPulse3 to keep selected number of
% harmonics

[n_omega_t,n_omega_v] = size(input_spectra);
n_omega_t_range = get_cental_block(n_omega_t,n_harm_left);
n_omega_v_range = get_cental_block(n_omega_v,n_harm_left);



rm = fftshift(input_spectra);
rm = fftshift(rm(n_omega_t_range,n_omega_v_range));

ft_reduced = fftshift(ft_signal);
ft_reduced = fftshift(ft_reduced(n_omega_t_range));

omega_vt  = fftshift(omega_v);
omega_vt  = fftshift(omega_vt(n_omega_v_range));
omega_tt  = fftshift(omega_t);
omega_tt  = fftshift(omega_tt(n_omega_t_range));


function block = get_cental_block(n_total,n_selected)
left_cent = floor(n_total/2);
if n_selected>= left_cent
    block = 1:n_total;
    return;
end
start_s  = left_cent-n_selected+1;
if rem(n_total,2) > 0
    end_s    = left_cent+1+n_selected;
else
    end_s    = left_cent+n_selected;
end
block = start_s:end_s;
