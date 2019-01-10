function [res_matrix,ft_reduced,omega_vt,omega_t] = p_filter3a(res_matrix,ft_signal,omega_v,omega_t,n_harm_left)
% P-filter to work together with  InvertPulse3 to keep selected number of
% harmonics

n_omega_v = numel(ft_signal);
n_omega_v_range = get_edge_block(n_omega_v,n_harm_left);


ft_reduced = zeros(size(ft_signal));
ft_reduced(n_omega_v_range) = ft_signal(n_omega_v_range);


omega_vt  = omega_v;


function block = get_edge_block(n_total,n_selected)
left_cent = floor(n_total/2);
if n_selected>= left_cent
    block = 1:n_total;
    return;
end

if rem(n_total,2) > 0
    block = [1:n_selected+1,n_total-n_selected+1:n_total];
else
    block = [1:n_selected,n_total-n_selected+1:n_total];
end

