function [rm,ft_reduced,omega_vt,omega_tt] = p_filter3(res_matrix,ft_signal,omega_v,omega_t,n_harm_left)
% P-filter to work together with  InvertPulse3 to keep selected number of
% harmonics

[n_omega_t,n_omega_v] = size(res_matrix);
[nt_block,nv_block] = p_filter_block(n_omega_t,n_omega_v,50000,n_harm_left);    



rm = res_matrix(nt_block,nv_block);
ft_reduced = ft_signal(nt_block);


omega_vt  = omega_v(nv_block);
omega_tt  = omega_t(nt_block);



