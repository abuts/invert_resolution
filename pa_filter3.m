function [rm,ft_reduced,omega_vt,omega_tt] = pa_filter3(res_matrix,ft_signal,omega_v,omega_t,eps,sam_spec)
% P-filter to work together with  InvertPulse3 to keep harmonics with the ampliture
% higher then the requested one

%[n_omega_t,n_omega_v] = size(res_matrix);

% max_el = max(abs(res_matrix),[],1);
% max_max = max(max_el);
% gr_eps_v = max_el >= eps*max_max;
% rm = res_matrix(:,gr_eps_v);
% omega_vt = omega_v(gr_eps_v);


max_el = max(abs(sam_spec),[],2);
max_max = max(max_el);
gr_eps_t = max_el >= eps*max_max;
rm = res_matrix(gr_eps_t,:);
omega_tt = omega_t(gr_eps_t);
ft_reduced = ft_signal(gr_eps_t);


max_el = max(abs(sam_spec),[],1);
gr_eps_v = max_el >= eps*max_max;

rm = rm(:,gr_eps_v);
omega_vt = omega_v(gr_eps_v);

