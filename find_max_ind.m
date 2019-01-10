function [n_max,m_max]  = find_max_ind(resolution_matrix,eps)
[N,M] = size(resolution_matrix);

[max_el,m_index] = max(abs(resolution_matrix),[],2);
max_max = max(max_el);

gr_eps = max_el >= eps*max_max;
n_max  = sum(gr_eps);
m_max  = m_index(n_max);