function [nt_block,nv_block] = p_filter_block(Nt,Nv,n_left_t,n_left_v)
% P-filter to work together with  InvertPulse3 to keep selected number of
% harmonics


nt_block = get_edge_block(Nt,n_left_t);
nv_block = get_edge_block(Nv,n_left_v);




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

