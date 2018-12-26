function [Sm,omega] = p_filter1(Sm,omega,n_harm_left)
% P-filter to work together with  InvertPulse2 to keep selected number of
% harmonics

n_harm_exisist = numel(Sm);
if  rem(n_harm_exisist,2) >0
    left = n_harm_left+1;
    right= n_harm_left;
    
    piece_to_leave = 2*n_harm_left+1;
else
    left = n_harm_left;
    right= n_harm_left;
    
    piece_to_leave = 2*n_harm_left;
end
if piece_to_leave >= n_harm_exisist
    return;
end


Sm = [Sm(1:left);Sm(end-right-1:end)];
omega = [omega(1:left),omega(end-right-1:end)];


