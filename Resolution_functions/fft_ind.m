function [ind,Np2] = fft_ind(Np)
% function returns indexes of the fourier frequencies for a transformation
% calculated using fft with NP points
%
% a correspondent frequency would be 2*pi*ind/DT where DT is the interval
% of the transformation.
%
%ind2 = 0:Np-1;
% if  rem(Np,2) >0
%     Ne=Np;
% else
%     Ne= Np-1;
% end
% ind = [1:Ne;-(1:Np-1)];
% ind = [0,reshape(ind,1,numel(ind))];

if  rem(Np,2) >0
    Np2 = floor((Np-1)/2);
    ind = [0:Np2,-Np2:1:-1];
else
    Np2 = floor(Np/2);
    ind = [0:Np2-1,-Np2:1:-1];
end





