function ind = fft_ind(Np)
% function returns indexes of the fourier frequencies for a transformation
% calculated using fft with NP points
%
% a correspondent frequency would be 2*pi*ind/DT where DT is the interval
% of the transformation. 
%
Np2 = floor(Np/2);
if  rem(Np,2) >0
    ind = [0:Np2,-Np2:1:-1];
else
    ind = [0:Np2-1,-Np2:1:-1];
end





