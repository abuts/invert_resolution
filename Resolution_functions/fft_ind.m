function ind = fft_ind(Np)
% function returns indexes of the fourier components calculated usign fft 
% for the transformation with NP points
%
Np2 = floor(Np/2);
if  rem(Np,2) >0
    ind = [0:Np2,-Np2:1:-1];
else
    ind = [0:Np2-1,-Np2:1:-1];
end





