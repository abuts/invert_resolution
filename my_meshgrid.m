function [XX,YY] = my_meshgrid(x,y)


Nx = numel(x);
Ny = numel(y);
if size(x,2)>1
    x = x';
end
if size(y,1)>1
    y = y';
end

XX= repmat(x,1,Ny);
YY = repmat(y,Nx,1);

