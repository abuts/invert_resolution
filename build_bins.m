function [bin_edges,bin_sizes] = build_bins(bin_centers)
% Auxiliary function to generate bin edges given bin centers
if size(bin_centers,1)~=1 % move bin centers along rows
    bin_centers = bin_centers';
end

bin_edges = 0.5*(bin_centers(1:end-1)+bin_centers(2:end));

bin_edges = [2*bin_centers(1)-bin_edges(1),bin_edges,2*bin_centers(end)-bin_edges(end)];

if nargout > 1
    bin_sizes = bin_edges(2:end)-bin_edges(1:end-1);
end


