function distr = build_distribution(bin_centers,f_distr,num_points)


[~,bins] = build_bins(bin_centers);
norma   = abs(f_distr)*bins';
f_distr = abs(f_distr)/norma;
prod_function = cumsum(f_distr.*bins);
%Nip = numel(f_distr);

norm_dist = rand([num_points,1]);
[prod_function,ind_uniq] = unique(prod_function);
bin_centers = bin_centers(ind_uniq);
distr = interp1(prod_function,bin_centers,norm_dist,'linear',0);

%histogram(distr);