function test_fftf_1D
% testing filter, related to standard fft transformation.

t = (1:1000)*0.01;

fd = fdistr(t);
plot(t,fd);
distr = build_distribution(t,fd,1000);

[bedg,bin] = build_bins(t);
fc = histcounts(distr,bedg);
fc = fc./bin;
hold on;
plot(t,fc);


[X, f, y, y2] = fftf_1D(t, fc, [],20);

[X, f, y, y2] = fftf_1D(t, fc, [],15);
[X, f, y, y2] = fftf_1D(t, fc, [],50);








function fd = fdistr(t)

t_min = min(t);
t_max = max(t);
dT = t_max-t_min;
t0 =[0.1,0.55,0.7]*dT;
Sigma = [sqrt(0.0002),sqrt(0.0010),sqrt(0.0050)]*dT;
Sigma2 =Sigma.*Sigma;
A = [2,3,5];
fd = sum(exp(-(t-t0').^2./(2*Sigma2')).*(A'./Sigma')/sqrt(2*pi),1);
