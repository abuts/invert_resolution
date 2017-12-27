function TestFFT_convolution()
dt = 0.01;
t_min= -10;
t_max = 5;
dT = t_max-t_min;
t=t_min:dt:t_max;
s = signal(t);
Nt = numel(t);

plot(t,s);
hold on
pls = pulse(t);
plot(t,pls);

% fs = fft(s,2*Nt);
% fp = fft(pls,2*Nt);
fs = fft(s);
fp = fft(pls);

find = fft_ind(numel(fs));
phase = 1i*pi*find*(1-(t_min+t_max)/dT); % 

fc = fs.*fp.*exp(phase);
conv = ifft(fc)/dT;
plot(t,conv);
hold off

function pls = pulse(t)
t0 = 2;
sig1 = 0.05;
sig12 = sig1*sig1;
pls = abs(t).*exp(-(t-t0).^2/(2*sig12))/sqrt(2*pi*sig1);
%pls = fliplr(pls);



function sig = signal(t)

sig1 = 0.3;
sig12 = sig1*sig1;
x1 = 1;
sig2 = 0.1;
sig22 = sig2*sig2;

sig = ... %exp(-(t-x1).^2/(2*sig12))/sqrt(2*pi*sig1)+...
    exp(-(t+x1).^2/(2*sig12))/sqrt(2*pi*sig1)+exp(-(t).^2/(2*sig22))/sqrt(2*pi*sig2);







