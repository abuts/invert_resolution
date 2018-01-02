function TestFFT_convolution()
dt = 0.01;
t_min= -3;
t_max = 20;
dT = t_max-t_min;
t=t_min:dt:t_max;
s = signal(t);
Nt = numel(t);
DN = Nt;

plot(t,s);
hold on
pls = pulse(t);
plot(t,pls);

Ntp = Nt+DN;
fs = fft(s,Ntp);
fp = fft(pls,Ntp);

% i = 0:Ntp-1;
% te = t_min + i*dt;

% fst = ifft(fs);
% fpt = ifft(fp);
% hold on
% plot(te,fst)
% plot(te,fpt);


% fs = fft(s);
% fp = fft(pls);

find = fft_ind(numel(fs));
phase = 1i*pi*find*(1-(2*t_min+dt*(Ntp-1))/(dt*(Ntp-1))); %

fc = fs.*fp.*exp(phase);
conv = ifft(fc)/dT;
if numel(t) ~= numel(conv)
    i = 0:Ntp-1;
    t = t_min + i*dt;
end
plot(t,real(conv));
if ~any(isreal(conv))
    imp = imag(conv);
    realmax = max(real(conv));
    imgmax = max(imp);
    if imgmax<10*realmax
        fprintf('imaginary part renormalized from %e to %e\n',imgmax,0.1*realmax);
        imp = imp*(0.1*realmax/imgmax);
    end
    plot(t,imp);
end
hold off

function pls = pulse(t)
t0 = 10;
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







