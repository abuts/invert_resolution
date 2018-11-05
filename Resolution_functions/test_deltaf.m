function test_deltaf
% let's test equations, used in fourier transformation of delta-functions
%

% check fft(Delta(v-V0)) = Sum_n{exp(-i*Omega_n*V0)}
Np=32;
iV = 3;
RangeMin = 2;
RangeMax =  3;
Range0 = 0.5*(RangeMax+RangeMin);
Range = RangeMax-RangeMin;
dx = Range/Np;
V0 = RangeMin+(iV-1)*dx;
xx = RangeMin+(0:Np-1)*dx;
%ind = 0:Np-1;
ind = fft_ind(Np);
Omega_n = ind*(2*pi/(Range));


DeltaF1 = zeros(2,Np);
DeltaF1(1,iV) = 1;
DeltaF1(2,iV+4) = 1;
%omg = fftshift(Omega_n);
%DeltaF1 = exp(1i*omg(iV)*xx);

fft_spec = fft(DeltaF1);
rev_fft_data = ifft(fft_spec);
theor_sp = exp(-1i*(V0-Range0)*Omega_n);
%theor_sp = zeros(size(DeltaF1)); %
%theor_sp(iV) = Np*exp(1i*omg(iV)*Range0);
ref_theor_spec = ifft(theor_sp);
rev_theor_spec = (sum((exp(1i*(xx-V0).*Omega_n')),1))/Np;
[sft_omega,sft_spec] = sft(xx,DeltaF1,ind);
[sft_t,rev_sft_spec] = isft(sft_omega,sft_spec,xx);
%theor_rev = (sum((exp(1i*(xx).*Omega_n')),1))/(Np-1);
assertElementsAlmostEqual(DeltaF1,rev_sft_spec)
plot(xx,real(DeltaF1),':', xx,real(rev_fft_data),'*',xx,real(rev_theor_spec),'o',xx,real(ref_theor_spec),sft_t,real(rev_sft_spec),'+');
legend('original','reversed fft','reversed theor','rev theor using ifft','reversed sft');
% if any(any(abs(imag(DeltaF1'))>1.e-9)) || any(any(abs(imag(rev_fft_data'))>1.e-9)) || any(abs(imag(ref_theor_spec))>1.e-9)
%     hold on;
%     plot(xx,imag(DeltaF1),':',xx,imag(rev_fft_data),'*',xx,imag(ref_theor_spec));
%     hold off;
% end



expT = IX_dataset_1d(fftshift(Omega_n/(2*pi)),(180/pi)*fftshift(atan2(imag(fft_spec(1,:)),real(fft_spec(1,:)))));
%expT = IX_dataset_1d((Omega_n/(2*pi)),(180/pi)*(atan2(imag(spec(1,:)),real(spec(1,:)))));
acolor('r');
aline('-');
dl(expT);

theorR = IX_dataset_1d(fftshift(Omega_n/(2*pi)),(180/pi)*fftshift(atan2(imag(theor_sp),real(theor_sp))));
%theorR = IX_dataset_1d((Omega_n/(2*pi)),(180/pi)*(atan2(imag(theor_sp),real(theor_sp))));
acolor('b');
aline(':');
pl(theorR);

sft_s = IX_dataset_1d(fftshift(sft_omega/(2*pi)),(180/pi)*fftshift(atan2(imag(sft_spec(1,:)),real(sft_spec(1,:)))));
%sft_s = IX_dataset_1d((omg1/(2*pi)),(180/pi)*(atan2(imag(sf1(1,:)),real(sf1(1,:)))));
acolor('g');
aline('-');
pl(sft_s );



