function test_deltaf
% let's test equations, used in fourier transformation of delta-functions
%

% check fft(Delta(v-V0)) = Sum_n{exp(-i*Omega_n*V0)}
Np=32;
iV = 5;
RangeMin = 5;
RangeMax =  8;
Range = RangeMax-RangeMin;
dx = Range/Np;
V0 = RangeMin+(iV-1)*dx;
xx = RangeMin+(0:Np-1)*dx;
x_edges = build_bins(xx);
Range = (max(x_edges)-min(x_edges));
%ind = 0:Np-1;
ind = fft_ind(Np);
Omega_n = ind*(2*pi/(Range));


DeltaF1 = zeros(1,Np);
DeltaF1(1,iV) = 1;
%DeltaF1(2,iV+4) = 1;
%omg = fftshift(Omega_n);
%DeltaF1 = exp(1i*omg(iV)*xx);

fft_spec = fft(DeltaF1');
rev_fft_data = ifft(fft_spec);
theor_sp = exp(-1i*(V0)*Omega_n);
%theor_sp = zeros(size(DeltaF1)); %
%theor_sp(iV) = Np*exp(1i*omg(iV)*Range0);
ref_theor_spec = ifft(theor_sp);
rev_theor_spec = (sum((exp(1i*(xx-V0).*Omega_n')),1))/Np;
[sft_omega,sft_spec] = sfft(xx,DeltaF1);
[sft_t,rev_sft_spec] = isfft(sft_omega,sft_spec,xx);
%theor_rev = (sum((exp(1i*(xx).*Omega_n')),1))/(Np-1);
assertElementsAlmostEqual(DeltaF1,rev_sft_spec')
plot(xx,real(DeltaF1),':p', xx,real(rev_fft_data),'--',xx,real(rev_theor_spec),'o',xx,real(ref_theor_spec),sft_t,real(rev_sft_spec),'>');
legend('original','rv fft','rv theor','rv theor using ifft','rv sft');
% if any(any(abs(imag(DeltaF1'))>1.e-9)) || any(any(abs(imag(rev_fft_data'))>1.e-9)) || any(abs(imag(ref_theor_spec))>1.e-9)
%     hold on;
%     plot(xx,imag(DeltaF1),':',xx,imag(rev_fft_data),'*',xx,imag(ref_theor_spec));
%     hold off;
% end



expT = IX_dataset_1d(fftshift(Omega_n'/(2*pi)),(180/pi)*fftshift(atan2(imag(fft_spec(:,1)),real(fft_spec(:,1)))));
%expT = IX_dataset_1d((Omega_n/(2*pi)),(180/pi)*(atan2(imag(spec(1,:)),real(spec(1,:)))));
acolor('r');
aline('-');
dl(expT);

theorR = IX_dataset_1d(fftshift(Omega_n/(2*pi)),(180/pi)*fftshift(atan2(imag(theor_sp),real(theor_sp))));
%theorR = IX_dataset_1d((Omega_n/(2*pi)),(180/pi)*(atan2(imag(theor_sp),real(theor_sp))));
acolor('b');
aline(':');
pl(theorR);

sft_s = IX_dataset_1d(fftshift(sft_omega'/(2*pi)),(180/pi)*fftshift(atan2(imag(sft_spec(:,1)),real(sft_spec(:,1)))));
%sft_s = IX_dataset_1d((omg1/(2*pi)),(180/pi)*(atan2(imag(sf1(1,:)),real(sf1(1,:)))));
acolor('g');
aline('-');
pd(sft_s );

[sft_spec,sft_omega]=p_filter1(sft_spec,sft_omega,Np/4);
[sft_tr,rev_sft_specr] = isfft(sft_omega,sft_spec,xx);
plot(sft_tr,rev_sft_specr);


L=10;
iV = 10;
V_i = L/RangeMax;
V_a = L/RangeMin;
dV = (V_a-V_i)/Np;
vi = V_i+(0:Np)*dV;
V0 = V_i+(iV-1)*dV;
vu = L./vi;
vu = sort(vu);

f = exp(-(vi-V0).^2/0.1);
fu =exp(-(vu-L/V0).^2/0.01);
fn0 = sprintf('direct and inverse distribution');
fh = findobj('type','figure', 'Name', fn0);
if  isempty(fh)
    figure('Name',fn0);
else
    hold off
    figure(fh);
end
%plot(vi,f,'r-o',vu,fu,'g-+');
%plot(vi,f,'r-o');
plot(vu,fu,'g-+');

[omg_dir,sft_dir_spec] = sfft(vi,f);
[omg_rev,sft_rev_spec] = sfft(vu,fu);
fn = sprintf('direct and inverse spectra');
fh1 = findobj('type','figure', 'Name', fn);
if  isempty(fh1)
    figure('Name',fn);
else
    hold off
    figure(fh1);
end
%plot(omg_dir/2*pi,abs(sft_dir_spec),'r-o',omg_rev,abs(sft_rev_spec),'g:+');
%plot(omg_dir/2*pi,abs(sft_dir_spec),'r-o');
plot(omg_rev,abs(sft_rev_spec),'g:+')

[vi,fdr] = isfft(omg_dir,sft_dir_spec,vi);
[vu,frr] = isfft(omg_rev,sft_rev_spec,vu);

fh = findobj('type','figure', 'Name', fn0);
figure(fh);
hold on
%plot(vi,fdr,'r:<',vu,frr,'g:>');
%plot(vi,fdr,'r:<');
plot(vu,frr,'g:>');
