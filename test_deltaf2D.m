function test_deltaf2D
% let's test equations, used in fourier transformation of delta-functions
%

% check fft(Delta(v-V0)) = Sum_n{exp(-i*Omega_n*V0)}
Nt=14;
Nv = 20;
IT = 5;
IV = 4;

TRangeMin = 1;
TRangeMax =  4;
VRangeMin = 10;
VRangeMax = 20;


dT = TRangeMax-TRangeMin;
dt = dT/Nt;
t = TRangeMin+(0:Nt-1)*dt;
T0 = TRangeMin+(IT-1)*dt; 
dV = VRangeMax-VRangeMin;
dv = dV/Nv;
v = VRangeMin+(0:Nv-1)*dv;
V0 = VRangeMin+(IV-1)*dv; 

DeltaF2 = zeros(Nt,Nv);
DeltaF2(IT,IV) = 1;



[omegaT,omegaV,spec] = sfft2(t,v,DeltaF2);
[omegaT1,omegaV1,spec1] = sft2(t,v,DeltaF2);
assertElementsAlmostEqual(spec,spec1);
assertElementsAlmostEqual(omegaT,omegaT1);
assertElementsAlmostEqual(omegaV,omegaV1);

[t1,v1,f] = isfft2(omegaT,omegaV,spec,t,v);

assertElementsAlmostEqual(DeltaF2,f)
assertElementsAlmostEqual(t,t1');
assertElementsAlmostEqual(v,v1);

L=10;

Tii = L/TRangeMax;
Tia = L/TRangeMin;
%dti = (Tia-Tii)/Nt;
%ti = Tii+(0:Nt-1)*dti;
ti = L./t;


Vii = L/VRangeMax;
Via = L/VRangeMin;
%dvi = (Via -Vii)/Nv;
%vi = Vii+(0:Nv-1)*dvi;
vi = L./v;

[xi,yi] = meshgrid(ti,vi);
f = exp(-(xi'-L/T0).^2/0.4 -(yi'-L/V0).^2/0.002);

[omegaT,omegaV,spec] = sfft2(ti,vi,f);
[omegaTp,omegaVp,spec1] = sft2(ti,vi,f);

[t1,v1,fr] = isfft2(omegaT,omegaV,spec,ti,vi);
assertElementsAlmostEqual(f,fr);
