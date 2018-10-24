function test_deltaf2D
% let's test equations, used in fourier transformation of delta-functions
%

% check fft(Delta(v-V0)) = Sum_n{exp(-i*Omega_n*V0)}
Nt=7;
Nv = 8;
IT = 2;
IV = 4;

TRangeMin = 1;
TRangeMax =  3;
VRangeMin = 20;
VRangeMax = 40;


dT = TRangeMax-TRangeMin;
dt = dT/Nt;
t = TRangeMin+(0:Nt-1)*dt;
dV = VRangeMax-VRangeMin;
dv = dV/Nv;
v = VRangeMin+(0:Nv-1)*dv;
DeltaF2 = zeros(Nt,Nv);
DeltaF2(IT,IV) = 1;



[omegaT,omegaV,spec] = sft2(t,v,DeltaF2);


[t1,v1,f] = isft2(omegaT,omegaV,spec,t,v);

assertElementsAlmostEqual(DeltaF2,f)
assertElementsAlmostEqual(t,t1');
assertElementsAlmostEqual(v,v1);


