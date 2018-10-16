function co = test_direct_conditionality

Np=11; % do it aways odd
Np2=floor(Np/2);
v_min = 1.5;
v_max = 2;
dV = (v_max-v_min);
dv = dV/Np;
vi = v_min+dv*(0:Np-1);
Omega_n = (2*pi/dV)*(-Np2:Np2);

N_vi = 5;
V0 = vi(N_vi);

defmat = exp(1i*Omega_n.*(vi'-V0));
fi = sum(defmat,2);
co = cond(defmat);
[sv,U,V] = svd(defmat);

