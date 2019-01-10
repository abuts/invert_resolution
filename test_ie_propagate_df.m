function test_ie_propagate_df
%
% calculations according to T1B -pp1:3 document from 2018/10/01
%
v_min = 0.5;
v_max = 2;


dV = v_max - v_min;
v_avrg = 0.5*(v_min+v_max)/dV;
Nu_v = v_avrg*v_avrg-0.25;

Dp = .3;
Np=33;
%f_ind = fft_ind(Np);
f_ind = -floor(Np/2):floor(Np/2);

I_nm = zeros(Np,Np);
fn = zeros(Np,1);
sumN = zeros(1,Np);
Sm   = zeros(Np,1);
for n=1:Np
    fprintf('step: %d\n',n);
    for m=1:Np
        I_nm(m,n) = I_mk5(f_ind(m),f_ind(n),v_avrg);
        %I_nm(m) = I_mk5(f_ind(m),f_ind(n),v_avrg);
        Sm(m) = exp(-2i*m*(Dp));
    end
    sumN(n) = sum(I_nm(:,n).*Sm);
    fn(n) = exp(-2i*pi*n*(Nu_v/(v_avrg+Dp)));
end
figure
plot(f_ind,(180/pi)*atan2(imag(sumN),real(sumN)),f_ind,(180/pi)*atan2(imag(fn),real(fn)));
figure
plot(f_ind,abs(sumN),'o',f_ind,abs(fn));



f_vel_sp  = linsolve(I_nm,fn);
plot(f_ind,abs(Sm),'o',f_ind,abs(f_vel_sp));

dv = dV/(Np-1);
v = dv*f_ind;
Omega_v = (2*pi/dV)*f_ind;
Ex_dis  = exp(1i*(v).*Omega_v'); % this gives omega change along columns (first dimension)
vel_distR = (sum((repmat(f_vel_sp,1,Np).*Ex_dis),1))/(Np-1);
vel_distO = (sum((repmat(Sm,1,Np).*Ex_dis),1))/(Np-1);
plot(v,vel_distR,v,vel_distO);

