function [i0,i1] = TestExpansions()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


delta = 0.0001;
v_min = 0.11;
Lim = 10000;
Lim0 = 1/delta;
AbsErr = 1.e-12;
kv = 100;
nt = 100;

fb = @(v)(exp(2i*pi*(kv*v-v_min*nt./v)));

i1 = integral(fb,-0.5,-delta/(1+2*delta))+integral(fb,delta/(1+2*delta),0.5);

betta = v_min*nt*kv;
fn = @(u)(cos(2*pi*(u-betta./u)));

i2 = (2/abs(kv))*integral(fn,delta/(1+2*delta)*abs(kv),0.5*abs(kv));


fprintf(' zero difr: %d\n',i1-i2);
% ibes = ass_exp(kv,nt,v_min);
% fprintf(' bes difr: %d, i_num: %e, i_bes:%e \n',i2-ibes,i2,ibes);

%
betta = 4*v_min*kv*nt;
u0 = sqrt(abs(betta)); % the point where the module of the direct term equal to the module of the reverse term

N0 = ceil(u0);
if N0>abs(kv)
    N0 = abs(kv);
end
dir_ints = dir_int(kv,N0,betta,AbsErr);
inv_ints = inv_int(N0,betta,AbsErr);
i2p = (1/abs(kv))*(sum(dir_ints)+sum(inv_ints));


fprintf(' fin difr: %d\n',i2p-i2);

assertEqual(real(i1),i2p);

f0 = @(v)(v.*(1+1i*pi*alpha/2*v).*exp(-2i*pi*betta./v));

i0 = 2i*pi*alpha*(integral(f0,-delta,-1/Lim)+integral(f0,1/Lim,delta));




%i1 = 4*alpha/betta*cos(2*pi*betta*Lim0)/Lim0^3;
f1cos = @(u)(cos(2*pi*betta*u)./u.^4);
ic = integral(f1cos,Lim0,Lim);

albet =alpha/betta;
i1 = 2*albet*cos(2*pi*betta/delta)*delta^3+(4*albet/3-(2*pi*alpha)^2)*ic;
%
% fi1 = @(v)(v.*exp(-2i*pi*betta./v));
% ip1i = 2i*pi*alpha*(integral(fi1,-delta,-1/Lim)+integral(fi1,1/Lim,delta));
% ip1r = (alpha/betta)*(2*cos(2*pi*betta*Lim0)/Lim0^3+4/3*integral(f1cos,Lim0,Lim));

% fp2  = @(u)(sin(2*pi*betta*u)./u.^3);
% ip1  = 4*pi*alpha*integral(fp2,Lim0,Lim);

% fi1 = @(v)(v.^2.*exp(-2i*pi*betta./v));
% ip1 = (integral(fi1,-delta,-1/Lim,'RelTol',0,'AbsTol',1.e-12)+integral(fi1,1/Lim,delta,'RelTol',0,'AbsTol',1.e-12));
% ipp2 = 2*integral(f1cos,Lim0,Lim,'RelTol',0,'AbsTol',1.e-12);


%ip2 = -4*(alpha/betta)*cos(2*pi*betta*Lim0)/Lim0^3+12*integral(f1,Lim0,Lim);
function inv_ints = inv_int(N0,betta,AbsErr)

u0 = 1/N0;
phase = betta*u0;
step = 10/abs(betta);
N_max = 1000;
inv_ints = zeros(1,N_max);
finv = @(u,Del)(cos(pi*(1./(Del+u)-phase-betta*u))./(Del+u).^2);
for i=1:N_max
    Del_k = u0+(i-1)*step;
    inv_ints(i) = integral(@(u)finv(u,Del_k),0,step,'RelTol',0,'AbsTol',AbsErr);
    if abs(inv_ints(i))<0.1*AbsErr
        break;
    end
end



function ints = dir_int(kv,N0,betta,AbsErr)
Np = abs(kv)-N0;
if Np<=0
    ints = 0;
    return
end
ints = zeros(1,Np);
f_dir = @(u,Del)(cos(pi*(u-betta./(Del+u))));
for k=1:Np
    Del_k = N0+k-1;
    ints(k) = (-1)^(N0+k-1)*integral(@(u)f_dir(u,Del_k),0,1,'RelTol',0,'AbsTol',AbsErr);
end







function int = inv_int0(betta,AbsErr,delta)


if exist('delta','var')
    MaxN = ceil(1/delta/abs(betta));
else
    MaxN = ceil(1/AbsErr/abs(betta));
end
if MaxN > 1000000
    MaxN = 1000000;
end

step = 1;
phase_shift = 1/betta;
step0 = 100*phase_shift;
x_lim = 2.+step0;
f0inv = @(u)(cos(2*pi*(1./u-betta*u))./(u.^2));
int0  = integral(f0inv,2,x_lim ,'RelTol',0,'AbsTol',AbsErr);


sec_term = @(n)(abs(dec_sq(x_lim,betta,step,n))-AbsErr);
if sec_term(x_lim)>0 && sec_term(MaxN)<0
    nit = ceil(fzero(sec_term,[x_lim,MaxN]));
else
    nit = MaxN-floor(step0)+1;
end



f1inv = @(u,Dlt)(cos(2*pi*(1./(u+Dlt)-betta*(u+x_lim)))./((Dlt+u).^2));
if nit*betta>x_lim
    
    i = 1:nit;
    [sub_sum1,sub_sum2] = dec_0and1(x_lim,betta,step,i);
    sub_sumI =zeros(1,nit) ;
    for i=1:nit
        Dlt_k = x_lim+(i-1)*step/abs(betta);
        int1 = integral(@(u)f1inv(u,Dlt_k),0,step/abs(betta),'RelTol',0,'AbsTol',1.e-8);
        sub_sumI(i) = int1;
        %         %sum = sum+int1;
        %         %intL = psi(4,i+1)*(2*pi)^2*betta^3/24;
        %         % assymptotic value fr the
        %         arg = 2*pi/Dlt_k;
        %         Ass0 =step/pi/betta^2/Dlt_k^3*(sin(arg)+pi/Dlt_k*cos(arg ));
        %         %
        if abs(int1-sub_sum1(i)-sub_sum2(i))<0.1*AbsErr
            break
        end
    end
    ii = (int0+sum(sub_sumI));
    int = (int0+sum(sub_sum1+sub_sum2));
end
function fi = int_fun(u,Dlt_k,betta,x0)
fi = (cos(2*pi./(Dlt_k+u)).*cos(2*pi*betta*(u+x0))+sin(2*pi./(Dlt_k+u)).*sin(2*pi*betta*(u+x0)))./(Dlt_k+u).^2;
function fi = dec_fun(u,Dlt_k,betta,x0)
fi = (cos(2*pi./(Dlt_k+u)).*cos(2*pi*betta*(u+x0))+sin(2*pi./(Dlt_k+u)).*sin(2*pi*betta*(u+x0)))./(Dlt_k+u).^2;


function val = dec_0(x0,betta,step,niter)
Dlt_k = x0+niter*(step/abs(betta));
arg = 2*pi*(1./Dlt_k);
val = (step/pi/betta^2)*(sin(arg)+(pi./Dlt_k).* cos(arg))./(Dlt_k2.*Dlt_k);



function val = dec_sq(x0,betta,step,niter)
Dlt_k = x0+niter*(step/abs(betta));
arg = 2*pi*(1./Dlt_k);
sina = sin(arg);
cosa = cos(arg);
Dlt_k2 = Dlt_k.^2;

fs = sina.*(3-2*pi*pi./Dlt_k2)+6*pi*cosa./Dlt_k;
fc = cosa.*(3-2*pi*pi./Dlt_k2)-6*pi*sina./Dlt_k;

val =(0.5*step/pi/betta^3)*(fc/pi-step*fs)./(Dlt_k2.^2);

function [val1,val2] = dec_0and1(x0,betta,step,niter)
Dlt_k = x0+(niter-1)*(step/abs(betta));
arg = 2*pi./Dlt_k;
sina = sin(arg);
cosa = cos(arg);
Dlt_k2 = Dlt_k.^2;

fs = (3-2*pi*pi./Dlt_k2).*sina+6*pi*cosa./Dlt_k;
fc = (3-2*pi*pi./Dlt_k2).*cosa-6*pi*sina./Dlt_k;

val1 = (step/pi/betta^2)*(sina+(pi./Dlt_k).*cosa)./(Dlt_k2.*Dlt_k);

val2 =(0.5*step/pi/betta^3)*(fc/pi-step*fs)./(Dlt_k2.^2);

function val = ass_exp(kv,nt,v_min)

mu = 2/abs(nt); %+1.e-10i;
betta= 2*pi*abs(kv*nt)*sqrt(v_min*abs(nt/kv));
if nt/kv>0
    ival = -2i*betta*besselk(1,mu*betta); 
else
    ival = 2*betta*besselk(1,1i*mu*betta); 
end
val = real(ival)*(2/abs(kv*nt));
