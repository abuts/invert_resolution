function int = I_mk3(k_v,n_t,v_min,v_max,v_index,t_index)
% v_min == v_min[cm/sec]/DV[cm/sec]
% v_max == v_min[cm/sec]/DV[cm/sec]

kv = v_index(k_v);
nt = t_index(n_t);

if kv~=0
    if nt==0
        int = 0;
        %int = (exp(1i*pi*kv)-exp(-1i*pi*kv))/(2i*pi*kv);
        %         fun0 = @(v)(exp(2i*pi*(kv*v)));
        %         int0 = integral(fun0,-0.5,0.5);
        
    else
        
        %         funD = @(v)(cos(2*pi*(v/nt-(v_min/kv)./v)));
        %         funI = @(u)(cos(2*pi*(kv./u-v_min*nt*u))./(u.*u));
        %         int1 = 2*integral(funI,2,Inf);
        %         int0 = integral(funD,-0.5,0.5);
        int = interator(kv,nt,v_min);
    end
else % kv==0
    if nt == 0
        int = 1;
    else
        arg = 2*pi*nt*v_min;
        int = cos(2*arg)+1i*arg*(expint(-2i*arg)-expint(2i*arg));
        
        %
        %         fun0 = @(v)(exp(-2i*pi*v_min*nt./v));
        %         int0 = integral(fun0,-0.5,0.5);
        
    end
end

function int = interator(kv,nt,v_min)

AbsErr = 1.e-12;
betta = 4*v_min*kv*nt;
u0 = sqrt(abs(betta)); % the point where the module of the direct term equal to the module of the reverse term

N0 = ceil(u0);
if N0>abs(kv)
    N0 = abs(kv);
end
dir_ints = dir_int(kv,N0,betta,AbsErr);
inv_ints = inv_int(N0,betta,AbsErr);
int = (1/abs(kv))*(sum(dir_ints)+sum(inv_ints));


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



