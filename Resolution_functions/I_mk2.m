function int = I_mk2(k_v,n_t,v_min,v_max,v_index,t_index)
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
        AbsTol = 1.e-10;
        delta = 0.1*abs(nt)/2*pi;
        i_main = dir_int(kv,nt,v_min,delta,AbsTol);
        i_inv  = inv_int(kv,nt,v_min,delta,AbsTol);
        int =(i_main+i_inv)/(kv*nt);
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
function int1 = inv_int(kv,nt,v_min,delta,AbsTol)

lim_val = 0;
min_val = 1/delta;
if u0<max_lim
    funIS = @(u)(cos(2*pi*v_min*nt*u)./(u).^4);    
    intC = integral(funIS,u0,max_lim,'RelTol',0,'AbsTol',AbsTol);
else
    intC  = 0;
end
betta = v_min*nt;
alpha = kv;
albet =alpha/betta;
int22 = 2*albet*cos(2*pi*betta/delta)*delta^3+(4*albet/3-(2*pi*alpha)^2)*intC;
arg = 2*pi*nt*v_min;
int21 = 2*delta*cos(arg/delta)+1i*arg*(expint(-1i*arg/delta)-expint(1i*arg/delta));
int1   = 2*int20+int22+int21;
% else
%     funI = @(u)(cos(2*pi*(kv./u-v_min*nt*u))./(u.*u));
%     int = 2*integral(funI,2,Inf);
% end



function int0 = dir_int(kv,nt,v_min,delta,AbsTol)

funD = @(v)(cos(2*pi*(v/nt-(v_min/kv)./v)));
max_lim = 0.5*nt*kv;
fd = @(v)(v/nt-(v_min/kv)./v);
np = root_range(fd,delta,max_lim);
if ~isempty(np)
    det = sqrt(0.25*np.^2+nt*v_min/kv);
    shi = np*0.5;
    roots = [shi+det,shi-det];
    valid = roots>delta&roots<max_lim;
    wp = roots(valid);
    int0 =integral(funD,delta,0.5,'Waypoints',wp,'RelTol',0,'AbsTol',AbsTol);
else
    int0 = integral(funD,delta,0.5,'RelTol',0,'AbsTol',AbsTol);
end



function range = root_range(fun,min_val,max_val)

min_val = min(fun(min_val),fun(max_val));
max_val = max(fun(min_val),fun(max_val));

min_val_n = round(min_val);
if min_val_n<min_val
    min_val_n = min_val_n+1;
end
max_val_n = round(max_val);
if max_val_n>max_val
    max_val_n = max_val_n-1;
end
if min_val_n<=max_val_n
    range = min_val_n:1:max_val_n;
else
    range = [];
end

