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
        %         fp = 0.5*kv-2*nt*v_min;
        %         fun = @(u)(exp(2i*pi*u).*u./sqrt(u.*u+4*nt*kv*v_min));
        %         int = sin(2*pi*fp)/(2*pi) + integral(fun,-fp,fp)/kv;
        
%         funI = @(u)(cos(2*pi*(kv./u-v_min*nt*u))./(u.*u));
%         int1 = 2*integral(funI,2,Inf);
%         
%         funD = @(v)(cos(2*pi*(kv*v-v_min*nt./v)));
%         int0 = integral(funD,-0.5,0.5);
        AbsTol = 1.e-10;
        
        u0 = 2*pi*kv/v_min;
        delta = 1/u0;
        if delta<0.5
            fd = @(v)(kv*v-v_min*nt./v);
            np = root_range(fd,delta,0.5);
            if ~isempty(np)
                det = sqrt(0.25*np.^2/(kv*kv)+nt*v_min/kv);
                shi = np*0.25/kv;
                roots = [shi+det,shi-det];
                valid = roots>delta&roots<0.5;
                wp = roots(valid);
                int20 =integral(funD,delta,0.5,'Waypoints',wp,'RelTol',0,'AbsTol',AbsTol);
            else
                int20 = integral(funD,delta,0.5,'RelTol',0,'AbsTol',AbsTol);
            end
            
            
            funIS = @(u)(cos(2*pi*v_min*nt*u)./(u).^4);
            max_lim = (1/AbsTol)^(1/3);
            if u0<max_lim
                intC = integral(funIS,u0,max_lim);
            else
                intC  = 0;
            end
            betta = v_min*nt;
            alpha = kv;
            albet =alpha/betta;
            int21 = 2*albet*cos(2*pi*betta/delta)*delta^3+(4*albet/3-(2*pi*alpha)^2)*intC;
            
            int = 2*(int20+int21);
        else
            funI = @(u)(cos(2*pi*(kv./u-v_min*nt*u))./(u.*u));
            int = 2*integral(funI,2,Inf);            
        end
        
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

