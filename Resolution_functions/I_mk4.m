function int = I_mk4(k_v,n_t,v_min,v_max,v_index,t_index)
% v_min == v_min[cm/sec]/DV[cm/sec]
% v_max == v_min[cm/sec]/DV[cm/sec]

kv = v_index(k_v);
nt = t_index(n_t);
AbsErr = 1.e-12;
betta = v_min*v_max*nt;

if kv~=0
    if nt==0
        int = 0; % not obvious, but the row below is always = 0
        %int = (exp(2i*pi*kv*v_max)-exp(2i*pi*kv*v_min))/(2i*pi*kv);
        %         fun0 = @(v)(exp(2i*pi*(kv*v)));
        %         int0 = integral(fun0,-0.5,0.5);
        
    else
        
        
        %         funI = @(u)(cos(2*pi*(kv./u-v_min*nt*u))./(u.*u));
        %         int1 = 2*integral(funI,2,Inf);
        %         int0 = integral(funD,-0.5,0.5);
        funD = @(v)(exp(2i*pi*(kv*v-(betta./v))));
        int = integral(funD,v_min,v_max,'RelTol',0,'AbsTol',AbsErr);
        %int = interator(kv,nt,v_min,AbsErr);
    end
else % kv==0
    if nt == 0
        int = 1;
    else
        funI = @(u)(exp((-2i*pi*betta)*u)./u.^2);
        
        int = integral(funI,1/v_max,1/v_min,'RelTol',0,'AbsTol',AbsErr);
        %         betta = nt*v_min;
        %         ints = inv_int0(2,betta,AbsErr);
        %         int = sum(ints);
        %         %arg = 2*pi*nt*v_min;
        %int = cos(2*arg)+1i*arg*(expint(-2i*arg)-expint(2i*arg));
        
        %
        %         fun0 = @(v)(exp(-2i*pi*v_min*nt./v));
        %         int0 = integral(fun0,-0.5,0.5);
        
    end
end

