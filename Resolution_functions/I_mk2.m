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
        
        fun0 = @(v)(exp(2i*pi*(kv*v-v_min*nt./v)));
        int = integral(fun0,-0.5,0.5);
        
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

