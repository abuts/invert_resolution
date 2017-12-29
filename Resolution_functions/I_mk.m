function int = I_mk(k_v,n_t,v_min,v_max,v_index,t_index)
% v_min == v_min[cm/sec]/DV[cm/sec]
% v_max == v_min[cm/sec]/DV[cm/sec]
kv = v_index(k_v);
nt = t_index(n_t);

if kv~=0
    fun0 = @(v)(exp(2i*pi*(kv*v-v_min*nt./v)));
    fun = @(v)(exp(2i*pi*(kv*v-v_min*nt./v))./(v.*v));
    v_sq = v_min*nt/kv;
    if v_sq>0
        vm0 = sqrt(v_sq);
        if vm0>v_min && vm0<v_max
            single_integral = false;
        end
    end
else
    if nt == 0
        int = v_max-v_min;
    else
%         fun = @(v)(exp(-2i*pi*v_min*nt./v));    
%         int = integral(fun,v_min,v_max);
        v_r = v_min/v_max;
        int = v_max*exp(2i*pi*nt*v_r) - v_min*exp(2i*pi*nt)+...
            2i*pi*nt*v_min*(expint(-2i*pi*nt*v_r)...
            -expint(-2i*pi*nt));
    end
end

