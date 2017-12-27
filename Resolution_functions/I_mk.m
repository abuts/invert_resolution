function int = I_mk(k_v,n_t,v_min,v_max,v_index,t_index)
% v_min == v_min[cm/sec]/DV[cm/sec]
% v_max == v_min[cm/sec]/DV[cm/sec]
kv = v_index(k_v);
nt = t_index(n_t);
fun = @(v)(exp(2i*pi*(kv*v-v_min*nt./v)));
single_integral = true;

if kv~=0
    v_sq = v_min*nt/kv;
    if v_sq>0
        vm0 = sqrt(v_sq);
        if vm0>v_min && vm0<v_max
            single_integral = false;
        end
    end
end
if single_integral
    int = integral(fun,v_min,v_max);
else
    int = integral(fun,v_min,vm0)+integral(fun,vm0,v_max);
end
