function int = I_mk(k_v,n_t,v_min,v_max,v_index,t_index)
% v_min == v_min[cm/sec]/DV[cm/sec]
% v_max == v_min[cm/sec]/DV[cm/sec]
kv = v_index(k_v);
nt = t_index(n_t);

if kv~=0
    if nt==0
        int = (exp(2i*pi*kv*v_max)-exp(2i*pi*kv*v_min))/(2i*pi*kv);
        %         fun0 = @(v)(exp(2i*pi*(kv*v)));
        %         int0 = integral(fun0,v_min,v_max);
        
    else
        v_r = v_min/v_max;
        min_val = min(kv*v_min-v_min*nt/v_min,kv*v_max-v_r*nt);
        max_val = max(kv*v_min-v_min*nt/v_min,kv*v_max-v_r*nt);
        min_val_n = round(min_val);
        if min_val_n<min_val
            min_val_n = min_val_n+1;
        end
        max_val_n = round(max_val);
        if max_val_n>max_val
            max_val_n = max_val_n-1;
        end
        fun = @(v)(exp(2i*pi*(kv*v-v_min*nt./v))./(v.*v));
        int = (exp(2i*pi*(kv*v_max-nt*v_r))-exp(2i*pi*(kv*v_min-nt)))/(2i*pi*kv);
        if min_val_n<=max_val_n
            np = min_val_n:1:max_val_n;
            det = sqrt(0.25*np.^2/(kv*kv)+nt*v_min/kv);
            shi = np*0.25/kv;
            roots = [shi+det,shi-det];
            valid = roots>v_min&roots<v_max;
            wp = roots(valid);
            int = int-(nt*v_min/kv)*integral(fun,v_min,v_max,'Waypoints',wp,'RelTol',0,'AbsTol',1.e-10);
        else
            int = int-(nt*v_min/kv)*integral(fun,v_min,v_max,'RelTol',0,'AbsTol',1.e-10);
        end
        
        %                 fun0 = @(v)(exp(2i*pi*(kv*v-v_min*nt./v)));
        %                 int0 = integral(fun0,v_min,v_max);
        
    end
else
    if nt == 0
        int = v_max-v_min;
    else
        v_r = v_min/v_max;
        min_val = min(-nt,-v_r*nt);
        max_val = max(-nt,-v_r*nt);
        min_val_n = round(min_val);
        if min_val_n<min_val
            min_val_n = min_val_n+1;
        end
        max_val_n = round(max_val);
        if max_val_n>max_val
            max_val_n = max_val_n-1;
        end
          
   
        fun = @(v)(exp(-2i*pi*v_min*nt./v));        
        if min_val_n<=max_val_n
            np = min_val_n:1:max_val_n;
            roots = v_min*nt./np;
            valid = roots>v_min&roots<v_max;
            wp = roots(valid);
            int = integral(fun,v_min,v_max,'Waypoints',wp,'RelTol',0,'AbsTol',1.e-10);
        else
            int = integral(fun,v_min,v_max,'RelTol',0,'AbsTol',1.e-10);
        end        
        %         int = v_max*exp(2i*pi*nt*v_r) - v_min*exp(2i*pi*nt)+...
        %             2i*pi*v_min*nt.*(expint(-2i*pi*nt*v_r)...
        %             -expint(-2i*pi*nt));
    end
end

