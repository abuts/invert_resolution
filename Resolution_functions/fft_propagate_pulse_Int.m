function [f_out,t_out] = fft_propagate_pulse_Int(f_in,time_in,vel_in,L,t_char,v_char)
% Calculate interpolited time-velocity profile at the position L using fft
%
% f_in     -- 2D signal function in units tau(mks) vs
% time_in  -- time axis for signal  (sec)
% vel_in   -- velocity accis for signal (m/sec)
% L        -- distance to the target
% tau -- time at choper to shift function around (in chopper opening time
%
% Output:
%  Interpolated fime-velocity profile at position L
% f_out  --
% t_out -- time axis for the profile above (sec)
% v_out -- velocity axis for the profile above (in m/s)
v_min = min(vel_in);
v_max = max(vel_in);
DV0 = v_max - v_min;
tp_min = min(time_in);
tp_max = L/v_min+max(time_in);
DT0 = tp_max-tp_min;
dt = time_in(2) - time_in(1);
% if dt<2e-6 % debugging -- not needed in reality
%     dt = 2e-6;
%     t_int = min(time_in):dt:max(time_in);
%     [xi,yi]= meshgrid(t_int,vel_in);
%     [xb,yb] = meshgrid(time_in,vel_in);
%     f_in = interp2(xb,yb,f_in,xi,yi,'linear',0);
% end
t_out = tp_min:dt:tp_max;
Nt = numel(t_out);
Nv = numel(vel_in);


v_min_norm = v_min/DV0;
v_max_norm = v_max/DV0;

v_index = fft_ind(Nv);
t_index = fft_ind(Nt);

Imk = @(nv,kt)I_mk(nv,kt,v_min_norm,v_max_norm,v_index,t_index);


%I1 = Imk(5,5);
%
% [xi,yi]= meshgrid(t_out,vel_in);
% xi = xi-L./yi;
% [xb,yb] = meshgrid(time_in,vel_in);
% f_in_int = interp2(xb,yb,f_in,xi,yi,'linear',0);
%
%
% fft frequencies do not correspond to indexes as index m in fft array
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
v_phase =  exp(-1i*pi*v_index*((v_min+v_max)/(DV0)));

f_in_sp = bsxfun(@times,f_in_sp,v_phase');
%f_in_sp = bsxfun(@times,f_in_sp,t_phase);
f_in_t  = ifft2(f_in_sp);
figure(20)
[xi,yi] = meshgrid(t_out/t_char,vel_in/v_char);
surf(xi,yi,abs(f_in_t),'Edgecolor','None');
view(0,90);

persistent f_con_mat;
if isempty(f_con_mat)
    cn = sprintf('CashNv%dNt%d',Nv,Nt);
    if exist([cn,'.mat'],'file')
        disp(['***** loading ',cn]);
        cns = load(cn);
        f_con_mat = cns.f_con_mat;
    else
        disp(['***** processing ',cn]);
        f_con_mat = ones(Nv,Nt)*NaN;
    end
end

t_start=tic;
if any(isnan(reshape(f_con_mat,1,numel(f_con_mat))))
    call_count = 0;
    for m=1:Nv
        undef = isnan(f_con_mat(m,:));
        if ~any(undef)
            continue;
        end
        for n=1:Nt
            if undef(n)
                f_con_mat(m,n) = Imk(m,n);
            end
        end
        if call_count > 10
            fprintf('k=%d#%d\n',m,Nv);
            cn = sprintf('CashNv%dNt%d',Nv,Nt);
            save(cn,'f_con_mat');
            call_count = 0;
        end
        call_count = call_count+1;
        
    end
    fprintf('k=%d#%d\n',m,Nv);
    cn = sprintf('CashNv%dNt%d',Nv,Nt);
    save(cn,'f_con_mat');
end

f_in_sp = f_in_sp.*f_con_mat;
%v_phase =  exp(1i*pi*v_index);

% f_in_sp = bsxfun(@times,f_in_sp,v_phase');
f_out_sp = sum(f_in_sp.*f_con_mat,1);
t_phase = (DV0)*exp(-1i*pi*t_index*((tp_min+tp_max)/(DT0)-1)); %
%t_phase = (DV0)*exp(1i*pi*t_index); %
f_out_sp = f_out_sp.*t_phase;

time_c = toc(t_start)/60 % convert in minutes
f_out = ifft(f_out_sp);



