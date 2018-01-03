function [f_out,dv_out] = fft_invert_propagation(f_samp,t_samp,v_samp,f_det_vs_t,t_det,L_det,v_max)
% retrieve samle gain/loss velocity distribution given signal on the detector and
% and the smple incident beam time/velocity distribution.

t0 = min(t_samp);
t_det = t_det-t0;

dt = t_det(2)-t_det(1);
if dt<2e-6 % debugging -- not needed in reality
    dt = 2e-6;
end
tp_min = 0;
tp_max = max(t_det);
[t_in,dt,Nt] = adjust_step(tp_min,tp_max,dt);

f_det = interp1(t_det,f_det_vs_t,t_in,'linear',0);
if ~exist('v_max','var')
    i1 = find(f_det,1);
    v_max = L_det/t_in(i1);
end
v_min = L_det/(tp_max-tp_min);

dv = v_samp(2)-v_samp(1);
[dv_out,dv,Nv] = adjust_step(-v_max,v_max,dv);
[xb,yb] = meshgrid(t_samp-t0,v_samp);
[xi,yi] = meshgrid(t_in,dv_out);
f_in = interp2(xb,yb,f_samp,xi,yi,'linear',0); % expand sample pulse onto the whole integration range



% %
% % % [xi,yi]= meshgrid(t_out,vel_in);
% % [xb,yb] = meshgrid(t_int,vel_in);
% surf(xi/t_char,yi/v_char,f_in,'Edgecolor','None');
% view(0,90);

% f_in_int = interp2(xb,yb,f_in,xi,yi,'linear',0);
%
%
% fft frequencies do not correspond to indexes as index m in fft array
% has frequency m-1 if m<N/2 and end-m if bigger;
f_in_sp = fft2(f_in,Nv,Nt);
% v_phase =  exp(-1i*pi*v_index*DNv/(Nv-1)); % appears to shift zeros added at the end to zeros, added to the beginning
% %
% f_in_sp = bsxfun(@times,f_in_sp,v_phase');
%f_in_sp = bsxfun(@times,f_in_sp,t_phase);
% f_in_t  = ifft2(f_in_sp);
% figure(21)
% [xi,yi] = meshgrid(t_in,dv_out);
% surf(xi,yi,abs(f_in_t),'Edgecolor','None');
% view(0,90);



v_index = fft_ind(Nv);
t_index = fft_ind(Nt);
DV0 = max(dv_out)-min(dv_out);



f_con_mat = [];
if isempty(f_con_mat)
    cn = sprintf('CashNv%dNt%d',Nv,Nt);
    if exist([cn,'.mat'],'file')
        disp(['***** loading ',cn]);
        cns = load(cn);
        f_con_mat = cns.f_con_mat;
    else
        disp(['***** processing ',cn]);
        f_arr = cell(Nv,1);
        f_con_mat = ones(Nv,Nt)*NaN;
    end
end

t_start=tic;
if any(isnan(reshape(f_con_mat,1,numel(f_con_mat))))
    v_min_norm = v_min/(DV0);
    v_max_norm = v_max/(DV0);
    Imk = @(nv,kt)I_mk(nv,kt,v_min_norm,v_max_norm,v_index,t_index);
    
    for m=1:Nv
        undef = isnan(f_con_mat(m,:));
        if ~any(undef)
            continue;
        end
        
        vec_mat = f_con_mat(m,:);
        for n=1:Nt
            if undef(n)
                vec_mat(n) = DV0*Imk(m,n);
            end
        end
        f_arr{m} = vec_mat;
        if rem(m,10) == 0
            fprintf('k=%d#%d\n',m,Nv);
            %             cn = sprintf('CashNv%dNt%d',Nv,Nt);
            %             save(cn,'f_con_mat');
        end
    end
    f_con_mat = reshape([f_arr{:}],Nt,Nv)';
    fprintf('final write of %d\n',Nv);
    cn = sprintf('CashNv%dNt%d',Nv,Nt);
    save(cn,'f_con_mat');
end

f_in_sp = f_in_sp.*f_con_mat;
f_det_sp = fft(f_det,Nt);
%
f_vel_sp  = linsolve(f_in_sp,f_det_sp);


%v_phase =  exp(1i*pi*v_index);
% f_in_sp = bsxfun(@times,f_in_sp,v_phase');

%t_phase = (DV0)*exp(-1i*pi*t_index*((tp_min+tp_max)/(2*tp_max))); %
%t_phase = (DV0)*exp(1i*pi*t_index); %
%f_vel_sp   = f_vel_sp .*v_phase;
%
time_c = toc(t_start)/60; % convert in minutes
fprintf(' Calculations take: %f min\n',time_c);
f_out = ifft(f_vel_sp);



