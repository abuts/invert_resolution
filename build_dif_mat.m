function dif_mat  = build_dif_mat(betta0,v_min_norm,v_max_norm,funName,Nv,Nt)
%   Calculate diffraction matrix describing the propagation of a momentum
%   from sample to a detector
v_index = fft_ind(Nv);
t_index = fft_ind(Nt);

res_name = sprintf('Cash_%s_Nv%dNt%d',funName,Nv,Nt);
t_start=tic;


%persistent f_con_mat;
dif_mat = [];
if isempty(dif_mat)
    if exist([res_name,'.mat'],'file')
        disp(['***** loading ',res_name]);
        cns = load(res_name);
        if isfield(cns,'f_con_mat')
            dif_mat = cns.f_con_mat;
        else
            dif_mat = cns.dif_mat;
        end
    else
        disp(['***** processing ',res_name]);
        dif_mat = ones(Nv,Nt)*NaN;
    end
end
if any(isnan(reshape(dif_mat,1,numel(dif_mat))))
    %     ws = warning('off','MATLAB:integral:MaxIntervalCountReached');
    %     clob = onCleanup(@()warning(ws));
    %I1 = Imk(1,2);
    for n=1:Nt
        undef = isnan(dif_mat(:,n));
        if ~any(undef)
            continue;
        end
        
        vec_mat = dif_mat(:,n);
        parfor m=1:Nv
            if undef(m)
                vec_mat(m) = I_mk4(m,n,betta0,v_min_norm,v_max_norm,v_index,t_index);
            end
        end
        %         nt = t_index(n);
        %         if nt ~= 0
        %             vec_mat = int_fitter(vec_mat,v_steps_norm,nt,v_min_norm,v_max_norm);
        %         end
        %
        dif_mat(:,n)= vec_mat;
        if rem(n,5) == 0
            fprintf('writing step %d#%d\n',n,Nt);
            save(res_name,'dif_mat');
        end
    end
    %f_con_mat = reshape([f_arr{:}],Nt,Nv)';
    fprintf('final write of %d\n',Nv);
    save(res_name,'dif_mat');
end
time_c = toc(t_start)/60; % convert in minutes
fprintf(' Calculations take: %f min\n',time_c);
