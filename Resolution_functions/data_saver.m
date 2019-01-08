classdef data_saver
    %Helper class to save or load resolution pulse data
    
    
    
    properties
        ds
    end
    
    methods
        function obj = data_saver(vel_distr,tsample,vsample,fsample,...
                V_pulseI,t_chop,t_det,f_det_vs_t,L_det,L_samp,tau_char,V_char)
            obj.ds = struct('tsample',tsample,'vsample',vsample,'fsample',fsample,...
                'V_pulseI',V_pulseI,'t_chop',t_chop,...
                't_det',t_det,'f_det_vs_t',f_det_vs_t,...
                'L_det',L_det,'L_samp',L_samp,...
                'tau_char',tau_char,'V_char',V_char,...
                'vel_distr_fun',vel_distr);
            if ~isa(obj.ds.vel_distr_fun,'function_handle')
                obj.ds.vel_distr_fun= @vel_distribution0;
            end
        end
        function save_data(obj)
            ds_ = obj.ds;
            fn = obj.pulse_name();
            save(fn,'-struct','ds_');
        end
        function obj = load_data(obj,filename)
            ds_ = load(filename);
            obj.ds = ds_;
        end
        function name = pulse_name(obj)
            % function to generate pulse dependent file name to cache appropriate data
            %
            velocity = obj.ds.V_pulseI;
            fh = obj.ds.vel_distr_fun;
            addinfo = fh('Name');
            
            name = sprintf('Pulse_V%d_%s_input_data',floor(velocity*100),addinfo);
        end
        function varargout = get_data(obj)
            fn = fieldnames(obj.ds);
            for i=1:nargout
                varargout{i} = obj.ds.(fn{i});
            end
        end
    end
    
end

