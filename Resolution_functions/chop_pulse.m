function chop_shape = chop_pulse(velocity,time,t_opening,R_chop)
% Calculate chopper pulse shape as function of time and neutrons velocity 
% at chopper position
%
% Additional parameters are:
% t_opening -- chopper opening time.
% R_chop    -- chopper radius
%Time has to be units of chopper opening time e.g.
% units of Dt=H_chopper/(Omega*R_Chopper) where omega -- circular rotation
%       velocity
% velocity: in units of m/s*Dt above (so [m] left)
% [R_chop] = m
% 
% Implpementing forumala derived in:
% NI&M, v4 issue 3 pp140-150 April 1959
%

n_times = numel(time);
n_vel  = numel(velocity);
chop_shape = zeros(n_times,n_vel );
%

mult2 = bsxfun(@plus, 2*R_chop./velocity, abs(time-t_opening)'); % 
ff = 1 - bsxfun(@times,(velocity/(4*R_chop)),mult2.^2);
out = ff<0;
ff(out) = 0;

Gm = 4*R_chop./velocity;
Gm_lt_1  = reshape(repmat(Gm<1,n_times,1),n_times,n_vel);
ff_range = reshape(repmat(Gm<=4 & Gm>=1,n_times,1),n_times,n_vel);

chop_shape(ff_range) = ff(ff_range);
%chop_shape = ff;


f1 = reshape(repmat(1 - 2*abs(time - t_opening),1,n_vel)',n_times,n_vel);
f1_gt_0 = f1>0;
f2_range =Gm_lt_1&f1_gt_0; 
%f2_range =f1_gt_0; 


f1_gt_fv =bsxfun(@gt,f1,1-Gm);
f1_range= f2_range &  ~f1_gt_fv;

chop_shape(f1_range) = f1(f1_range);
ff_range= f2_range &  f1_gt_fv;
chop_shape(ff_range) = ff(ff_range);
chop_shape = chop_shape';



