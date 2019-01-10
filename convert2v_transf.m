function v_scale = convert2v_transf(time_scale,L_det,L_sam,t0,v0)
% convert the arrival time into velocity transfer. The function similar to
% convert_to_energy_transfer used in vanilla reduction
%
% time_scale -- the neutron arival time to the detector
% L_det -- distance from sample to detector
% L_sam -- distance from source to sample 
% t0 -- time of the pulse peak at moderator
% v0 -- average velocity in the input burst
%
t_det = time_scale-t0;

t_sample = L_sam/v0;
v_scale = L_det./(t_det-t_sample)-v0;
