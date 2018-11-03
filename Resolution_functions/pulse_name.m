function name = pulse_name(velocity,addinfo,Nt,Nv)
% function to generate pulse dependent file name to cache appropriate data
%
if ~exist('addinfo','var')
    addinfo = '';
end
if nargin == 2
    name = sprintf('Pulse_V%d_%s',floor(velocity*100),addinfo);
else
    name = sprintf('Pulse_V%d_Sz%dx%d_%s',floor(velocity*100),Nt,Nv,addinfo);
end

