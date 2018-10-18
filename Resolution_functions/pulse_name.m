function name = pulse_name(velocity,addinfo)
% function to generate pulse dependent file name to cache appropriate data
%
if ~exist('addinfo','var')
    addinfo = '';
end
name = sprintf('Pulse_V%d_%s',floor(velocity*100),addinfo);


