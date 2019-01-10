function [f_mod,t_chop,v_chop] = extract_moderator_f(in_mat,time_mod,V_mod,tau,L,tau_char,V_char,nip)
% in_mat -- 2D moderator function in units tau(mks) vs energy (mEv)
% time_mod - time axis at moderator  in mks
% lambda_mod - wavelength axis at moderator (in A)
% mod_Energy - energy at moderator array (in mEv)
% tau -- time at choper to shift function around (in chopper opening time
% units)
%
% Output:
% Interpolated fime-velocity profile at choper position
% t_chop -- time axis for the profile above (in units of chopper opening
%           time)
% v_chop -- velocity axis for the profile above (in m/s)

dt = 2*tau_char/(nip-1);

t_chop = (tau-tau_char:dt:tau+tau_char); % chopper opening time in sec



tau_mod = L./V_mod; %

tau_at_mod = bsxfun(@minus, t_chop, tau_mod'); % t-L/v;
%time_mod    = time_mod*(1.e-6/tau_char); % time at moderator in chopper opening times
tau_mod_min = min(time_mod);
tau_mod_max = max(time_mod);

non_zero = tau_at_mod>=tau_mod_min & tau_at_mod<=tau_mod_max;

nL = numel(V_mod);
V_mod_w   = reshape(repmat(V_mod,1,nip),nL,nip); % moderator velocity axis
v_chop    = V_mod_w(non_zero);
v_ch_min  = min(v_chop);
v_ch_max  = max(v_chop);
if v_ch_min == v_ch_max
    v_ch_min = v_ch_min-1;
    v_ch_max = v_ch_max+1;
    dV = 2/(0.5*nip-1);
else
    dV = (v_ch_max-v_ch_min)/(0.5*nip-1);
end
v_chop = v_ch_min:dV:v_ch_max;

[xb,yb] =meshgrid(time_mod,V_mod);
[xi,yi]= meshgrid(t_chop,v_chop);
xi = xi - L./yi;

f_mod = interp2(xb,yb,in_mat',xi,yi,'nearest',0);

fh = findobj('type','figure', 'Name', 'Moderator pulse');
if ~isempty(fh)
    figure(fh)
    hold on
    %[xi,yi]=meshgrid(tchop/tau_char,vchop/V_char);
    surf(xi/tau_char,yi/V_char,f_mod);
    view(0,90);
    hold off
end


%[tau1,tau2,t0,R,a0] = evaluate_tau(yi(:,1),f_mod(:,1));


function y = ikeda(x,pin)

a = pin(1);
b = pin(2);
R = pin(3);

t0 = pin(5);
xs = x-t0;


y = ((1-R)*(xs*a).^2.*exp(-xs*a)+...
    (R*a*a*b/(a-b)^3)*...
    (2*exp(-xs*b) - (2+2*(a-b)*xs + (a+b)^2*xs.*xs).*exp(-xs*a)));


function [aa,be,t0,R,a0] = evaluate_tau(tf,sf)


np = numel(tf);
nf = floor(0.8*np);
t_tail = tf(nf:np);
s_tail = log(sf(nf:np));
t0 = t_tail(1);
t_tail = t_tail-t0;
f_tail = polyfit(t_tail,s_tail,1);
be = -f_tail(1);
U  = exp(f_tail(2));
dt = tf(2)-tf(1);
Norma = sum(sf)*dt;
aa = 2/Norma*(1-2*U);

i1 = find(sf>1.e-6,1);
t_head = tf(i1:i1+7)-tf(i1);
s_head = sf(i1:i1+7);
f_head = polyfit(t_head,s_head,2);
C2 = f_head(1);
R = 1-C2/(2*aa*aa)-2*U;


t0=0;
a0=0;






%f_mod = interp2(yb',xb',in_mat,yi,xi,'nearest',0);

% convert velocity m/s to mEv
%vSq2mEv = 5.227e-6;  % s^2/m^2
%v_chop = v_chop.^2*vSq2mEv;
%
% Plot moderator
% persistent sf;
% if isempty(sf)
%     sf = figure;
%     contourf(xb,yb,in_mat');
% else
%     set(0, 'currentfigure', sf);
% end
%
% patch = ones(size(xi))*16000;
% hold('on');
% surf(xi,yi,patch);
% disp([xi(1,1),yi(1,1);xi(1,200),yi(1,200);xi(100,1),yi(100,1);xi(100,200),yi(100,200)])
%
%
%
%


