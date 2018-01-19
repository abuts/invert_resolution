function  test_fft_propagate_pulse()

%Nt =295;
%Nv = 220;
Nt =1066;
Nv = 636;

Nvs = 100;
Nts = 3;

iv = 1:Nv;
v_min = 1;
v_max = 10;
dv =  (v_max-v_min)/(Nv-1);
v_char = 1;
vel_in = v_min+(iv-1)*dv;
L=10;


Nt_in = 100; % number of points in initial resolution 
dt=(L/v_min-L/v_max)/(Nt-Nt_in);
t0=10;
t_char = 1;
it = 1:Nt_in;
time_in = t0+(it-1)*dt;


v_s = v_min +(Nvs-1)*dv;
t_s = t0+(Nts-1)*dt;
ta = t_s+L/v_s;
V_pulse = v_s ; % not used

f_in = zeros(Nv,Nt_in);
f_in(Nvs ,Nts) = Nv*Nt;



[f_out,t_out] = fft_propagate_pulse_IntW(f_in,time_in,vel_in,L,V_pulse,t_char,v_char);

disp(['arrival time: ',num2str(ta)]);
p1 = IX_dataset_1d(t_out/t_char,abs(f_out));
p1.x_axis = sprintf('Time/(%3.2g sec)',t_char);
p1.s_axis = 'Signal';
dl(p1);
keep_figure
p1.signal = imag(f_out);
p1.s_axis = 'Img error';
dl(p1);


end

