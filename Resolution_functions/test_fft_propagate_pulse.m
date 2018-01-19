function  test_fft_propagate_pulse()

Nt =295;
Nv = 220;

iv = 1:Nv;
v_min = 1;
dv = 10;
v_char = dv;
vel_in = v_min+(iv-1)*dv;
v_max = max(vel_in);
L=1;

Nt_in = 10;
dt=(L/v_min-L/v_max)/(Nt-Nt_in);
t0=0;
t_char = 1;
it = 1:Nt_in;
time_in = t0+(it-1)*dt;
V_pulse = 0; % not used

f_in = zeros(Nv,Nt_in);
f_in(180,2) = Nv*Nt;



[f_out,t_out] = fft_propagate_pulse_IntW(f_in,time_in,vel_in,L,V_pulse,t_char,v_char);

p1 = IX_dataset_1d(t_out/t_char,abs(f_out));
p1.x_axis = sprintf('Time/(%3.2g sec)',t_char);
p1.s_axis = 'Signal';
dl(p1);
keep_figure
p1.signal = imag(f_out);
p1.s_axis = 'Img error';
dl(p1);


end

