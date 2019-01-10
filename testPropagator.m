function testPropagator
t=0:0.2:100;
v=10:1:100;
f = @(v,t)(exp(-(t-0.1*v).^2-(t-5).^2));

[xi,yi]= meshgrid(v,t);
f_val = f(xi,yi);
surf(xi,yi,f_val,'Edgecolor','None');

v0=5;
f2 = @(v,t)(exp(-(t-0.1*(v-v0)).^2-(t-5).^2));


f_val2 = f2(xi,yi);
surf(xi,yi,f_val2,'Edgecolor','None');