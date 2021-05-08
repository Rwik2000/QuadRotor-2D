clc;
clear;

kp = [20, 80, 1000];
kd = [5, 20, 10];
ki = [0, 0, 0];
% kphi = [0, 0, 0];
% kp = [kp_y, kp_z, kp_phi];
% kd = [kd_y, kd_z, kd_phi];
% kphi = [kphi_y, kphi_z, kphi_phi];

del_t = 0.1;
g = 9.8;
L = 0.086;
I_xx = 2.5e-4;
m = 0.18;

y_init = 0;
z_init = 0;
phi_init = 0;

ydot = 0;
zdot = 0;
phidot = 0;

y_des = 7;
z_des = 12;
phi_des = 0;

t=0:del_t:100;

y = zeros(length(t),1);
z = zeros(length(t),1);
phi = zeros(length(t),1);

z(1) = z_init;
y(1) = y_init;
phi(1) = phi_init;

length(t);

for i = 1:length(t)-1
    del_z = z_des - z(i); 
    del_y = y_des - y(i);
    del_phi = phi_des - phi(i);
    
    u1 = m*(g+(kp(2)+kd(2)/del_t + ki(2)*del_t)*del_z);
    u2 =(kp(3)+kd(3)/del_t  + ki(3)*del_t)+del_phi*I_xx;
    phi_des = -1/g*(kp(1)+kd(1)/del_t + ki(1)*del_t)*del_y;
    
    yd_dot = -g*phi(i);
    zd_dot = -g +u1/m;
    phid_dot = u2/I_xx;
    
    y(i+1) = y(i) + ydot*del_t + 0.5*yd_dot*(del_t^2);
    ydot = ydot + yd_dot*del_t;
    
    z(i+1) = z(i) + zdot*del_t + 0.5*zd_dot*(del_t^2);
    zdot = zdot + zd_dot*del_t;
    
    phi(i+1) = phi(i) + phidot*del_t + 0.5*phid_dot*(del_t^2);
    phidot = phidot + phid_dot*del_t;
    
      
end

figure(1)
plot(t, y)
hold on;
plot(t, z)
hold on;
plot(t, phi)
legend('y','z','phi')

figure(2)
plot(y,z)
grid on