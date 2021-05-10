clc;
clear;

kp = [1, 1,1];
kd = [1, 100, 1];
y_init = 0;
z_init = 0;
phi_init = 0;
ydot = 0;
zdot = 0;
phidot = 0;

del_t = 0.01;
g = 9.8;
L = 0.086;
I_xx = 2.5e-4;
m = 0.18;

y_des = 7;
z_des = 0;

%%%%%%%%
phi_des = 0;

t=0:del_t:10;

y = zeros(length(t),1);
z = zeros(length(t),1);
phi = zeros(length(t),1);
prev_phi_c = 0;

z(1) = z_init;
y(1) = y_init;
phi(1) = phi_init;

rise_time_y = 0;
rise_time_z = 0;
check_y = 0;
check_z = 0;

for i = 2:length(t)-1
    del_y = y(i)-y(i-1);
    del_z = z(i)-z(i-1);
    del_phi = phi(i)-phi(i-1);
    
    u1 = m*(g + 0 +kd(2)*(2 - del_z/del_t) + kp(2)*(2*t(i) - z(i)));
    phi_c = -1/g*(0 + kd(1)*(1 - del_y/del_t) + kp(1)*(t(i) - y(i)));
    del_phidot = (phi_c - prev_phi_c)/del_t - (phi(i)-phi(i-1))/del_t;
    u2 = I_xx*(0 + kd(3)*(del_phidot) + kp(3)*(phi_c - phi(i)));
    
    yd_dot = -g*phi(i);
    zd_dot = -g +u1/m;
    phid_dot = u2/I_xx;
    
    y(i+1) = y(i) + ydot*del_t + 0.5*yd_dot*(del_t^2);
    ydot = ydot + yd_dot*del_t;
    
    z(i+1) = z(i) + zdot*del_t + 0.5*zd_dot*(del_t^2);
    zdot = zdot + zd_dot*del_t;
    
    phi(i+1) = phi(i) + phidot*del_t + 0.5*phid_dot*(del_t^2);
    phidot = phidot + phid_dot*del_t;
    prev_phi_c = phi_c;
    
    if phi(i+1)>2*pi
        phi(i+1) = rem(phi(i+1),2*pi);
    end
    
    if z(i+1) >  0.1*(z_des - z_init) && check_z == 0
        rise_time_z = t(i+1);
        check_z = 1;
    elseif z(i+1) >  0.9*(z_des - z_init) && check_z == 1
        rise_time_z = t(i+1)-rise_time_z;
        check_z = 2;
    end
    
    if y(i+1) >  0.1*(y_des - y_init) && check_y == 0
        rise_time_y = t(i+1);
        check_y = 1;
    elseif y(i+1) >  0.9*(y_des - y_init) && check_y == 1
        rise_time_y = t(i+1)-rise_time_y;
        check_y = 2;
    end
end
per_OS_z = ((max(z)- z_des)/z_des)*100
rise_time_z

per_OS_y = ((max(y)- y_des)/y_des)*100
rise_time_y
figure(1)
plot(t, y)
hold on;
plot(t, z)
hold on;
plot(t, phi)
legend('y','z','phi')
grid on