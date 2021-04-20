clc; 
clear;

syms y z phi y_dot z_dot phi_dot  
syms y_ddot z_ddot phi_ddot

g = 9.8;
m= 0.18; % mass in kg
I_xx = 2.5e-4; %moment of inertia in x dir
F2 = 0;
F4 = 0;
L = 0.086;

u1 = F2 + F4;
u2 = (F2 - F4)*L;

X = [y; z; phi; y_dot; z_dot; phi_dot];
X_dot = [y_dot; z_dot; phi_dot; 0; -g; 0] + [0 0; 0 0; 0 0; (-1/m)*sin(phi) 0; (1/m)*cos(phi) 0; 0 I_xx]*[u1; u2];


