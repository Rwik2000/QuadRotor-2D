%%
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


%% 

clc;
clear;

syms y(t) z(t) phi(t)

g = 9.8;
m= 0.18; % mass in kg
I_xx = 2.5e-4; %moment of inertia in x dir
F2 = 1;
F4 = 1;
L = 0.086;

u1 = F2 + F4;
u2 = (F2 - F4)*L;

ode1 = diff(y, t, 2) == -u1*sin(phi)/m;
ode2 = diff(z, t, 2) == -g +u1*cos(phi)/m;
ode3 = diff(phi, t, 2) == u2;

odes = [ode1; ode2; ode3];

Dy = diff(y,t);
Dz = diff(z,t);
Dphi = diff(phi,t);

cond1 = [y(0)==0, Dy(0)==0];
cond2 = [z(0)==0, Dz(0)==0];
cond3 = [phi(0)==0, Dphi(0)==0];
conds = [cond1; cond2; cond3];

S = dsolve(odes, conds)
% [ySol(t), zSol(t), phiSol(t)] = dsolve(odes,conds)

% [ySol(t), zSol(t), phiSol(t)] = dsolve(odes)



