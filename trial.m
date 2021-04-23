%%
clc; 
clear;

syms y z phi y_dot z_dot phi_dot  
syms y_ddot z_ddot phi_ddot

g = 9.8;
m= 0.18; % mass in kg
I_xx = 2.5e-4; %moment of inertia in x dir
F2 = 1;
F4 = 1;
L = 0.086;

u1 = F2 + F4;
u2 = (F2 - F4)*L;

X = [y; z; phi; y_dot; z_dot; phi_dot];
X_dot = [y_dot; z_dot; phi_dot; 0; -g; 0] + [0 0; 0 0; 0 0; (-1/m)*sin(phi) 0; (1/m)*cos(phi) 0; 0 I_xx]*[u1; u2];


%% maa chuda

clc;
clear;

syms y(t) z(t) phi(t) T Y

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
ode3 = diff(phi, t, 2) == u2/I_xx;

odes = [ode1; ode2; ode3];

Dy = diff(y,t);
Dz = diff(z,t);
Dphi = diff(phi,t);

% cond1 = [y(0)==0, Dy(0)==0];
% cond2 = [z(0)==0, Dz(0)==0];
% cond3 = [phi(0)==0, Dphi(0)==0];
conds = [y(0)==0 Dy(0)==0 z(0)==0, Dz(0)==0 phi(0)==0, Dphi(0)==0];

dsolve(odes, conds);
% [VF, Sbs] = odeToVectorField(odes);
% DEFcn = matlabFunction(VF,'Vars',{T,Y});
% tspan = [0,10];
% [T,Y] = ode45(DEFcn, tspan, [0 0 0 0 0 0]);
% plot(T,Y)

%% savage (step response done)
% Non  Linear system of equation
clc;
clear;

syms y(t) z(t) phi(t) T Y

g = 9.8;
m= 0.18; % mass in kg
I_xx = 2.5e-4; %moment of inertia in x dir
F2 = 1;
F4 = 1;
L = 0.086;

u1 = F2 + F4;
u2 = (F2 - F4)*L;

f = @(t,x)[x(2);-u1*sin(x(5))/m;x(4);-g + u1*cos(x(5))/m; x(6);u2/I_xx];
[t, xa] = ode45(f, [0,10], [0 0 0 0 0 0]);

plot(t,xa(:, 1))
hold on;
plot(t,xa(:, 3))
hold on;
plot(t,xa(:, 5))
hold on;
legend('y','z','phi')
%% Linearizing the above system and finding the transfer function:
clc;
clear;

syms y(t) z(t) phi(t) T Y
syms I_xx g s m

% g = 9.8;
% m= 0.18; % mass in kg
% I_xx = 2.5e-4; %moment of inertia in x dir
F2 = 1;
F4 = 1;
L = 0.086;

% State space eqution : x' = Ax + Bu ; y = Cx + Du
A = [0 1 0 0 0 0; 
     0 0 0 0 -g 0; 
     0 0 0 1 0 0; 
     0 0 0 0 0 0; 
     0 0 0 0 0 1; 
     0 0 0 0 0 0];
 
B = [0 0; 
     0 0; 
     0 0; 
     0 1/m; 
     0 0; 
     0 1/I_xx;];

C = [1 0 0 0 0 0; 
      0 0 1 0 0 0; 
      0 0 0 0 1 0];
  
D =[0 0; 
     0 0; 
     0 0];%zero matrix shape: 3x2

% Phi=inv(s*eye(6)-A);
% H=C*Phi*B+D
% [num12, den12] = numden(H(1,2))
% C = coeffs(den12, s, 'All')
% r = roots(C)

 
% direct with numerical values
%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %sys = ss(A,B,C,D);
 %  sys_as_tf = tf(sys);
%  P = pole(sys_as_tf);
%  Z = tzero(sys_as_tf);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi=inv(s*eye(6)-A);
H=C*Phi*B+D

pole = [];
zero = [];
for i=1:2
    for j=1:3

        [num12, den12] = numden(H(j,i));
        Cd = coeffs(den12, s, 'All');
        Cn = coeffs(num12, s, 'All');
        rd = roots(Cd);
        rn = roots(Cn);
        pole = [pole; rd];
        zero = [zero; rn];
        
    end
end

pole
zero