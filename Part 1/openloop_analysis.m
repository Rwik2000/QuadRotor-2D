%%
clc;
clear;

clc;
clear;

syms y(t) z(t) phi(t) T Y
% syms I_xx g s m

g = 9.8;
m= 0.18; % mass in kg
I_xx = 2.5e-4; %moment of inertia in x dir
l = 0.086;

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
     1/m 1/m; 
     0 0; 
     l/I_xx -l/I_xx;];

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
%  sys = ss(A,B,C,D);
%  sys_as_tf = tf(sys);
% step(sys_as_tf);
% [y_OL,t]=step(sys_as_tf);
% step(sys_as_tf);
% xlim([0 10])
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sys = ss(A,B,C,D);
sys_as_tf = tf(sys);
t = 0:0.04:10;
F1 = ones(length(t),1);
F2 = ones(length(t),1);

% F1 = 1;
% F2 = 1;
lsim(sys,[F1 F2],t);

