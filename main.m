%% 常量定义
m_0=23000;
dm=100;
J_2=1.08263e-3;
a_0_phi=0.5;
a_0_psi=0.5;
P_0=3e5;
m_z_alpha=-0.12*180/pi;
omega_e=7.29e5;
L=15;
m_zCy=0.1;
rho_e=7;
x_g=8;
S_e=0.1;
S_M=3.14;
phi_0=deg2rad(19+19/60);
lambda_0=deg2rad(109+48/60);
A_0=deg2rad(30);
a_e=6378145;
b_e=6356752;
p_0=101325;
fmM=3.986005e14;

%% 变量初始化
theta=pi/2;
V=[0;0.1;0];
r=[0;0;0];
m=23000;
alpha=0;
beta=0;
sigma=0;
phi=90;
psi=0;
delta_phi=0;
delta_psi=0;
h=0;
Phi=phi_0;
p_H=p_0;