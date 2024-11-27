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
C_y_alpha=m_z_alpha/m_zCy;
rho_e=7;
x_g=8;
S_e=0.1;
S_M=3.14;
phi_0=deg2rad(19+19/60);
lambda_0=deg2rad(109+48/60);
A_0=deg2rad(30);
a_e=6378137;
b_e=6356752;
p_0=101325;
fmM=3.986005e14;
R_0=earthR_0(A_0,rad2deg(Phi),2);
B_0=atan(a_e^2/b_e^2 * tan(Phi));
omega_e1=omega_e*[cos(B_0)*cos(A_0);sin(B_0);-cos(B_0)*sin(A_0)];
%% 变量初始化
theta=zeros(2150,1);
V=zeros(3,2150);
r=zeros(3,2150);
m=zeros(2150,1);
alpha=zeros(2150,1);
beta=zeros(2150,1);
sigma=zeros(2150,1);
phi=zeros(2150,1);
psi=zeros(2150,1);
delta_phi=zeros(2150,1);
delta_psi=zeros(2150,1);
h=zeros(2150,1);
Phi=zeros(2150,1);
R=zeros(2150,1);
p_H=zeros(2150,1);
rho=zeros(2150,1);


theta(1)=pi/2;
V(:,1)=[0;0.1;0];
r(:,1)=[0;0;0];
m(1)=23000;
phi(1)=90;
R(1)=norm(R_0);
[~,~,p_H(1),rho(1)]=atmosisa(h);
%% 主循环
dt=0.1;%时间步长
t1=10;
t2=100;
t3=215;

for i=1:2149
    q=0.5*rho*(V'*V);
    C_x=0.02+0.005*(rad2deg(alpha(i)))^2;
    G_B=(goaround(psi(i),2)*goaround(phi(i),3));

    P=[P_0+S_e*(p_0-p_H(i));0;0];%推力
    R1=[-C_x*q*S_M;C_y_alpha*q*S_M*alpha(i);-C_y_alpha*q*S_M*beta(i)];%气动力
    F_c=[0;1/sqrt(2)*delta_phi(i);-1/sqrt(2)*delta_psi(i)]*norm(P);%控制力
    g1_r=-fmM/(r(:,i)'*r(:,i))*(1+1.5*J_2*(a_e/norm(r(:,i)))^2*(1-5*sin(phi(i))^2));
    gwe=-3*fmM/(r(:,i)'*r(:,i))*J_2*(a_e/norm(r(:,i)))^2*sin(phi(i));
    g=g1_r*(R_0+r(:,i))/norm(r(:,i))+gwe*omega_e1/omega_e;%引力加速度
    F_e=-m(i)*cross(omega_e1,cross(omega_e1,r(:,i)+R_0));%离心惯性力
    F_k=-2*m(i)*cross(omega_e1,V(:,i));%哥氏惯性力

    

end
