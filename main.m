clear
%% 时间步长设置
dt=0.1;
t1=10;
t2=100;
t3=215;
n=t3/dt;
%% 常量定义
m_0=23000;
dm=100;
J_2=1.08263e-3;
a_0_phi=0.5;
a_0_psi=0.5;
P_0=3e5;
m_z_alpha=-0.12*180/pi;
omega_e=7.292e-5;
L=15;
m_zCy=0.1;
C_y_alpha=-m_z_alpha/m_zCy;
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
R_0=earthR_0(A_0,rad2deg(phi_0),2);
B_0=atan(a_e^2/b_e^2 * tan(phi_0));
omega_e1=omega_e*[cos(B_0)*cos(A_0);sin(B_0);-cos(B_0)*sin(A_0)];

%% 变量初始化
theta=zeros(n,1);
V=zeros(3,n);
r=zeros(3,n);
m=zeros(n,1);
alpha=zeros(n,1);
beta=zeros(n,1);
sigma=zeros(n,1);
phi=zeros(n,1);
psi=zeros(n,1);
delta_phi=zeros(n,1);
delta_psi=zeros(n,1);
h=zeros(n,1);
Phi=zeros(n,1);
R=zeros(n,1);
p_H=zeros(n,1);
rho=zeros(n,1);
phi_pr=zeros(n,1);
t=zeros(n,1);
delta_v_1k=zeros(n,1);

theta(1)=pi/2;
V(:,1)=[0;0.1;0];
r(:,1)=[0;0;0];
m(1)=23000;
phi(1)=pi/2;
R(1)=norm(R_0);
[~,~,p_H(1),rho(1)]=atmosisa(h(1));



%% 程序角配置
for i=1:n
    t(i)=i*dt;
    if t(i)<t1
        phi_pr(i)=pi/2;
    elseif t(i)<t2
        phi_pr(i)=pi/2+(pi/2-deg2rad(66))*(((t(i)-10)/90)^2-2*(t(i)-10)/90);
    else
        phi_pr(i)=phi_pr(t2/dt-1)/(t2-t3)*(t(i)-t3);
    end
end
%% 主循环
for i=1:n-1
    q=0.5*rho(i)*(V(:,i)'*V(:,i));
    C_x=0.02+0.005*(rad2deg(alpha(i)))^2;
    G_E=goaround(-pi/2-A_0,2)*goaround(phi_0,1)*goaround(lambda_0-pi/2,3);
    G_B=goaround(psi(i),2)*goaround(phi(i),3);
    V_G=goaround(sigma(i),2)*goaround(theta(i),3);
    B_V=goaround(alpha(i),2)*goaround(beta(i),3);
    
    

    P=[P_0+S_e*(p_0-p_H(i));0;0];%推力
    R1=[-C_x*q*S_M;C_y_alpha*q*S_M*alpha(i);-C_y_alpha*q*S_M*beta(i)];%气动力
    F_c=[0;1/sqrt(2)*delta_phi(i);-1/sqrt(2)*delta_psi(i)]*norm(P);%控制力
    r1=r(:,i)+R_0;
    g1_r=-fmM/(r1'*r1)*(1+1.5*J_2*(a_e/norm(r1))^2*(1-5*sin(Phi(i))^2));
    gwe=-3*fmM/(r1'*r1)*J_2*(a_e/norm(r1))^2*sin(Phi(i));
    g=g1_r*r1/norm(r1)+gwe*omega_e1/omega_e;%引力加速度
    F_e=-m(i)*cross(omega_e1,cross(omega_e1,r(:,i)+R_0));%离心惯性力
    F_k=-2*m(i)*cross(omega_e1,V(:,i));%哥氏惯性力


    M_z1_alpha=m_z_alpha*q*S_M*L;
    M_y1_beta=-m_z_alpha*q*S_M;%气动力矩系数

    M_z1_delta=-sqrt(0.5)*norm(P)*rho_e;
    M_y1_delta=sqrt(0.5)*norm(P)*rho_e;%控制力矩系数

    A_phi=a_0_phi*M_z1_delta/(M_z1_alpha+a_0_phi*M_z1_delta);
    A_psi=a_0_phi*M_y1_delta/(M_y1_beta+a_0_phi*M_y1_delta);
    
    dV=1/m(i) * (G_B'*(P+F_c)+V_G'*R1+F_e+F_k)+g;
    V(:,i+1)=V(:,i)+dt*dV;
    r(:,i+1)=r(:,i)+dt*V(:,i+1);

    theta(i+1)=asin(V(2,i+1)/norm(V(:,i+1))); 
    sigma(i+1)=-asin(V(3,i+1)/norm(V(:,i+1)));

    alpha(i+1)=A_phi*(phi_pr(i+1)-omega_e1(3)*t(i+1)-theta(i+1));
    phi(i+1)=theta(i+1)+alpha(i+1);
    beta(i+1)=A_psi*((omega_e1(1)*sin(phi(i+1))-omega_e1(2)*cos(phi(i+1)))*t(i+1)-sigma(i+1));
    psi(i+1)=sigma(i+1)+beta(i+1);

    delta_phi(i+1)=-a_0_phi*(phi_pr(i+1)-omega_e1(3)*t(i+1)-theta(i+1));
    delta_psi(i+1)=a_0_psi*(psi(i+1)+(omega_e1(1)*sin(phi(i+1))-omega_e1(2)*cos(phi(i+1)))*t(i+1));

    m(i+1)=m_0-dm*t(i+1);
    Phi(i+1)=asin((r(:,i+1)+R_0)'*omega_e1/(norm((r(:,i+1)+R_0))*omega_e));
    R(i+1)=a_e*b_e/sqrt(b_e^2*cos(Phi(i+1))^2+a_e^2*sin(Phi(i+1))^2);
    h(i+1)=norm(r(:,i+1)+R_0)-R(i+1);
    [~,~,p_H(i+1),rho(i+1)]=atmosisa(h(i+1));
end




%% 数据后处理
u_r=2701.325;

plot(t,[phi,phi_pr,theta,delta_phi,alpha]);