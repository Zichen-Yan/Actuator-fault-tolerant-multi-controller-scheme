clc;close all;clear all
%% 初始参数设置
A=[0.9931 0.0035;0.0068 0.9823];
B=[0.0081 -0.0032 -0.0034;0.0000 0.0032 0.0034];
E=-[0.9966;0.0034]*10^(-3);
w=Polyhedron('lb',-10^(-3),'ub',10^(-3));

F2_s1=diag([1 0 1]);
F2_s2=diag([1 0.6 1]);
F3_s1=diag([1 1 0]);
F3_s2=diag([1 1 0.6]);
F4_star=eye(3);

F2=diag([1 0 1]);
F3=diag([1 1 0]);
F4=eye(3);

C=eye(2);
eta=Polyhedron('lb',[-10^(-5), -10^(-5)],'ub',[10^(-5), 10^(-5)]);
Cv=[0 1];
H=Cv;
L = (place(A',C',[0.1,0.05]))';
K2=[74.0487 144.0117 41.4957;0 0 0;1.3575 350.8722 103.8830];
K21=K2(:,1:2);
K22=K2(:,3);
K3=[49.4359 102.8133 20.8364;1.8759 246.6231 47.2899;0 0 0];
K31=K3(:,1:2);
K32=K3(:,3);
K4=[74.0487 144.0117 41.4957;0.6775 175.1142 51.8462;0.7199 186.0589 55.0866];
K41=K4(:,1:2);
K42=K4(:,3);

d_ref=[0;-0.2185;0.2056];
u_ref2=[0.0356;0;0.1053];
u_ref3=[0.0356;0.1119;0];
u_ref4=[0.0356;0.0526;0.0558]+d_ref;

%% 构建初始不变集X
v=[10^(-3);10^(-5);10^(-5)];
mid1=A-L*C;
mid2=[E -L];
X_init=construct(mid1,mid2,v);
Epsilon_ll=C*X_init+eta;

%%
A_F4_K4=[A zeros(2,1);Cv 1]-[B*F4;zeros(1,3)]*K4;
A_F3_K3=[A zeros(2,1);Cv 1]-[B*F3;zeros(1,3)]*K3;
A_F2_K2=[A zeros(2,1);Cv 1]-[B*F2;zeros(1,3)]*K2;

%% 构造不变集zeta2 3 4
mid1=max(X_init.V(:,1));
mid2=max(X_init.V(:,2));
v=[mid1;mid2;10^(-3);10^(-5);10^(-5)];
mid1=A_F4_K4;
mid2=[B*F4*K41 E zeros(2,2);zeros(1,3) H];
Zeta4=construct(mid1,mid2,v);

mid1=A_F3_K3;
mid2=[B*F3*K31 E zeros(2,2);zeros(1,3) H];
Zeta3=construct(mid1,mid2,v);

mid1=A_F2_K2;
mid2=[B*F2*K21 E zeros(2,2);zeros(1,3) H];
Zeta2=construct(mid1,mid2,v);
%% 发生故障时的不变集 2->3 3->2
Epsilon3_f2=C*B*(F3-F2_s1*F3)*K3*Zeta3+C*(A-L*C+B*(F2_s1*F3-F3)*K31)*X_init+C*B*(F2_s1*F3-F3)*u_ref3+C*E*w+(-C*L)*eta+eta;
mid=C*B*(F3-F2_s2*F3)*K3*Zeta3+C*(A-L*C+B*(F2_s2*F3-F3)*K31)*X_init+C*B*(F2_s2*F3-F3)*u_ref3+C*E*w+(-C*L)*eta+eta;
P(1)=mid;
P(2)=Epsilon3_f2;
Epsilon3_f2=PolyUnion(P).convexHull;

Epsilon2_f3=C*B*(F2-F3_s1*F2)*K2*Zeta2+C*(A-L*C+B*(F3_s1*F2-F2)*K21)*X_init+C*B*(F3_s1*F2-F2)*u_ref2+C*E*w+(-C*L)*eta+eta;
mid=C*B*(F2-F3_s2*F2)*K2*Zeta2+C*(A-L*C+B*(F3_s2*F2-F2)*K21)*X_init+C*B*(F3_s2*F2-F2)*u_ref2+C*E*w+(-C*L)*eta+eta;
M(1)=mid;
M(2)=Epsilon2_f3;
Epsilon2_f3=PolyUnion(M).convexHull;

% figure()
% plot(Epsilon3_f2)
% hold on
% plot(Epsilon_ll)
% figure()
% plot(Epsilon2_f3)
% hold on
% plot(Epsilon_ll)

%%  44->2、3
Epsilon4_f2=C*B*(F4-F2_s1*F4)*K4*Zeta4+C*(A-L*C+B*(F2_s1*F4-F4)*K41)*X_init+C*B*(F2_s1*F4-F4)*u_ref4+C*E*w+(-C*L)*eta+eta;
mid=C*B*(F4-F2_s2*F4)*K4*Zeta4+C*(A-L*C+B*(F2_s2*F4-F4)*K41)*X_init+C*B*(F2_s2*F4-F4)*u_ref4+C*E*w+(-C*L)*eta+eta;
P(1)=mid;
P(2)=Epsilon4_f2;
Epsilon4_f2=PolyUnion(P).convexHull;

Epsilon4_f3=C*B*(F4-F3_s1*F4)*K4*Zeta4+C*(A-L*C+B*(F3_s1*F4-F4)*K41)*X_init+C*B*(F3_s1*F4-F4)*u_ref4+C*E*w+(-C*L)*eta+eta;
mid=C*B*(F4-F3_s2*F4)*K4*Zeta4+C*(A-L*C+B*(F3_s2*F4-F4)*K41)*X_init+C*B*(F3_s2*F4-F4)*u_ref4+C*E*w+(-C*L)*eta+eta;
P(1)=mid;
P(2)=Epsilon4_f3;
Epsilon4_f3=PolyUnion(P).convexHull;

% figure()
% plot(Epsilon4_f2)
% hold on
% plot(Epsilon4_f3)
% hold on
% plot(Epsilon_ll)

%% 仿真
iter=2000;

x=zeros(2,iter+1);
y=zeros(2,iter+1);
x_ref=zeros(2,iter+1);
y_ref=zeros(2,iter+1);
x_obs2=zeros(2,iter+1); %2
y_obs2=zeros(2,iter+1);
x_obs3=zeros(2,iter+1); %3
y_obs3=zeros(2,iter+1);
x_obs4=zeros(2,iter+1); %4
y_obs4=zeros(2,iter+1);

sigma=zeros(1,iter+1);

x_init=[0;0];
x(:,1)=x_init;
x_ref(:,1)=x_init;
x_obs2(:,1)=x_init;
x_obs3(:,1)=x_init;
x_obs4(:,1)=x_init;

u=zeros(3,iter+1);
eta=-10^(-5) + 2*10^(-5).*rand([2 iter+1]);
w=-10^(-3) + 2*10^(-3).*rand([1 iter+1]);
decision=4;
t=0;

F2_s=diag([1 0.5 1]);
F3_s=diag([1 1 0.3]);

u_ref=u_ref4;
F_real=F4;
F_l=F4;
K_p1=K41;
K_p2=K42;

decision_list=zeros(1,iter);

for k=1:iter
%% 故障发生时间点
    if k==230
        F_real=F2_s;
    end
    if k==446
        F_real=F4;
    end
    if k==714
        F_real=F3_s;
    end
    if k==980
        F_real=F2_s;
    end
    if k==1245
        F_real=F3_s;
    end
    if k==1480
        F_real=F4;
    end
    if k==1750
        F_real=F3_s;
    end
%% 模型迭代
    y(:,k) = C*x(:,k)+eta(:,k);
    
    if decision==4
        x_obs=x_obs4(:,k);
    end
    if decision==3
        x_obs=x_obs3(:,k);
    end
    if decision==2
        x_obs=x_obs2(:,k);
    end
    sigma(k+1)=sigma(k)+H*y(:,k)-Cv*x_ref(:,k); %控制部分
    mid=x_obs-x_ref(:,k);
    u(:,k+1)=F_l*(-K_p1*mid-K_p2*sigma(k)+u_ref);
    
    x(:,k+1)=A*x(:,k)+B*F_real*u(:,k+1)+E*w(k);
    x_ref(:,k+1) = A*x_ref(:,k)+B*F_l*u_ref;
    
    y_obs4(:,k) = C*x_obs4(:,k);
    x_obs4(:,k+1)=(A-L*C)*x_obs4(:,k)+B*F4*u(:,k+1)+L*y(:,k);

    y_obs2(:,k) = C*x_obs2(:,k);
    x_obs2(:,k+1)=(A-L*C)*x_obs2(:,k)+B*F2*u(:,k+1)+L*y(:,k);
    
    y_obs3(:,k) = C*x_obs3(:,k);
    x_obs3(:,k+1)=(A-L*C)*x_obs3(:,k)+B*F3*u(:,k+1)+L*y(:,k);
%% 误差选择
    if decision==2
        error=y(:,k)-y_obs2(:,k);
    end
    if decision==3
        error=y(:,k)-y_obs3(:,k);
    end
    if decision==4
        error=y(:,k)-y_obs4(:,k);
    end
%% FDI
    if t>100  
        if decision==4
            if Epsilon4_f2.contains(error)
                decision=2;
                t=0;
            end
            if Epsilon4_f3.contains(error)
                decision=3;
                t=0;
            end
        end
        
        if decision==3
            if Epsilon3_f2.contains(error)
                decision=2;
                t=0;
            end
        else
            if decision==2
                if Epsilon2_f3.contains(error)
                    decision=3;
                    t=0;
                end
            end
        end  
        if t==0 %% 模型调整
            if decision==3
                u_ref=u_ref3;
                F_real=F3;
                F_l=F3;
                K_p1=K31;
                K_p2=K32;
            else if decision==2
                    u_ref=u_ref2;
                    F_real=F2;
                    F_l=F2;
                    K_p1=K21;
                    K_p2=K22;
                end
            end
        end
    end
    t=t+1;
    decision_list(k)=decision;           
end

%% plot
Fault_list=4*ones(1,iter);
for i=230:445
    Fault_list(i)=2;
end
for i=446:713
    Fault_list(i)=4;
end
for i=714:979
    Fault_list(i)=3;
end
for i=980:1244
    Fault_list(i)=2;
end
for i=1245:1479
    Fault_list(i)=3;
end
for i=1480:1749
    Fault_list(i)=4;
end
for i=1750:2000
    Fault_list(i)=3;
end
figure(1)
subplot(311)
plot(decision_list,'black','LineWidth',2)
hold on
plot(Fault_list,'-.r')
xlim([0,2000]);
ylim([0,5]);
ylabel('Indices')
legend('FDI','Fault')

subplot(312)
plot(x(1,:)+0.4,'b','LineWidth',1)
xlim([0,2000]);
ylim([0.395,0.405]);
ylabel('Tank 1 Level')
legend('Plant')

subplot(313)
plot(x(2,:)+0.06,'b','LineWidth',1)
xlim([0,2000]);
ylim([0.06,0.1]);
xlabel('Time')
ylabel('Tank 2 Level')
legend('Plant')

ref=ones(3,iter+1);
for i=1:230
    ref(:,i)=u_ref4-d_ref;
end
for i=231:714
    ref(:,i)=u_ref2;
end
for i=715:980
    ref(:,i)=u_ref3;
end
for i=981:1245
    ref(:,i)=u_ref2;
end
for i=1246:2001
    ref(:,i)=u_ref3;
end

figure(2)
subplot(311)
plot(u(1,:),'b','LineWidth',1);
hold on
plot(ref(1,:),'-.g');
xlim([0,2000]);
ylim([-0.5,0.5]);
ylabel('U1')
legend('Plant input','Reference','Location','SouthEast')

subplot(312)
plot(u(2,:),'b','LineWidth',1);
hold on
plot(ref(2,:),'-.g');
xlim([0,2000]);
ylim([-0.5,0.5]);
ylabel('U2')
legend('Plant input','Reference','Location','SouthEast')

subplot(313)
plot(u(3,:),'b','LineWidth',1);
hold on
plot(ref(3,:),'-.g');
xlim([0,2000]);
ylim([-0.5,0.5]);
xlabel('Time')
ylabel('U3')
legend('Plant input','Reference','Location','SouthEast')



