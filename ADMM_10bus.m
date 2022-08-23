clc;close all;clear;
%% 参数设置
SOC0 = [15,13,0,5,0,10,20,2,0,10];
Dmax = [100,200,150,90,180,200,210,190,220,130];
Dmin = [80,140,120,80,150,160,190,185,192,98];
Smax = [150,150,150,85,175,185,196,180,210,150];
SOCmax = [100,90,80,40,55,35,85,40,45,100];
Qchmax = [6,5,4,5,6,4,5,3,5,6];
Qdismax = [6,5,4,5,6,4,5,3,5,6];
Emax = [80,80,40,30,40,60,45,60,39,50];
a = [2,0.08,0.07,0.5,0.02,1.5,0.06,1.2,0.05,0.8];
b = [0.001,0.009,0.005,0.001,0.002,0.003,0.007,0.006,0.009,0.002];
alpha = [0.03,0.006,0.05,0.06,0.002,0.001,0.007,0.04,0.006,0.05];
beta = [0.02,0.001,0.01,0.02,0.005,0.01,0.002,0.003,0.001,0.02];
cch = [5,0.0008,0.003,0.2,0.02,0.007,0.02,5,0.7,0.002];
cdis = [0.02,0.04,0.01,0.001,0.02,0.08,0.1,0.002,0.001,0.004];
dch = [0.0005,0.06,0.003,0.01,0.08,0.0006,0.001,0.002,0.0003,0.01];
ddis = [0.003,0.0001,0.002,0.0008,0.004,0.0006,0.0007,0.0003,0.005,0.002];
gamma = 0.0003;
%% ADMM——参与共享——大配电网
k = 1;%迭代次数
r(k) = 0.1;%初始共享电价
p = 0.01;%步长
epri = 0.0001;
edual = 0.0001;
aE = [0,0,0,5,-5,0,0,0,0,0];
for k = 2:10000
    % 子问题1计算
    D = sdpvar(1,1);
    S = sdpvar(1,1);
    Qch = sdpvar(1,1);
    Qdis = sdpvar(1,1);
    E = sdpvar(1,1);
    SOC = sdpvar(1,1);
    uch = binvar(1,1);
    C = [];
    C = [C;Qch + D==Qdis + S + E;...
        SOC == SOC0(1) + 0.95*Qch - Qdis/0.95;...
        SOC <=SOCmax(1);SOC>=0;...
        D <= Dmax(1);D >=Dmin(1);S <= Smax(1);S >=0;...
        Qch <= uch*Qchmax(1);Qch >= 0;...
        E>=-Emax(1);E<=Emax(1);...
        Qdis <= ((1-uch)*Qdismax(1));Qdis >= 0];
    z = -(a(1)*D - 0.5*b(1)*(D^2)-alpha(1)*S - 0.5*beta(1)*(S^2)-...
        gamma * (Qch + Qdis)+cch(1)*Qch - ...
        0.5*dch(1)*(Qch^2)-cdis(1)*Qdis - 0.5*ddis(1)*(Qdis^2)-r(k-1)*E)...
        +p/2*(sum(aE)-aE(1)+E+r(k-1)/p)^2;
    ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
    reuslt=optimize(C,z,ops);
    if reuslt.problem==0
        az(1)=-value(z);
        aE(1) = value(E);
        aD(1) = value(D);
        aS(1) = value(S);
        aQch(1) = value(Qch);
        aQdis(1) = value(Qdis);
        aSOC(1) = value(SOC);
    else
        disp('求解出错')
    end

    % 子问题2计算
    for j = 2:10
        C = [];
        C = [C;Qch + D==Qdis + S + E;...
            SOC == SOC0(j) + 0.95*Qch - Qdis/0.95;...
            SOC <=SOCmax(j);SOC>=0;...
            D <= Dmax(j);D >=Dmin(j);S <= Smax(j);S >=0;...
            Qch <= uch*Qchmax(j);Qch >= 0;...
            E>=-Emax(j);E<=Emax(j);...
            Qdis <= ((1-uch)*Qdismax(j));Qdis >= 0];
        z = -(a(j)*D - 0.5*b(j)*(D^2)-alpha(j)*S - 0.5*beta(j)*(S^2)-...
            gamma * (Qch + Qdis)+cch(j)*Qch - ...
            0.5*dch(j)*(Qch^2)-cdis(j)*Qdis - 0.5*ddis(j)*(Qdis^2)-r(k-1)*E)...
            +p/2*(sum(aE)-aE(j)+E+r(k-1)/p)^2;
        ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
        reuslt=optimize(C,z,ops);
        if reuslt.problem==0
            az(j)=-value(z);
            aE(j) = value(E);
            aD(j) = value(D);
            aS(j) = value(S);
            aQch(j) = value(Qch);
            aQdis(j) = value(Qdis);
            aSOC(j) = value(SOC);
        else
            disp('求解出错')
        end
    end
    % 更新步长
    r(k) = r(k-1) + p*(sum(aE));
    fprintf('第%d次迭代\n',k);
    fprintf('epri = %f\n',sqrt((sum(aE))^2));
    if (sqrt((sum(aE))^2)<=epri)
        disp('迭代结束')
        break;
    end
end
Z = (a(1)*aD(1,k) - 0.5*b(1)*(aD(1,k)^2)-alpha(1)*aS(1,k) ...
    - 0.5*beta(1)*(aS(1,k)^2)-...
        gamma * (aQch(1,k) + aQdis(1,k))+cch(1)*aQch(1,k) - ...
        0.5*dch(1)*(aQch(1,k)^2)-cdis(1)*aQdis(1,k) ...
        - 0.5*ddis(1)*(aQdis(1,k)^2)-r(k)*aE(1,k))+...
        (a(2)*aD(2,k) - 0.5*b(2)*(aD(2,k)^2)-alpha(2)*aS(2,k) ...
    - 0.5*beta(2)*(aS(2,k)^2)-...
        gamma * (aQch(2,k) + aQdis(2,k))+cch(2)*aQch(2,k) - ...
        0.5*dch(2)*(aQch(2,k)^2)-cdis(2)*aQdis(2,k) ...
        - 0.5*ddis(2)*(aQdis(2,k)^2)-r(k)*aE(2,k));