clc;close all;clear;
%% 参数设置
SOC0 = [15,13];
Dmax = [100,200];
Dmin = [80,140];
Smax = [150,150];
SOCmax = [100,90];
Qchmax = [6,5];
Qdismax = [6,5];
Emax = [80,80];
a = [2,0.008];
b = [0.001,0.009];
alpha = [0.03,0.006];
beta = [0.02,0.001];
cch = [5,0.0008];
cdis = [0.02,0.04];
dch = [0.0005,0.06];
ddis = [0.003,0.0001];
gamma = 0.0003;
%% 集中式——参与共享
D = sdpvar(1,2);
S = sdpvar(1,2);
Qch = sdpvar(1,2);
Qdis = sdpvar(1,2);
E = sdpvar(1,2);
SOC = sdpvar(1,2);
uch = binvar(1,2);
C1 = [];
C1 = [C1;sum(E)==0];
for i = 1:1:2
   C1 = [C1;Qch(i) + D(i)==Qdis(i) + S(i) + E(i);...
       SOC(i) == SOC0(i) + 0.95*Qch(i) - Qdis(i)/0.95;...
       SOC(i) <=SOCmax(1);SOC(i)>=0;...
       D(i) <= Dmax(i);D(i) >=Dmin(i);S(i) <= Smax(i);S(i) >=0;...
       Qch(i) <= uch(i)*Qchmax(i);Qch(i) >= 0;...
       Qdis(i) <= ((1-uch(i))*Qdismax(i));Qdis(i) >= 0];
end
z1 = -sum(a.*D - 0.5*b.*(D.^2)-alpha.*S - 0.5*beta.*(S.^2)-...
    gamma * (Qch + Qdis)+cch.*Qch - ...
    0.5*dch.*(Qch.^2)-cdis.*Qdis - 0.5*ddis.*(Qdis.^2));
ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
reuslt=optimize(C1,z1,ops);
if reuslt.problem==0 
    -value(z1)
    res1_E = value(E);
    res1_D = value(D);
    res1_S = value(S);
    res1_Qch = value(Qch);
    res1_Qdis = value(Qdis);
    res1_SOC = value(SOC);
else
    disp('求解出错')
end
%% 集中式——不参与共享
C2 = [];
for i = 1:1:2
   C2 = [C2;Qch(i) + D(i)==Qdis(i) + S(i);...
       SOC(i) == SOC0(i) + 0.95*Qch(i) - Qdis(i)/0.95;...
       SOC(i) <=SOCmax(1);SOC(i)>=0;...
       D(i) <= Dmax(i);D(i) >=Dmin(i);S(i) <= Smax(i);S(i) >=0;...
       Qch(i) <= uch(i)*Qchmax(i);Qch(i) >= 0;...
       Qdis(i) <= ((1-uch(i))*Qdismax(i));Qdis(i) >= 0];
end
z2 = -sum(a.*D - 0.5*b.*(D.^2)-alpha.*S - 0.5*beta.*(S.^2)-...
    gamma * (Qch + Qdis)+cch.*Qch - ...
    0.5*dch.*(Qch.^2)-cdis.*Qdis - 0.5*ddis.*(Qdis.^2));
ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
reuslt=optimize(C2,z2,ops);
if reuslt.problem==0 
    -value(z2)
    res2_D = value(D);
    res2_S = value(S);
    res2_Qch = value(Qch);
    res2_Qdis = value(Qdis);
    res2_SOC = value(SOC);
else
    disp('求解出错')
end
%% ADMM——参与共享
k = 1;%迭代次数
r(k) = 0;%初始共享电价
p = 0.1;%步长
epri = 0.0001;
edual = 0.0001;
aD(:,1) = [90,150];
aS(:,1) = [100,100];
aQch(:,1) = [0,0];
aQdis(:,1) = [0,0];
aE(:,1) = [0,0];
aSOC(:,1) = [15,13];
auch(:,1) = [0,0];
for k = 2:10000
    % 子问题1计算
    D = sdpvar(1,1);
    S = sdpvar(1,1);
    Qch = sdpvar(1,1);
    Qdis = sdpvar(1,1);
    E = sdpvar(1,1);
    SOC = sdpvar(1,1);
    uch = binvar(1,1);
    C3 = [];
    C3 = [C3;Qch + D==Qdis + S + E;...
            SOC == SOC0(1) + 0.95*Qch - Qdis/0.95;...
            SOC <=SOCmax(1);SOC>=0;...
            D <= Dmax(1);D >=Dmin(1);S <= Smax(1);S >=0;...
            Qch <= uch*Qchmax(1);Qch >= 0;...
            Qdis <= ((1-uch)*Qdismax(1));Qdis >= 0];
    z3 = -(a(1)*D - 0.5*b(1)*(D^2)-alpha(1)*S - 0.5*beta(1)*(S^2)-...
        gamma * (Qch + Qdis)+cch(1)*Qch - ...
        0.5*dch(1)*(Qch^2)-cdis(1)*Qdis - 0.5*ddis(1)*(Qdis^2)-r(k-1)*E)...
        +p/2*(aE(2,k-1)+E+r(k-1)/p)^2;
    ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
    reuslt=optimize(C3,z3,ops);
    if reuslt.problem==0
        az(1,k)=-value(z3);
        aE(1,k) = value(E);
        aD(1,k) = value(D);
        aS(1,k) = value(S);
        aQch(1,k) = value(Qch);
        aQdis(1,k) = value(Qdis);
        aSOC(1,k) = value(SOC);
    else
        disp('求解出错')
    end
    % 子问题2计算
    C4 = [];
    C4 = [C4;Qch + D==Qdis + S + E;...
            SOC == SOC0(2) + 0.95*Qch - Qdis/0.95;...
            SOC <=SOCmax(2);SOC>=0;...
            D <= Dmax(2);D >=Dmin(2);S <= Smax(2);S >=0;...
            Qch <= uch*Qchmax(2);Qch >= 0;...
            Qdis <= ((1-uch)*Qdismax(2));Qdis >= 0];
    z4 = -(a(2)*D - 0.5*b(2)*(D^2)-alpha(2)*S - 0.5*beta(2)*(S^2)-...
        gamma * (Qch + Qdis)+cch(2)*Qch - ...
        0.5*dch(2)*(Qch^2)-cdis(2)*Qdis - 0.5*ddis(2)*(Qdis^2)-r(k-1)*E)...
        +p/2*(aE(1,k)+E+r(k-1)/p)^2;
    ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
    reuslt=optimize(C4,z4,ops);
    if reuslt.problem==0
        az(2,k)=-value(z4);
        aE(2,k) = value(E);
        aD(2,k) = value(D);
        aS(2,k) = value(S);
        aQch(2,k) = value(Qch);
        aQdis(2,k) = value(Qdis);
        aSOC(2,k) = value(SOC);
    else
        disp('求解出错')
    end
    % 更新步长
    r(k) = r(k-1) + p*(aE(1,k) + aE(2,k));
    fprintf('第%d次迭代\n',k);
    fprintf('epri = %f\n',sqrt((aE(1,k) + aE(2,k))^2));
    if (sqrt((aE(1,k) + aE(2,k))^2)<=epri&&sqrt((aE(2,k)-aE(2,k-1))^2)<=edual)
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
figure(1)
edualiteration(1) = 0;
for i = 2:33
    edualiteration(i) = sqrt((aE(2,i)-aE(2,i-1))^2);
end
x = 1:33;
plot(x,sqrt((aE(1,:) + aE(2,:)).^2),'-*',x,edualiteration,'-*','LineWidth',1)
xlim([1 33])
set(gca,'FontName','Times New Roman','FontSize',14,'box','off');
xlabel('Iterations','FontSize',14,'FontName','Times New Roman')
ylabel('Residual','FontSize',14,'FontName','Times New Roman')
legend('Primal Residual','Dual Residual')
g1 = sum(res1_S)/sum(Smax);
g2 = sum(res2_S)/sum(Smax);



