clc;close all;clear;
%% 参数设置
for k = 1:10

SOC0 = [15,13,0,5,0,10,20,2,0,10];
Dmax = [100,200,150,90,180,200,210,190,220,130];
Dmin = [80,140,120,80,150,160,190,185,192,98];
Smax = [150,150,150,85,175,185,196,190,210,150];
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
for j = 1:k
    SOC0(j) = 0;
    SOCmax(j) = 0;
end
%% 集中式——参与共享
D = sdpvar(1,10);
S = sdpvar(1,10);
Qch = sdpvar(1,10);
Qdis = sdpvar(1,10);
E = sdpvar(1,10);
SOC = sdpvar(1,10);
uch = binvar(1,10);
C1 = [];
C1 = [C1;sum(E)==0];
for i = 1:1:10
   C1 = [C1;Qch(i) + D(i)==Qdis(i) + S(i) + E(i);...
       SOC(i) == SOC0(i) + 0.95*Qch(i) - Qdis(i)/0.95;...
       SOC(i) <=SOCmax(i);SOC(i)>=0;...
       D(i) <= Dmax(i);D(i) >=Dmin(i);S(i) <= Smax(i);S(i) >=0;...
       E(i)>=-Emax(i);E(i)<=Emax(i);...
       Qch(i) <= uch(i)*Qchmax(i);Qch(i) >= 0;...
       Qdis(i) <= ((1-uch(i))*Qdismax(i));Qdis(i) >= 0];
end
z1 = -sum(a.*D - 0.5*b.*(D.^2)-alpha.*S - 0.5*beta.*(S.^2)-...
    gamma * (Qch + Qdis)+cch.*Qch - ...
    0.5*dch.*(Qch.^2)-cdis.*Qdis - 0.5*ddis.*(Qdis.^2));
ops=sdpsettings('solver','Gurobi','verbose',0,'debug',1);
reuslt=optimize(C1,z1,ops);
if reuslt.problem==0 
    ZZZ(k) = -value(z1);
    res1_E = value(E);
    res1_D = value(D);
    res1_S = value(S);
    res1_Qch = value(Qch);
    res1_Qdis = value(Qdis);
    res1_SOC = value(SOC);
else
    disp('求解出错')
end
     
end
ZZZ = [-287.05 ZZZ];
Z = fliplr(ZZZ);
figure(1)
x = 1:11;
plot(x,Z,'-*','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',14,'box','off');
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlabel('Number of Prosumers with Householde Energy Storage','FontSize',14,'FontName','Times New Roman')
ylabel('Social Welfare/$','FontSize',14,'FontName','Times New Roman')
text(8.8,Z(11),'-287.05','FontSize',14,'FontName','Times New Roman');