clear all;clc
g = 1;
lp = 1;
mp =0.2;
mk = 1;
mt = 1;
a = g/(lp*(4.0/3 - mp/(mp+mk)));
A = [0 1 0 0;0 0 a 0;0 0 0 1;0 0 a 0];
b = -1/(lp*(4.0/3 - mp/(mp+mk)));
B = [0;1/mt;0;b];
C = [0 1 0 1];
D = 0;
sys = ss(A,B,C,D,0.01);
N = 100;
u = 0;
%%  Apply Disturbance 
w = wgn(N,1,1);
v = w;
%%  Open Loop System
x = .1*ones(4,1);
for i = 1:N
    x(:,i+1) = A*x(:,i)+B*u+w(i);
    y(i) = C*x(:,i)+v(i);
    x(:,i) = x(:,i+1);
end
figure(1);clf
plot(x');
grid on;
xlim([.1 N])
xlabel('time step')
ylabel('amplitude')
legend('x_1','x_2','x_3','x_4')
title('Open Loop System with Disturbance at Process and Measurement Signal')
% %%  LQE manual 
% [S,kk,ll] = idare(A',C',1e-5,1)
%%  LQG 
R = 1;
for i = 1:5
    Q(i) = 10^-i;
    [~,L(:,i),~] = kalman(sys,Q(i),R);
    [K(i,:),~,cl(:,i)] = lqr(sys,Q(i),R);
end
Rww = 1e-5;
Rvv = 1;
[LL,S] = KalmanConventional(sys,N,Rww,Rvv);
%%  Observer-based-controller 
x_hat = ones(4,1,5);
for i = 1:N
    for j = 1:5
    u_c(j,i) = -K(j,:)*x_hat(:,i,j);
    y_hat(j,i) = C*x_hat(:,i,j);
    x_hat(:,i+1,j) = (A-B*K(j,:))*x_hat(:,i,j)+L(:,j)*(y_hat(i)-y(i));
    end
end
%%  Plot 
figure(2);clf
lgn1 = sprintf('Q=1.00E-1');
lgn2 = sprintf('Q=1.00E-2');
lgn3 = sprintf('Q=1.00E-3');
lgn4 = sprintf('Q=1.00E-4');
lgn5 = sprintf('Q=1.00E-5');
figure(2);
for i = 1:4
    for j = 1:5
        subplot(2,2,i);
        plot(x_hat(i,:,j)');hold on
        xlabel('time step')
        grid on
        ylabel(sprintf('x_%d',i))
    end
end
legend(lgn1,lgn2,lgn3,lgn4,lgn5,'Location','SouthWest')
for i = 1:5
    figure(3)
    plot(y_hat(i,:));hold on
end
xlabel('time step')
ylabel('degree')
grid on
legend(lgn1,lgn2,lgn3,lgn4,lgn5,'Location','NorthWest')
title('Output Response with Various Matrices Q')
%%  Observer-based-controller via KalmanNet and RL  
x_nw = ones(4,1);
for i = 1:N
    x_nw(:,i+1) = A*x_nw(:,i);
    y_nw(i) = C*x_nw(:,i);
    x_nw(:,i) = x_nw(:,i+1);
end
x_nw = x_nw(:,1:N);
[x_hat_net,y_hat_net,KG] = KalmanNet(x_nw,y_nw,A,C);
%%
[P_vi,K_vi,G_vi] = value_iteration(A,N,B,Q(5),R,0.2);
x_hat_vi = ones(4,1);
for i = 1:N
    u_vi(j,i) = -K_vi*x_hat_vi(:,i);
%     y_hat_vi(i) = C*x_hat_vi(:,i);
    x_hat_vi(:,i+1) = (A-B*K_vi)*x_hat_vi(:,i)+KG(i)*(normalize(y_hat_net(i))-y(i));
end
%%  Plot 
figure(4);clf
subplot(211)
plot(x');hold on
subplot(212)
plot(x_hat_net')
figure(5);clf
plot(y_hat(5,:),'r')
hold on
plot(normalize(y_hat_net),'b')
xlim([5 N])
grid on
xlabel('time step')
ylabel('degree')
legend('LQG','KalmanNet-VI')
title('Comparison Output Response')