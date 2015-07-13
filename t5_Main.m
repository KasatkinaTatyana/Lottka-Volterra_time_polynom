clear all
close all
clc
%% Ограничения
%  constr1 < Psi(t) < constr2
%  dy_constr1 < dPsi(t) / dt < dy_constr2
dy_constr1 = -31;
dy_constr2 = 36;
constr2 = 27;

% начальные условия для исходной системы 
% система Лоттки-Вольтерра X = (x y)
% X_0 = [14 1];
% X_end = [4.5 10];

% коэффициенты системы
global c b r1 r2
c = 2;
b = 0.1;
r1 = 4;
r2 = 8;

% начальные условия для канонической системы
% [y0, dy0] = YTrans(X_0);
% [yend, dyend] = YTrans(X_end);
y0 = 1;
dy0 = 10;
yend = 10;
dyend = 20;

X_0 = [1/c*(dy0/y0 + r2) y0];
X_end = [1/c*(dyend/yend + r2) yend];


u0 = 0.5;
uend = 1;
ddy0 = f_kan([y0 dy0]) + g_kan([y0 dy0])*u0;
ddyend = f_kan([yend dyend]) + g_kan([yend dyend])*uend;

Y_0 =   [y0    dy0    ddy0  ];
Y_end = [yend  dyend  ddyend];

t_0 = 0;
t_end = 2;

% t5_trajectory_synthesis(Y_0, Y_end, 0, 100, 0.1);

c_t = t5_calc_coefs(Y_0,Y_end,t_0,t_end);
c0 = c_t(1);
c1 = c_t(2);
c2 = c_t(3);
c3 = c_t(4);
c4 = c_t(5);
c5 = c_t(6);
%% Массив, в котором будут содержаться коэффициенты для управлений
global control_arr
control_arr = [t_0 t_end c0 c1 c2 c3 c4 c5 0 0 t_0 t_end];
% 1 : t_0 - начальная точка кривой;
% 2 : t_end - конечная точка кривой;
% 3-8 : коэффициенты;
% 9 : d - параметр;
% 10 :  тип кривой
% 11 : t_0_c - это значение, которое вычитается из t 
% 12 : t_end_c  

% типы будут следующие:
% 0 - кривая вида Psi_0(t) = c0 + 1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - 
% t_0_c)^2*(t - t_0_c)^2 + c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3;
% 1 - кривая вида 
% 2 - кривая вида Psi_2(t) = Psi_0(t) + d*t_norm.^2.*(3 - 2*t_norm), 
% где y_norm = (t - t_0) / (t_end - t_0);

%%
figure(1);
hold on; grid on;
% title('Psi(t)');
xlabel('t, c');
ylabel('y');

N = 1000;
dt = t_end / N;
t = t_0:dt:t_end;

Psi = c0 + c1/(t_end - t_0)*(t - t_0) + c2/(t_end - t_0)^2*(t - t_0).^2 + c3/(t_end - t_0)^3*(t - t_0).^3 + ...
      c4/(t_end - t_0)^4*(t - t_0).^4 + c5/(t_end - t_0)^5*(t - t_0).^5;  % y(t)
dPsi = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t - t_0) + 3*c3/(t_end - t_0)^3*(t - t_0).^2 + ...
       4*c4/(t_end - t_0)^4*(t - t_0).^3 + 5*c5/(t_end - t_0)^5*(t - t_0).^4;               % \dot{y}  
ddPsi = 2*c2/(t_end - t_0)^2 + 6*c3/(t_end - t_0)^3*(t - t_0) + ...
        12*c4/(t_end - t_0)^4*(t - t_0).^2 + 20*c5/(t_end - t_0)^5*(t - t_0).^3; % \ddot{y}

F = zeros(1, N+1);
G = zeros(1, N+1);
for i=1:length(Psi)
    F(i) = f_kan([Psi(i) dPsi(i)]);
    G(i) = g_kan([Psi(i) dPsi(i)]);
end
u = (ddPsi - F)./ G;

plot(t,Psi);

figure(2);
hold on; grid on;
xlabel('t, c');
ylabel('dy / dt');
% title('dPsi');
plot(t,dPsi,'b');

figure(3);
hold on; grid on;
xlabel('t, c');
ylabel('u');
% title('Управление u(t)');
plot(t,u);

%% Убираю выход за ограничение dy / dt < dy_constr2
% [data1, data2, par, coefs] = t5_up_dPsi(Y_0, Y_end, t_0, t_end, dt, dy_constr2);
% t_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% ddPsi_1 = data1(4,:);
% N1 = length(Psi_1);
% F = []; G = [];
% F = zeros(1, N1);
% G = zeros(1, N1);
% for i=1:N1
%     F(i) = f_kan([Psi_1(i) dPsi_1(i)]);
%     G(i) = g_kan([Psi_1(i) dPsi_1(i)]);
% end
% u_1 = (ddPsi_1 - F)./G;
% 
% t_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% ddPsi_2 = data2(4,:);
% N2 = length(Psi_2);
% F = []; G = [];
% F = zeros(1, N2);
% G = zeros(1, N2);
% for i=1:N2
%     F(i) = f_kan([Psi_2(i) dPsi_2(i)]);
%     G(i) = g_kan([Psi_2(i) dPsi_2(i)]);
% end
% u_2 = (ddPsi_2 - F)./G;
% 
% set(0,'CurrentFigure',1);
% plot(t_1, Psi_1, 'g','LineWidth',1);
% plot(t_2, Psi_2, 'm','LineWidth',1);
% 
% set(0,'CurrentFigure',2);
% plot(t_1, dPsi_1, 'g','LineWidth',1);
% plot(t_2, dPsi_2, 'm','LineWidth',1);
% 
% set(0,'CurrentFigure',3);
% plot(t_1, u_1, 'LineStyle','-','Color','g','LineWidth', 1); % color g
% plot(t_2, u_2, 'LineStyle','-','Color','m','LineWidth', 1); % color m
% 
% row = [t_1(1) t_1(end) coefs(1,:) par(1) 2 par(2) par(3)];
% 
% t5_add_control_branch(row);
%  
% row = [t_2(1) t_2(end) coefs(2,:) 0 0 t_2(1) t_2(end)];
%  
% t5_add_control_branch(row);
%% еще раз
% [data1, data2, par, coefs] = t5_up_dPsi([Psi_1(end) dPsi_1(end) ddPsi_1(end)], Y_end, t_1(end), t_end, dt, dy_constr2);
% t_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% ddPsi_1 = data1(4,:);
% u_1 = ddPsi_1 - sin(Psi_1);
% 
% t_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% ddPsi_2 = data2(4,:);
% u_2 = ddPsi_2 - sin(Psi_2);
% 
% set(0,'CurrentFigure',1);
% plot(t_1, Psi_1, 'LineStyle','--','Color', 'g','LineWidth',1);
% plot(t_2, Psi_2, 'LineStyle','--','Color', 'm','LineWidth',1);
% 
% set(0,'CurrentFigure',2);
% plot(t_1, dPsi_1, 'LineStyle','--','Color', 'g','LineWidth',1);
% plot(t_2, dPsi_2, 'LineStyle','--','Color', 'm','LineWidth',1);
% 
% set(0,'CurrentFigure',3);
% plot(t_1, u_1, 'LineStyle','--','Color','g','LineWidth', 1); % color g
% plot(t_2, u_2, 'LineStyle','--','Color','m','LineWidth', 1); % color m
% 
% row = [t_1(1) t_1(end) coefs(1,:) par(1) 2 par(2) par(3)];
% 
% t5_add_control_branch(row);
%  
% row = [t_2(1) t_2(end) coefs(2,:) 0 0 t_2(1) t_2(end)];
%  
% t5_add_control_branch(row);
%% Убираю выход за ограничение dy / dt > dy_constr1
% [data1, data2, par, coefs] = t5_down_dPsi(Y_0, Y_end, t_0, t_end, dt, dy_constr1);
% % [data1, data2, par, coefs] = t5_down_dPsi([Psi_1(end) dPsi_1(end) ddPsi_1(end)],...
% %                                         Y_end, t_1(end), t_end, dt, dy_constr1);
% t_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% ddPsi_1 = data1(4,:);
% N1 = length(Psi_1);
% F = []; G = [];
% F = zeros(1, N1);
% G = zeros(1, N1);
% for i=1:N1
%     F(i) = f_kan([Psi_1(i) dPsi_1(i)]);
%     G(i) = g_kan([Psi_1(i) dPsi_1(i)]);
% end
% u_1 = (ddPsi_1 - F)./G;
% 
% t_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% ddPsi_2 = data2(4,:);
% N2 = length(Psi_2);
% F = []; G = [];
% F = zeros(1, N2);
% G = zeros(1, N2);
% for i=1:N2
%     F(i) = f_kan([Psi_2(i) dPsi_2(i)]);
%     G(i) = g_kan([Psi_2(i) dPsi_2(i)]);
% end
% u_2 = (ddPsi_2 - F)./G;
% 
% set(0,'CurrentFigure',1);
% plot(t_1, Psi_1, 'LineStyle','-','Color','g','LineWidth',2);
% plot(t_2, Psi_2, 'LineStyle','-','Color','m','LineWidth',2);
% 
% set(0,'CurrentFigure',2);
% plot(t_1, dPsi_1, 'LineStyle','-','Color','g', 'LineWidth',2);
% plot(t_2, dPsi_2, 'LineStyle','-','Color','m','LineWidth',2);
% 
% set(0,'CurrentFigure',3);
% plot(t_1, u_1, 'LineStyle','-','Color','g','LineWidth', 2); % color g
% plot(t_2, u_2, 'LineStyle','-','Color','m','LineWidth', 2); % color m
% 
% row = [t_1(1) t_1(end) coefs(1,:) par(1) 2 par(2) par(3)];
% 
% t5_add_control_branch(row);
%  
% row = [t_2(1) t_2(end) coefs(2,:) 0 0 t_2(1) t_2(end)];
%  
% t5_add_control_branch(row);
%% еще раз
% [data1, data2, par, coefs] = t5_down_dPsi([Psi_1(end) dPsi_1(end) ddPsi_1(end)],...
%                                           Y_end, t_1(end), t_end, dt, dy_constr1);
% t_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% ddPsi_1 = data1(4,:);
% u_1 = ddPsi_1 - sin(Psi_1);
% 
% t_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% ddPsi_2 = data2(4,:);
% u_2 = ddPsi_2 - sin(Psi_2);
% 
% set(0,'CurrentFigure',1);
% plot(t_1, Psi_1, 'LineStyle','--','Color','g','LineWidth',1);
% plot(t_2, Psi_2, 'LineStyle','--','Color','m','LineWidth',1);
% 
% set(0,'CurrentFigure',2);
% plot(t_1, dPsi_1, 'LineStyle','--','Color','g', 'LineWidth',1);
% plot(t_2, dPsi_2, 'LineStyle','--','Color','m','LineWidth',1);
% 
% set(0,'CurrentFigure',3);
% plot(t_1, u_1, 'LineStyle','--','Color','g','LineWidth', 1); % color g
% plot(t_2, u_2, 'LineStyle','--','Color','m','LineWidth', 1); % color m
% 
% row = [t_1(1) t_1(end) coefs(1,:) par(1) 2 par(2) par(3)];
% 
% t5_add_control_branch(row);
%  
% row = [t_2(1) t_2(end) coefs(2,:) 0 0 t_2(1) t_2(end)];
%  
% t5_add_control_branch(row);
%% Убираю выход за ограничение y > constr2
[data1, par, coefs] = t5_upper_Psi(Y_0, Y_end, t_0, t_end, dt, constr2);
t_1 = data1(1,:);
Psi_1 = data1(2,:);
dPsi_1 = data1(3,:);
ddPsi_1 = data1(4,:);
N1 = length(Psi_1);
F = []; G = [];
F = zeros(1, N1);
G = zeros(1, N1);
for i=1:N1
    F(i) = f_kan([Psi_1(i) dPsi_1(i)]);
    G(i) = g_kan([Psi_1(i) dPsi_1(i)]);
end
u_1 = (ddPsi_1 - F)./G;

set(0,'CurrentFigure',1);
plot(t_1, Psi_1, 'm','LineWidth',1);

set(0,'CurrentFigure',2);
plot(t_1, dPsi_1, 'm','LineWidth',1);

set(0,'CurrentFigure',3);
plot(t_1, u_1, 'LineStyle','-','Color','m','LineWidth', 1); % color c

row = [t_1(1) t_1(end) coefs(1,:) par(1) 1 par(2) par(3)];

t5_add_control_branch(row);
%% интегрирование системы по времени
[time, traj, U] = t5_modeling(X_0, t_0, t_end);
% [time, traj, U] = t_modeling([Psi_1(1) dPsi_1(1)],t_1(1),t_1(end));

% figure(4);
% hold on; grid on
% xlabel('y');
% ylabel('dy / dt');
% plot(traj(:,1), traj(:,2),'b');

figure(4);
subplot(2,1,1);
plot(time,traj(:,1));
xlabel('t, c'); ylabel('x_1');

subplot(2,1,2);
plot(time,traj(:,2));
xlabel('t, c'); ylabel('x_2');

figure(5);
hold on; grid on;
xlabel('t, c');
ylabel('u(t)');
title('Стабилизирующее правление u(t)');
plot(U(:,1), U(:,2));


%% -----------------------------------------------------------------------
% c = 1;
% b = 1;
% r1 = 0.5;
% r2 = 0.5;
% 
% t = 0;
% i = 1;
% X0 = [1 2];
% Z0 = zeros(1,2);
% [Z0(1) Z0(2)] = YTrans(X0);
% X = X0;
% Z = Z0;
% T = 10;
% dt = 0.01;
% u = [];
% while (t <= T)
%     u(i) = 0.1;
%     
%     Xz(i,1) = 1/c*(Z(i,2) / Z(i,1) + r2);
%     Xz(i,2) = Z(i,1);
%     
%     Z(i+1,1) = Z(i,1) + Z(i,2)*dt;
%     Z(i+1,2) = Z(i,2) + (f_kan(Z(i,:)) + ...
%     + g_kan(Z(i,:))*u(i))*dt;
% 
%     X(i+1,1) = X(i,1) + X(i,1)*(r1 - u(i) - b*X(i,2))*dt;
%     X(i+1,2) = X(i,2) + X(i,2)*(- r2 +c*X(i,1))*dt;
%     
%     t(i + 1) = t(i) + dt;
%     i = i + 1;
% end
% 
% Xz(i,1) = 1/c*(Z(i,2) / Z(i,1) + r2);
% Xz(i,2) = Z(i,1);
% 
% figure()
% hold on;
% plot(t, X(:,1), 'r');
% plot(t, Xz(:,1), 'b');
% 
% figure()
% hold on;
% plot(t, X(:,2), 'r');
% plot(t, Xz(:,2), 'b');