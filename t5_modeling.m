function [time, F, U] = t5_modeling(X_0, t_0, t_end)
% В данной функции производится интегрирование системы по времени
global control_arr
global count

count = length(control_arr(:,1));

control_arr = sortrows(control_arr,1);

options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
global Mas_u % массив управлений

% коэффициенты системы
global c b r1 r2

[T, X] = ode45(@t5_control_system, [t_0, t_end], [X_0(1) X_0(2)], options);

F = X;
time = T;
U = Mas_u;
