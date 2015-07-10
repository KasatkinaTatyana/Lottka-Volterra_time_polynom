function dx = t5_control_system(t,x)
dx = zeros(2,1);

global control_arr
global count

global Mas_u

% коэффициенты системы
global c b r1 r2

y = zeros(1,2);
[y(1) y(2)] = YTrans(x);

% y = x;

u = control(t,y);

% dx(1) = y(2);
% dx(2) = f_kan(y) + g_kan(y)*u;

dx(1) = x(1)*(r1 - u - b*x(2)); 
dx(2) = x(2)*(-r2 + c*x(1));