function F = g_kan(y)
% коэффициенты системы
global c b r1 r2

F = -(y(2) + r2*y(1));