function F = f_kan(y)
% ������������ �������
global c b r1 r2

F = y(2)^2 / y(1) + (y(2) + r2*y(1))*(r1 - b*y(1));