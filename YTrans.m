function [y, dy] = YTrans(X)
% ������� � ������������ ����������

% ������������ �������
global c b r1 r2

y = X(2);
dy = X(2)*(-r2 + c*X(1));