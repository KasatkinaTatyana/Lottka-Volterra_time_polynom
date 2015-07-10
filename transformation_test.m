clear all
close all
clc


X0 = [1 2];
% x0 y0

c = 1;
b = 1;
r1 = 0.5;
r2 = 0.5;

Z0 = [X0(2) X0(2)*(-r2 + c*X0(1))];

T = 10;
dt = 0.01;

t = 0;
i = 1;
X = X0;
Z = Z0;
while (t <= T)
    u(i) = 0.1;
    
    Xz(i,1) = 1/c*(Z(i,2) / Z(i,1) + r2);
    Xz(i,2) = Z(i,1);
    
    Z(i+1,1) = Z(i,1) + Z(i,2)*dt;
    Z(i+1,2) = Z(i,2) + (Z(i,2)^2 / Z(i,1) + (Z(i,2) + r2*Z(i,1))*(r1 - b*Z(i,1)) - ...
    (Z(i,2) + r2*Z(i,1))*u(i))*dt;

    X(i+1,1) = X(i,1) + X(i,1)*(r1 - u(i) - b*X(i,2))*dt;
    X(i+1,2) = X(i,2) + X(i,2)*(- r2 +c*X(i,1))*dt;
    
    t(i + 1) = t(i) + dt;
    i = i + 1;
end

Xz(i,1) = 1/c*(Z(i,2) / Z(i,1) + r2);
Xz(i,2) = Z(i,1);

figure()
hold on;
plot(t, X(:,1), 'r');
plot(t, Xz(:,1), 'b');

figure()
hold on;
plot(t, X(:,2), 'r');
plot(t, Xz(:,2), 'b');