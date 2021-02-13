% Discrete system with long discrete delays
%plotting and testing


Tf = 30;
h = 6;
k0 = 0.25;
kh = 0.75;

coeff = zeros(1, h);
coeff(1) = k0;
coeff(end) = kh;

A = [-coeff;
    eye(h-1), zeros(h-1, 1)];
C = [1 zeros(1, h-1)];
sys = ss(A, [], C, [], 1);
x0 = ones(1, h);

% [y, t, x] = initial(sys, x0, Tf);
figure(9)
clf
initial(sys, x0, Tf);
abs(eig(A))