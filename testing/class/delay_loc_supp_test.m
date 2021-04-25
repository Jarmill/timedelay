mset clear

lags = [1; 2];

n = 3;
mpol('t', 1, 1);
mpol('x', n, 1);
mpol('xd', n, 2);

vars = struct('t', t, 'x', x, 'x_lag', xd);

X = [sum(x) <= 1; x >= 0];

x0 = [0.8; 0; 0.2];
X_init = (x == x0);
X_history = (x == x0);
Tmax = 20;

%SLIR model
%susceptible, latent, infectious, removed

beta = 5;
gamma = 1;
alpha = 20;

%xd(3, 1) infected 
f = [-beta*x(1)*x(3);
     beta*xd(1, 1)*xd(3, 1) - alpha*xd(2, 2);
     alpha*xd(2, 2) - gamma*x(3)];
 
lsupp = delay_support(vars);
lsupp.lags = lags;
lsupp.vars = vars;
lsupp.X = X;
lsupp.X_init = X_init;
lsupp.X_history = X_history;

mpol('tq', 1, 1);
mpol('xq', n, 1);
mpol('xdq', n, 2);

varsq = struct;
varsq.t = tq;
varsq.x = xq;
varsq.x_lag = xdq;

qsupp = delay_support(varsq, lsupp);