mset clear
clear

lags = [1; 2];

n = 3;
mpol('t', 1, 1);
mpol('x', n, 1);
mpol('xd', n, 2);

vars = struct('t', t, 'x', x, 'x_lag', xd);

X = [sum(x) <= 1; x >= 0];

x0 = [0.8; 0; 0.2];
% X_init = (x == x0);
X_init = [x(3) <= 0.1; sum(x) <= 1; x >= 0];
X_history = X_init;
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

%% relaxation information
order = 1;
d = 2*order;

objective = x(3);

loc = delay_location_base(lsupp, f, x(3));
% sm = meas_time_slack(lsupp);
% s1 = sm.mom_index(1, d)
% % s2 = sm.mom_index(2, d)

% lag_span = cm.lag_intervals()

[objective, cons_eq, cons_ineq, len_dual] = loc.all_cons(d);
