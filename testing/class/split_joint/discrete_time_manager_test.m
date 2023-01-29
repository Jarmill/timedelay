mset clear

lags = 3;

n = 2;
mpol('t', 1, 1);
mpol('x', n, 1);
mpol('xd', n, 1);

vars = struct('t', t, 'x', x, 'x_lag', xd);

X = [x.^2 <= 2];

x0 = [0.6; 0];
R0 = 0.2;
X_init = (sum((x-x0).^2) <= R0^2);
X_history = X_init;
Tmax = 10;
objective = x(2);

%simple linear system
% f = [x(2); 0.5];
% f = [x(2); -0.4*x(1) - 0.2*x(2)];
f = [x(2); -0.4*x(1) - 0.2*x(2)+0.5*xd(1)];
% f = [0; 0];
 
lsupp = delay_support(vars);
lsupp.lags = lags;
lsupp.vars = vars;
lsupp.X = X;
lsupp.X_init = X_init;
lsupp.X_history = X_history;
lsupp.Tmax = Tmax;
lsupp.DISCRETE_TIME = 1;

%% relaxation information
order = 4;
d = 2*order;

%% split occupation measure
PM = peak_delay_manager_base_split(lsupp, f, objective);
sol = PM.run(order);
