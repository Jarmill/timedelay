mset clear

lags = [0.5, 0.8];

n =2;
mpol('t', 1, 1);
mpol('x', n, 1);
mpol('xd', n, 2);

vars = struct('t', t, 'x', x, 'x_lag', xd);

X = x.^2 <= 2;

x0 = [0.8; 0];
X_init = (x == x0);
Tmax = 10;

%SLIR model
%susceptible, latent, infectious, removed

f = [x(2); -x(1) - x(2) + 0.5*xd(1, 1) - 0.3*xd(2, 2)];
 
 
% above is a dumb model, but it is good for debugging 
%  %latent: L' = beta S(t) I(t) - beta S(t-tau) I(t-tau)
% %xd(3, 1) infected 
% f = [-beta*x(1)*x(3);
%      beta*xd(1, 1)*xd(3, 1) - gamma*x(3)];
 
lsupp = delay_support(vars);
lsupp.lags = lags;
lsupp.vars = vars;
lsupp.X = X;
lsupp.X_init = X_init;

objective = x(2);

%% relaxation information
order = 2;
d = 2*order;

%% split occupation measure

loc = delay_location_prop_split(lsupp,f, objective);

% [cons, len_consistency] = loc.consistency_con(d);

[objective, cons_eq, cons_ineq, len_dual]  = loc.all_cons(d);