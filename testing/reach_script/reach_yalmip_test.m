%% problem properties
tau = 0.3;
K0 = 0.5;
K1 = 4;
T = 1;
boxlim = 1;
order = 2;
d = 2*order;

xh0_max = -0.6;
xh0_min = -1;

%% sdpvariables
t = sdpvar(1,1);
x = sdpvar(1, 1);
x0 = sdpvar(1,1);
x1 = sdpvar(1,1);

f = -K0*x0-K1*x1;

[v, cv] = polynomial([t; x], d);
[w, cw] = polynomial([x], d);
[phi0, cphi0] = polynomial([t; x], d);
[phi1, cphi1] = polynomial([t; x], d);

v0 = replace(v, x, 0);
vT = replace(v, x, T);

cons = [];
coeff = [cv; cw; cphi0; cphi1];


%% support sets
X_base = [boxlim^2-x^2];
XT = struct('ineq', Xbase, 'eq', []);
X0 = struct('ineq', [(xh0_max-x)*(x-xh0_min)], 'eq', []);

Omega0 = struct('ineq', [t*(T-tau - t); X_base], 'eq', []);
Omega1 = struct('ineq', [(t-T-tau)*(T-t); X_base], 'eq', []);
Omeganeg1 = struct('ineq', [(0-t)*(t+tau); X_base], 'eq', []);

Xall = struct('ineq', [t*(T-t); boxlim^2-[x0^2; x1^2]], 'eq', []);

%% cost

%% constraints
