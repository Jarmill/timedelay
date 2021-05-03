%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
% Tmax = 5;      %time horizon
Tmax = 3;
% tau = 1;       % lag
tau = 0.75;       % lag
C0 = [1.5; 0];
R0 = 0.4;

% x0 = C0 - [R0; 0];
x0 = C0;

xh0 = x0;   %constant history x(t) = xh0 for times [-tau, 0]
x00 = xh0;  %no discontinuity at time 0

order = 3;
% d = order*2;

%% system variables
mpol('t', 1, 1);
mpol('x', 2, 1);
mpol('x_lag', 2, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

% f = -K0*x -K1*x_lag;
f = [x(2); -x_lag(1) + (1/3).*x(1).^3- x(2)];

lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
% lsupp.X = (x.^2 <= 4);
lsupp= lsupp.set_box([-1, 1.5; -1, 1]);
lsupp.X_init = x00;
lsupp.X_history = xh0;
% lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];

% p = -x;
p = -x(2);

PM = peak_delay_manager_base(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)
sol.obj_rec
mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

[-sol.obj_rec, opt_term]
