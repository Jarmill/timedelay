%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
x00 = xh0;  %no discontinuity at time 0
% x00 = 0;  %discontinuity at time 0

tau = 0.25;
K0 = 3;
K1 = 5;

order = 6;
% d = order*2;

%% system variables
mpol('t', 1, 1);
mpol('x', 1, 1);
mpol('x_lag', 1, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

f = -K0*x -K1*x_lag;

lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
lsupp.X = (x^2 <= xh0^2);
lsupp.X_init = (x==x00);
lsupp.X_history = (x == xh0);
lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];

p = x;

WM = peak_delay_manager_base(lsupp, f, p);

sol = WM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)
sol.obj_rec
mv_sys = double(WM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));