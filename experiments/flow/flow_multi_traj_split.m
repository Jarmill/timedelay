%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 5;      %time horizon
% tau = 1;       % lag
tau = 0.75;       % lag
C0 = [1.5; 0];
R0 = 0.4;

%JOINT+COMPONENT VALUES 
% order = 1; %-1.2500 , time: 0.7818 sec
% order = 2; %-1.2230, time: 0.9914 sec
% order = 3; %-1.1937 , time: 5.2708 sec
% order = 4;  %-1.1751, time:  31.8847 sec
% order = 5;  %-1.1636, time:  336.5093  sec
%SPLIT-JOINT values (always better than joint+component at the cost of more
%higher computational complexity)
% order = 1; %-1.2500 , time: 0.9372 sec
% order = 2; %-1.2183, time:  1.1900 sec
% order = 3; %-1.1913, time: 9.5082 sec
% order = 4;  %-1.1727, time: 105.7765 sec
order = 5; %-1.1630 
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
lsupp= lsupp.set_box([-1.25, 2.5; -1.25, 1.5]);
lsupp.X_init = sum((x-C0).^2) <= R0^2;
% lsupp.X_history = (x-x00).^2 <= Rh;
% lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];

% p = -x;
p = -x(2);

PM = peak_delay_manager_base_split(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)

mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

[-sol.obj_rec, opt_term, sol.solver_time]
