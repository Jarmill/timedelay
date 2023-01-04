%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 6;      %time horizon
% tau = 1;       % lag
tau = 1;       % lag
C0 = -1;
R0 = 0.2;

%JOINT+COMPONENT VALUES 
order=1; %[1.50000000092443]
% order=2; %[1.50000000185651]
% order=3; %[1.48438788315234]
% order=4;[1.45612716685445]
% order=5; %[1.45107294132040]
% order=6; %[1.44799863924451]
% order=7; %[1.44619915775023]
% order=8;[1.44471530756278]
% order=9; %[1.44267033567268]
order=10; %[1.43345049576323]
% order=11;
% order=12;

%% system variables
mpol('t', 1, 1);
mpol('x', 1, 1);
mpol('x_lag', 1, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

% f = -K0*x -K1*x_lag;
f =0.5*x(1) - x_lag(1);

lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
lsupp= lsupp.set_box([-1.5, 1.5]);
lsupp.X_init = sum((x-C0).^2) <= R0^2;
% lsupp.X_history = (x-x00).^2 <= Rh;
% lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];

% p = -x;
p = x(1);

PM = peak_delay_manager_base_split(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)

mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

[sol.obj_rec, opt_term]
