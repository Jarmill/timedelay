%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 5;      %time horizon
% tau = 1;       % lag
tau = 0.5;       % lag
C0 = [1.5; 0];
R0 = 0.4;

order = 4; %success at order 4: -0.0351
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


%unsafe set 

%circle
% Cu = -0.8*[1; 1];
Cu = [-0.5; -1];
% Cu = [; -0.5];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

%line
%tilt angle
%order 5: 

theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1], safe

w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

%objective to maximize
p = [c1f; c2f];


PM = peak_delay_manager_base(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)

mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

[sol.obj_rec, opt_term]
