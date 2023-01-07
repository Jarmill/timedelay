%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 5;      %time horizon
% tau = 1;       % lag
tau = 0.75;       % lag
C0 = [-0.75; 0];
R0 = 0.3;




%% system variables
mpol('t', 1, 1);
mpol('x', 2, 1);
mpol('x_lag', 2, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

% f = -K0*x -K1*x_lag;
f = [x(2)*t - 0.1*x(1) - x_lag(1)*x_lag(2);
    -x(1)*t - x(2) + x(1)*x_lag(1)];

lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
% lsupp.X = (x.^2 <= 4);
lsupp= lsupp.set_box([-1.25, 1.25; -0.75, 1.25]);
lsupp.X_init = sum((x-C0).^2) <= R0^2;
% lsupp.X_history = (x-x00).^2 <= Rh;
% lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];

% p = -x;
% p = -x(2);

p = x(1);
% order=1; %[1.25000000472619]
% order=2; %[1.25000000472619]
% order=3; %[1.19781822680836]
% order=4; %[0.854326158386573]
% order=5; %[0.718264618058180]

% order = 5; %

PM = peak_delay_manager_base_split(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)

mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

sol.obj_rec
