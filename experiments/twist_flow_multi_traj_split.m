%test the peak estimation manager
%compare againt single_delay/discrete_peak

clear all
mset clear

%% system properties
Tmax = 8;
% T = 5;
% tau = 0.5;
tau = 0.5;

C0 = [-0.5; 0; 0];
R0 = 0.2;


%SPLIT-JOINT values (always better than joint+component at the cost of more
%higher computational complexity)
% order = 1; %1.200000000949830
% order = 2; %1.200000000949830
% order = 3; %1.200000000949830
order = 4;  %
% order = 5; % 


%% system variables
mpol('t', 1, 1);
mpol('x', 3, 1);
mpol('x_lag', 3, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f = A_true*x - B_true*(4*(x).*(x_lag.^2) - 3*x);
      
lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
% lsupp.X = (x.^2 <= 4);
lsupp= lsupp.set_box([-1.1, 0.75; -1.5, 2; -1, 1.2]);
lsupp.X_init = sum((x-C0).^2) <= R0^2;
% lsupp.X_history = (x-x00).^2 <= Rh;
% lsupp.X_history = (x^2 <= 0.25);
% lsupp.X_history = [x >= xh0; x <= 0];


p = x(3);

PM = peak_delay_manager_base_split(lsupp, f, p);

sol = PM.run(order);
% [objective, mom_con, supp_con, len_dual] = WM.cons(d, Tmax)

mv_sys = double(PM.loc.sys{1}.meas_occ.mom_monom_marg(0, 2*order));
[opt_term, mom_term, corner_term]=PM.loc.term.recover();

[sol.obj_rec]
save('twist_4.mat', 'sol')
