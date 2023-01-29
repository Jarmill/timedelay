mset clear

%% parameters
T = 1;
xh0 = -1;
kappa = 0.6;
K0 = 1;
K1 = 5;



%% solve
n =1;
mpol('t', 1, 1);
mpol('x', n, 1);
mpol('xd', n, 1);

f = -K0*x(1) - K1*xd(1);
vars = struct('t', t, 'x', x, 'x_lag', xd);
lags = kappa;


lsupp = delay_support(vars);
lsupp.lags = lags;
lsupp.vars = vars;
lsupp = lsupp.set_box([-1, 0.5]);
lsupp.X_init = xh0;
lsupp.PROPORTIONAL = 1;
lsupp.Tmax = 1;

p = x(1);

%% relaxation information
order = 7;
d = 2*order;

PM = peak_delay_manager_base_split(lsupp, f, p);

sol = PM.run(order);
disp(sol.obj_rec);

%% visualize

options = ddeset('AbsTol', 1e-10, 'RelTol', 1e-7);



% kappa = 0.2;
% K0 = 1;
% K1 = 2;

traj = ddesd(@(t, y, z) -K0*y - K1*z,@(t, y) kappa*t,@(t) xh0,[0,T], options);

figure(1)
clf
hold on
plot(traj.x, traj.y)
plot([0, T], [1,1]*sol.obj_rec, '--r', 'LineWidth', 3);
plot([0, T], [0, 0], ':k')

xlabel('$t$', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(', num2str(kappa),'t)$'], 'interpreter', 'latex', 'fontsize', 16)
hold off
