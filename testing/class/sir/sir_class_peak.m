%test the weak solution manager
%compare againt single_delay/discrete_weak

%% system properties
order = 4;

POINT_HISTORY = 0;

%% system variables
mpol('t', 1, 1);
mpol('x', 2, 1);
mpol('x_lag', 2, 1);

vars = struct('t', t, 'x', x, 'x_lag', x_lag);

beta = 0.4;
gamma = 0.1;

Tmax = 30;  %max time of simulation (days)
tau = 9;    %incubation period (days)

I0 = 0.2;   %initial infection rate, constant history
if POINT_HISTORY
x00 = [1-I0; I0];
xh0 = x00;
else
    x00 = [sum(x) >= 0; x>= 0; x(2) <= I0];
    xh0 = x00;
end
% xh0 = [1; 0];

% f = -K0*x -K1*x_lag;
f = [-beta*x(1)*x(2); beta*x_lag(1)*x_lag(2) - gamma*x(2)];

lsupp = delay_support(vars);
lsupp.lags = tau;
lsupp.Tmax = Tmax;
lsupp.vars = vars;
lsupp.X = [x >= 0; sum(x) <= 1];
lsupp.X_init = x00;
lsupp.X_history = xh0;
lsupp.SCALE_TIME = 1;

PM = peak_delay_manager_base(lsupp, f, x(2));

sol = PM.run(order);

[opt, mom_out, corner] = PM.recover();
peak_val = sol.obj_rec;

%order 3 in free history: 0.6032
%% empirical moments
if POINT_HISTORY
sir_delay = @(t,y,Z) Tmax * [-beta*y(1)*y(2);
            beta*Z(1)*Z(2) - gamma*(y(2))];
sir_history = @(t) xh0 * (t ~= 0) + x00 * (t == 0);

Trange = linspace(0, Tmax, 1000);

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);
traj = dde23(sir_delay, tau/Tmax, sir_history, Trange/Tmax, options);
end
% d = 2*order;
% dv3 = genPowGlopti(3,d);
% m_traj3 = monom_int(traj.x, traj.y, dv3);
% dv4 = genPowGlopti(3,d+2);
% m_traj4 = monom_int(traj.x, traj.y, dv4);
% 
% save('sir_weak.mat', 'm_traj3', 'm_traj4', 'mv_sys3', 'mv_sys4', 'traj');
