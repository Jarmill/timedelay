PLOT = 1;
SOLVE = 1;

T = 10;
Tsmall = 5;
tau = 0.5;

C0 = [-0.5; 0; 0];
% R0 = 0.4;

% x0 = C0 - [R0; 0];
x0 = C0;

A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

sigma = 0.1;
f_ode = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);
f_dde = @(t,x, Z) A_true*x - B_true*(4*(x).*(Z.^2) - 3*x);

%start with constant history
f_history = @(t) x0;

ode_options =   odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.01);
[t_ode, x_ode] = ode45(f_ode, [0, T], x0, ode_options);

dde_options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);
sol_dde = dde23(f_dde, tau, f_history, [0, T], dde_options);
t_dde = sol_dde.x;
x_dde = sol_dde.y';

t_dde_small = t_dde(t_dde <= Tsmall);
x_dde_small = x_dde(t_dde <= Tsmall, :);

if PLOT
figure(50)
clf
hold on
cl = linspecer(4);
plot3(x_ode(:, 1), x_ode(:, 2), x_ode(:, 3), 'LineWidth', 2, 'color', cl(1, :))
plot3(x_dde(:, 1), x_dde(:, 2), x_dde(:, 3), 'LineWidth', 2, 'color', cl(4, :))
plot3(x_dde_small(:, 1), x_dde_small(:, 2), x_dde_small(:, 3), 'LineWidth', 2, 'color', cl(2, :))
scatter(C0(1), C0(2), 200, 'ok')
legend({'\tau = 0', ['\tau=', num2str(tau, 3)]}, 'location', 'northeast',...
    'fontsize', 12)
pbaspect([diff(xlim), diff(ylim), diff(zlim)])
% plot(t_dde, x_dde)
title('Single Point, Constant History', 'FontSize', 16)
end