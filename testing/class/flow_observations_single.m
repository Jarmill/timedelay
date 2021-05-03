PLOT = 1;
SOLVE = 1;

T = 40;
tau = 0.75;

C0 = [1.5; 0];
R0 = 0.4;

% x0 = C0 - [R0; 0];
x0 = C0;


f_ode = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

f_dde = @(t, x, Z) [x(2); -Z(1) + (1/3).*x(1).^3- x(2)];

%start with constant history
f_history = @(t) x0;

ode_options =   odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.01);
[t_ode, x_ode] = ode45(f_ode, [0, T], x0, ode_options);

dde_options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);
sol_dde = dde23(f_dde, tau, f_history, [0, T], dde_options);
t_dde = sol_dde.x;
x_dde = sol_dde.y';


if PLOT
figure(50)
clf
hold on
cl = linspecer(4);
plot(x_ode(:, 1), x_ode(:, 2), 'LineWidth', 2, 'color', cl(1, :))
plot(x_dde(:, 1), x_dde(:, 2), 'LineWidth', 2, 'color', cl(4, :))
scatter(C0(1), C0(2), 200, 'ok')
legend({'\tau = 0', ['\tau=', num2str(tau, 3)]}, 'location', 'northeast',...
    'fontsize', 12)
pbaspect([diff(xlim), diff(ylim), 1])
% plot(t_dde, x_dde)
title('Single Point, Constant History', 'FontSize', 16)
end