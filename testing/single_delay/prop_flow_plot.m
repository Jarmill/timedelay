options = ddeset('AbsTol', 1e-10, 'RelTol', 1e-7);

rng(4, 'twister');
T = 10;
Ntraj = 150;

R0 = 0.3;
C0 = [1; 0];

xball = ball_sample(Ntraj, 2);
xh0 = C0 + R0*xball';
sol = cell(Ntraj, 1);

Ntheta = 200;
theta = linspace(0, 2*pi, Ntheta);

%% sample
kappa = 0.6;
% kappa = 0.5;
K0 = 1;
K1 = 2;

f_true = @(t, x, z) [x(2); -x(1) + (1/3).* x(1).^2 .* z(1) - z(2)];
% f_true = @(t, x, z) [x(2); -x(1) + (1/3).* x(1).^3 - z(2)];


for i = 1:Ntraj
    
    sol{i} = ddesd(@(t, y, z) f_true(t, y, z), @(t, y) kappa*t,@(t) xh0(:, i),[0,T], options);
end

%% plot
figure(1)
clf
hold on
for i = 1:Ntraj
    plot(sol{i}.y(1, :), sol{i}.y(2, :), 'c')
    scatter(sol{i}.y(1, 1), sol{i}.y(2, 1), 100, 'k')
end
% plot([0, T], [0, 0], ':k')
plot(cos(theta)*R0 + C0(1), sin(theta)*R0 + C0(2), 'k', 'LineWidth', 2)
xlabel('x_1(t)')
ylabel('x_1(t)')
title('Proportional Delayed Flow System', 'interpreter', 'latex', 'fontsize', 16)
hold off
pbaspect([diff(xlim), diff(ylim), 1])