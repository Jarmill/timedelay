beta = 5;
gamma1 = 2;
gamma2 = 3;


xh0 = [0.9; 0.1];
x0 = xh0;
Tmax = 10;

f1 = @(t, x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma1*x(2)];
f2 = @(t, x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma2*x(2)];

p1 = struct;
p1.beta = 0.5;
p1.gamma = 0.1;

p2 = struct;
p2.beta = 0.5;
p2.gamma = 0.15;

sir_delay_1= @(t,x,Z) Tmax * sir_delay(t,x,Z,p1);
sir_delay_2 = @(t,x,Z) Tmax * sir_delay(t,x,Z,p2);
sir_history = @(t) x0*(t ==0) + xh0*(t ~= 0);

opts = odeset('reltol', 1e-6);
tau_1 = 1;
tau_2 = 3;

sol_1 = dde23(sir_delay_1, tau_1, sir_history,  [0, T], opts);
sol_2 = dde23(sir_delay_2, tau_2, sir_history,  [0, T], opts);

%% plot
FS_title = 24;

figure(5)
clf

tiledlayout(1, 2)
ax1=nexttile;
hold on
plot(sol_1.x, sol_1.y(2, :), 'linewidth', 3)

[m1, i1] = max(sol_1.y(2, :));
scatter(sol_1.x(i1), sol_1.y(2, i1), 400, '*k', 'Linewidth', 3)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$I(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('Peak Estimation', 'fontsize', FS_title)

ylim([0, m1+0.05])


ax2=nexttile;

hold on
plot(sol_1.x, sol_1.y(2, :), 'linewidth', 3)
plot(sol_2.x, sol_2.y(2, :), 'linewidth', 3)
[m2, i2] = max(sol_2.y(2, :));
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$I(t)$', 'interpreter', 'latex', 'fontsize', 14)
title('Peak-Minimizing Control', 'fontsize', FS_title)
ylim([0, m1+0.05])

scatter(sol_2.x(i2), sol_2.y(2, i2), 200, 'ok', 'filled', 'Linewidth', 3)
plot(xlim, [m2, m2], '--k', 'linewidth', 2)

function dydt = sir(t, y, p)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*y(1)*y(2) - p.gamma*(y(2))];
end

function dydt = sir_delay(t, y, Z, p)
    %Z = y(t - lag)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*Z(1)*Z(2) - p.gamma*(y(2))];
end
% 
% function s = sir_history_ref(t, x0)
%     s = x0;
% end