%show a set of trajectories with constant histories

SAMPLE = 1;
PLOT = 1;

% K0 = 0;
% K1 = 1;
K0 = 2;
K1 = 2;


hlim = [0.5, 1.5];
Nh = 7;

xh0_list = linspace(hlim(1), hlim(2), Nh);
xh0 = mean(hlim);
tau = 1;
Tmax = 5;

exp_delay = @(t,y,Z) -K0*y - K1*Z;
exp_history = @(t) xh0;
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'MaxStep', 0.1);
if SAMPLE
    out_sim = cell(Nh, 1);
    for i = 1:Nh
        xh0_curr = xh0_list(i);
        exp_history_curr = @(t) t*((xh0_curr-xh0)/tau) + xh0;
        traj = dde23(exp_delay, tau, exp_history_curr, [0,Tmax], options);
        out_sim{i} = traj;
    end
end

cl = linspecer(1);

figure(1)
clf
hold on 
for i = 1:Nh
    plot([-tau, out_sim{i}.x], [xh0_list(i), out_sim{i}.y], 'color', cl)
end
plot([0, 0], ylim, ':k', 'LineWidth', 2)
% patch([-tau, 0, 0, -tau, -tau], [hlim(1), hlim(1), hlim(2), hlim(2), hlim(1)], ...
%     'k', 'LineWidth', 3,  'FaceColor', 'None')

% patch([
% axis off

xlabel('time')
ylabel('x(t)')
title('Same Initial Point', 'FontSize', 14)
