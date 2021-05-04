%show a set of trajectories with constant histories

SAMPLE = 1;
PLOT = 1;

% K0 = 0;
% K1 = 1;
K0 = 2;
K1 = 2;


% hlim = [0.5, 1.5];
xh0 = 0.5;
x0 = 1;
tau = 1;
Tmax = 5;

exp_delay = @(t,y,Z) -K0*y - K1*Z;
exp_history = @(t) xh0;
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'MaxStep', 0.1);
if SAMPLE
%     out_sim = cell(Nh, 1);
        
        exp_history_curr = @(t) t*((xh0-x0)/tau) + x0;
        traj = dde23(exp_delay, tau, exp_history_curr, [0,Tmax], options);
%         out_sim{i} = traj;
end

% cl = linspecer(3);
cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cz = [0.915294117647059,0.281568627450980,0.287843137254902];
cp = [0.441568627450980,0.749019607843137,0.432156862745098];
FS = 16;
figure(1)
clf
hold on 


tn = [-tau; 0];
xn = [xh0; x0];

tz = linspace(0, Tmax-tau, 100*(Tmax-tau));
xz = deval(traj, tz);

tp = linspace(Tmax-tau, Tmax, 100*tau);
xp = deval(traj, tp);

tiledlayout(2, 1)
ax1 = nexttile;
hold on
plot(tz, xz, 'Color', cz, 'LineWidth', 3)
plot(tp, xp, 'Color', cp, 'LineWidth', 3)
ylabel('x(t)')
% xlabel('t')
% ylim([-2, 3])
title('Current Trajectory', 'FontSize', FS)


% plot([-tau, traj.x], [xh0, traj.y], 'color', cl(2, :), 'LineWidth', 2)

ax2 = nexttile;
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz + tau, xz, 'Color', cz, 'LineWidth', 3)
title('Delayed Trajectory', 'FontSize', FS)
ylabel('x(t-\tau)')

linkaxes([ax1, ax2])



xlabel('time')
% ylabel('x(t)')
% title('Same Initial Point', 'FontSize', 14)
