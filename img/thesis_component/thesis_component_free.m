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
% cn = [0,0,1];
% cz = [1,0,0];
% cp = [0,1,0];
cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cz = [0.915294117647059,0.281568627450980,0.287843137254902];
cp = [0.441568627450980,0.749019607843137,0.432156862745098];
FS = 16;


%% Set up the trajectory 
tn = [-tau; 0];
xn = [xh0; x0];

tz = linspace(0, Tmax-tau, 100*(Tmax-tau));
xz = deval(traj, tz);

tp = linspace(Tmax-tau, Tmax, 100*tau);
xp = deval(traj, tp);


Tstop = 3;
xstop = deval(traj, Tstop);
xstop_delay = deval(traj, Tstop-tau);

%% Stacked Components
figure(1)
clf
hold on 

tiledlayout(2, 1)
ax1 = nexttile;
hold on
plot(tz(tz <= Tstop), xz(tz <= Tstop), 'Color', cz, 'LineWidth', 3)
scatter(Tstop, xstop, 200, '*k', 'LineWidth', 2)
% plot(tp, xp, 'Color', cp, 'LineWidth', 3)
ylabel('x(t)')
% xlabel('t')
% ylim([-2, 3])
title('Current Trajectory', 'FontSize', FS)


% plot([-tau, traj.x], [xh0, traj.y], 'color', cl(2, :), 'LineWidth', 2)

ax2 = nexttile;
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz(tz <= Tstop) + tau, xz(tz <= Tstop), 'Color', cz, 'LineWidth', 3)
scatter(Tstop + tau, xstop, 200, '*k', 'LineWidth', 2)
title('Delayed Trajectory', 'FontSize', FS)
ylabel('x(t-\tau)')
xlim([0, Tmax])

linkaxes([ax1, ax2])
xlabel('time')


%% Stacked Stopped Components
figure(2)
clf
hold on 

tiledlayout(2, 1)
ax1 = nexttile;
hold on
plot(tz(tz <= Tstop), xz(tz <= Tstop), 'Color', cz, 'LineWidth', 3)
scatter(Tstop, xstop, 200, '*k', 'LineWidth', 2)
% plot(tp, xp, 'Color', cp, 'LineWidth', 3)
ylabel('x(t)')
% xlabel('t')
% ylim([-2, 3])
title('Current Trajectory', 'FontSize', FS)


% plot([-tau, traj.x], [xh0, traj.y], 'color', cl(2, :), 'LineWidth', 2)

ax2 = nexttile;
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz(tz <= Tstop-tau) + tau, xz(tz <= Tstop-tau), 'Color', cz, 'LineWidth', 3)
plot(tz((tz <= Tstop) & (tz >= Tstop-tau)) + tau, xz((tz <= Tstop) & (tz >= Tstop-tau)), ':', 'Color', cz, 'LineWidth', 3)
scatter(Tstop + tau, xstop, 200, '*k', 'LineWidth', 2)
scatter(Tstop, xstop_delay, 200, 'ok', 'LineWidth', 2)
title('Delayed Trajectory (Stopped)', 'FontSize', FS)
ylabel('x(t-\tau)')
xlim([0, Tmax])

linkaxes([ax1, ax2])



xlabel('time')


% ylabel('x(t)')
% title('Same Initial Point', 'FontSize', 14)

% %% Components in time
% figure(2)
% clf
% hold on
% plot(tn, xn, 'Color', cn, 'LineWidth', 3)
% plot(tz, xz, 'Color', cz, 'LineWidth', 3)
% plot(tp, xp, 'Color', cp, 'LineWidth', 3)
% plot([0, 0], ylim, ':k', 'LineWidth', 2)
% ylabel('x(t)')
% xlabel('t')
% title('Component Decomposition', 'FontSize', FS)

% %% Delay Embedding
% figure(3)
% clf
% hold on
% shift = 1.5;
% plot3(tz, xz, -shift+zeros(size(tz)), 'color', cz, 'LineWidth', 2)
% plot3(tp, xp, -shift+zeros(size(tp)), 'color', cp, 'LineWidth', 2)
% plot3(tz+tau, shift+zeros(size(tz)), xz, 'color', cz, 'LineWidth', 2)
% plot3(tn+tau, shift+zeros(size(tn)), xn, 'color', cn, 'LineWidth', 2)
% 
% 
% 
% %now plot the delay embedding
% N = 300;
% t = linspace(0, Tmax, N); 
% x0 = deval(traj, t);
% % x1 = f(t - tau);
% x1 = zeros(size(x0));
% x1(t>=tau) = deval(traj, t(t>=tau)-1);
% x1(t < tau) = exp_history_curr(tau- t(t<tau));
% plot3(t, x0, x1, 'k', 'LineWidth', 4)
% view(3)
% zlabel('x(t-\tau)', 'FontSize', FS)
% ylabel('x(t)', 'FontSize', FS)
% xlabel('t', 'FontSize', FS)
% title('Delay Embedding', 'FontSize', 1.5*FS)
% pbaspect([8 3 3])