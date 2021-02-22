%visualize component measures


%colors from linspecer
cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cz = [0.915294117647059,0.281568627450980,0.287843137254902];
cp = [0.441568627450980,0.749019607843137,0.432156862745098];

T = 6;
tau = 1;
N = 300;

f = @(x) 1*(0.5 * cos(3*x) -  sin(5.5*x) + 0.3*sin(5*x) + sin(0.1*x) - cos(10*x));

%full time
t = linspace(0, T, N); 
x0 = f(t);
x1 = f(t - tau);

%components
tn = linspace(-tau, 0, N/4);
tz = linspace(0, T - tau, N);
tp = linspace(T-tau, T, N/4);

xn = f(tn);
xz = f(tz);
xp = f(tp);

%% Plotting
FS = 14; %font size of titles

figure(1)
clf
hold on
plot3(tz, xz, -3+zeros(size(tz)), 'color', cz)
plot3(tp, xp, -3+zeros(size(tp)), '-.','color', cp)
plot3(tz+tau, 3+zeros(size(tz)), xz, 'color', cz)
plot3(tn+tau, 3+zeros(size(tn)), xn, ':', 'color', cn)
plot3(t, x0, x1, 'k', 'LineWidth', 2)

view(3)
zlabel('x(t-\tau)')
ylabel('x(t)')
xlabel('t')
title('Delay Embedding', 'FontSize', FS)
pbaspect([6 3 3])
figure(3)
clf
hold on
plot(tn, xn,  ':','Color',  cn, 'LineWidth', 3)
plot(tz, xz, 'Color', cz, 'LineWidth', 3)
plot(tp, xp, '-.','Color', cp, 'LineWidth', 3)
% plot([0, 0], ylim, ':k', 'LineWidth', 2)
ylabel('x(t)')
xlabel('t')
title('Full Trajectory', 'FontSize', FS)


figure(4)
clf
subplot(2, 1, 1)
hold on
plot(tz, xz, 'Color', cz, 'LineWidth', 3)
plot(tp, xp, '-.', 'Color', cp, 'LineWidth', 3)
ylabel('x(t)')
xlabel('t')
title('Current Trajectory', 'FontSize', FS)

subplot(2, 1, 2)
hold on
plot(tn + tau, xn,  ':', 'Color', cn,'LineWidth', 3)
plot(tz + tau, xz, 'Color', cz, 'LineWidth', 3)
title('Delayed Trajectory', 'FontSize', FS)
ylabel('x(t-\tau)')
xlabel('t')