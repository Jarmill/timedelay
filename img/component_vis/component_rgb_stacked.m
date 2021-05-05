%visualize component measures

%position of plot in monitor
pos = [744   628   691   422];
pos_half = [744   628   691   422/2];

%colors from linspecer
% cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cn = [0,0,1];
cz = [1,0,0];
cp = [0,1,0];
% cz = [0.915294117647059,0.281568627450980,0.287843137254902];
% cp = [0.441568627450980,0.749019607843137,0.432156862745098];

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

%peak
%for this data, the maximum occurs in the 'z' [0, T-tau] = [0, 5] region

[xpeak, ipeak] = max(xz);
tpeak = tz(ipeak);

%% Plotting
FS = 16; %font size of titles


%% marginals plotted independently




%% Stacked Plot
fig = figure(4);
set(fig, 'position', pos);
subplot(2, 1, 1)
clf
hold on
plot(tz, xz, 'Color', cz, 'LineWidth', 3)
plot(tp, xp, 'Color', cp, 'LineWidth', 3)
ylabel('x(t)')
xlabel('t')
ylim([-2, 3])
title('Current Trajectory', 'FontSize', FS)

% fig = figure(5);
set(fig, 'position', pos_half);
clf
subplot(2, 1, 2)
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz + tau, xz, 'Color', cz, 'LineWidth', 3)
title('Delayed Trajectory', 'FontSize', FS)
ylabel('x(t-\tau)')
ylim([-2, 3])
xlabel('t')
