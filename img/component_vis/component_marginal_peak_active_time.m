%visualize component measures

%position of plot in monitor
% pos = [744   628   691   422];
pos = [744   590   824   460];
pos_half = [744   590   824   460/2];

%colors from linspecer
cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cz = [0.915294117647059,0.281568627450980,0.287843137254902];
cp = [0.441568627450980,0.749019607843137,0.432156862745098];

T = 6;
tau = 1;
N = 400;

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
fig = figure(14);
set(fig, 'position', pos);
clf
subplot(2, 1, 1)
hold on
plot(tz(1:ipeak), xz(1:ipeak), 'Color', cz, 'LineWidth', 3)
% plot(tz(ipeak:end), xz(ipeak:end), ':', 'Color', cz, 'LineWidth', 3)
% plot(tp, xp, ':', 'Color', cp, 'LineWidth', 3)
scatter(tpeak, xpeak, 200, '*k', 'LineWidth', 2)
ylabel('x(t)')
% xlabel('t')
title(['Current Elapsed Time = ', num2str(tpeak, 3)], 'FontSize', FS)
xlim([0, T])

subplot(2, 1, 2)
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz(1:ipeak) + tau, xz(1:ipeak), 'Color', cz, 'LineWidth', 3)
% plot(tz(ipeak:end) + tau, xz(ipeak:end), ':', 'Color', cz, 'LineWidth', 3)
scatter(tpeak + tau, xpeak, 200, '*k', 'LineWidth', 2)
% title('Delayed Trajectory', 'FontSize', FS)
title(['Delayed Elapsed Time = ', num2str(tpeak + tau, 3)], 'FontSize', FS)
ylabel('x(t-\tau)')
xlabel('t')
xlim([0, T])

%% with absolute continuity
fig = figure(15);
set(fig, 'position', pos);
clf
subplot(2, 1, 1)
hold on
plot(tz(1:ipeak), xz(1:ipeak), 'Color', cz, 'LineWidth', 3)
% plot(tz(ipeak:end), xz(ipeak:end), ':', 'Color', cz, 'LineWidth', 3)
% plot(tp, xp, ':', 'Color', cp, 'LineWidth', 3)
scatter(tpeak, xpeak, 200, '*k', 'LineWidth', 2)
ylabel('x(t)')
xlim([0, T])
% xlabel('t')
title(['Current Elapsed Time = ', num2str(tpeak, 3)], 'FontSize', FS)


[~, ispeak] = min(abs(tz -(tpeak - tau)));
tspeak = tz(ispeak);
xspeak = xz(ispeak);

subplot(2, 1, 2)
hold on
plot(tn + tau, xn, 'Color', cn, 'LineWidth', 3)
plot(tz(1:ispeak) + tau, xz(1:ispeak), 'Color', cz, 'LineWidth', 3)
% plot(tz(ispeak:ipeak) + tau, xz(ispeak:ipeak), '-.k', 'LineWidth', 3)
% plot(tz(ipeak:end) + tau, xz(ipeak:end), ':', 'Color', cz, 'LineWidth', 3)
scatter(tpeak, xspeak, 200, 'ok', 'LineWidth', 2)
% scatter(tpeak + tau, xpeak, 200, '*k', 'LineWidth', 2)
% title('Delayed Trajectory', 'FontSize', FS)
title(['Delayed Elapsed Time = ', num2str(tpeak , 3)], 'FontSize', FS)
ylabel('x(t-\tau)')
xlabel('t')
xlim([0, T])
