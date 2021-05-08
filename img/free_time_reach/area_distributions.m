%plots of t-marginals associated with terminal and occupation measures of
%free-time reachable set systems

%colors from linspecer
cn = [0.363921568627451,0.575529411764706,0.748392156862745];
cz = [0.915294117647059,0.281568627450980,0.287843137254902];
cp = [0.441568627450980,0.749019607843137,0.432156862745098];


rng(30, 'twister')

T = 5;

%define the cross-sectional areas
alim = [1; 5];
sample_x = @() diff(alim)*(2*rand()-1)/2 + mean(alim);

Npts = 15;
[hh, t_pts, x_pts, pp] = area_handle(sample_x, T, Npts);

%hh(time) = cross-sectional area of reachable set at time t

%sample the areas
Nsample = 600;
t = linspace(0, T, Nsample);
x = hh(t);

%simpson's rule is exact for polynomials of degree 3 or less


V = simps(t, x); %volume
mass_term = V;

x_normalized = x/V;
%mean terminal time of this process
%equal to the mass of the occupation measure/mass of terminal measure

mass_occ = simps(t, x.*t);
avg_time = mass_occ/mass_term;

% avg_time = simps(t, x_normalized.*t);
% mass_occ = mass_term * avg_time;



%find the cumulative integral of area from T back to 0

%probability that a trajectory starting at X0 and ending in XT will pass
%through time [tmin, tmax]
% x_prob_occ = 1 - cumsimps(t, x_normalized);
x_occ = V - cumsimps(t, x);

%% one time delay
tau = 2;
t_lag = linspace(0, T-tau, Nsample*(T-tau)/T);
x_lag = hh(t_lag);
x_lag_occ = V - cumsimps(t_lag, x_lag);

t_history = linspace(0, tau, Nsample*tau/T);
x_history_occ = mass_term * ones(size(t_history));

t_delay = [t_history, tau + t_lag];
x_delay_occ = [x_history_occ, x_lag_occ];

x_slack = x_delay_occ - x_occ;


ind_z = 1:((T-tau)/T)*Nsample;
ind_p = (1:((tau/T)*Nsample)) + ind_z(end);


tn = t_history;
tz = t(ind_z);
tp = t(ind_p);
xn = x_history_occ;
xz = x_occ(ind_z);
xp = x_occ(ind_p);


%% Plots

FS = 16;

figure(1)
clf
tiledlayout(4, 1)
nexttile
plot(t, x, 'LineWidth', 3)
title('free-time terminal measure', 'fontsize', FS)
xlabel('time')
ylabel('cross sectional area')

ax2 = nexttile;
% plot(t, x_occ)
hold on

plot(tz, xz, 'Color', cz, 'LineWidth', 3)
plot(tp, xp, 'Color', cp, 'LineWidth', 3)
title('component measure current time', 'fontsize', FS)
xlabel('time')
ylabel('density')

ax3 = nexttile;
% plot(t_delay, x_delay_occ);
hold on

plot(tn, xn, 'Color', cn, 'LineWidth', 3)
plot(tz+tau, xz, 'Color', cz, 'LineWidth', 3)
% plot(tp, xp, 'Color', cp, 'LineWidth', 3)
title('component measure delayed', 'fontsize', FS)
xlabel('time')
ylabel('density')

ax4 = nexttile;
plot(t_delay, x_slack, 'Color', 'k', 'LineWidth', 3);
title('free-time slack measure', 'fontsize', FS)
xlabel('time')
ylabel('density')


linkaxes([ax2, ax3,ax4])

figure(2)
clf
x_rad = sqrt(x)/pi;
[X, Y, Z] = cylinder(x_rad, 40);
surf(T*Z, X, Y)
title('Reachable Set Example (surface of revolution)', 'fontsize', FS)
xlabel('time')
ylabel('x_1')
zlabel('x_2')
function [hh, t_pts, x_pts, pp] = area_handle(sampler, T, Npts)

t0 = linspace(0, T, Npts);
t_pts = t0;
t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*2))*(2*rand(1, Npts-2)-1);

% t_pts = rand(Npts-2, 1)*(-tau);
% t_pts = [-tau; sort(t_pts); 0];
x_pts = zeros(Npts, 1);
for i = 1:Npts-1
    x_pts(i) = sampler();
end
x_pts(Npts) = x_pts(Npts-1);
% hh = @(t) zoh(t_pts, x_pts, t);

pp = spline(t_pts, x_pts);
hh = @(t) ppval(pp, t);

end