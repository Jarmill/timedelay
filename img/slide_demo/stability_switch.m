%show a set of trajectories with constant histories

rng(20, 'twister')

SAMPLE = 1;
PLOT = 1;

K0 = 0;
K1 = 1;
xh0 = 1;

%system is stable for tau < pi/2
%unstable for tau > pi/2

hlim = [0.5, 1.5];
Ntau = 5;

tau_low = linspace(1e-3, pi/2-1e-3, Ntau);
tau_high = linspace(pi/2+1e-3, pi, Ntau);

tau = 1;
Tmax = 10;
% Npts = 3;

exp_delay = @(t,y,Z) -K0*y - K1*Z;
exp_history = @(t) xh0;



sample_x = @() diff(hlim)*(2*rand()-1) + mean(hlim);

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'MaxStep', 0.1, 'Jumps', [0]);
if SAMPLE
    out_sim_low = cell(Ntau, 1);
    out_sim_high = cell(Ntau, 1);
    
    for i = 1:Ntau
        out_sim_low{i}= dde23(exp_delay, tau_low(i), exp_history, [0,Tmax], options);            
        out_sim_high{i}= dde23(exp_delay, tau_high(i), exp_history, [0,Tmax], options);            
        
    end
end

cl = linspecer(1);

%% Plot
if PLOT

cl = linspecer(Ntau*2);
figure(1)
clf
tiledlayout(2, 1)
ax1 = nexttile;
hold on 



for i = Ntau:-1:1
    plot(ax1,[-tau_low(i), 0], [1,1]*xh0, 'color', cl(i, :), 'LineWidth', 2)
    plot(ax1,out_sim_low{i}.x, out_sim_low{i}.y, 'color', cl(i, :), 'LineWidth', 2)
end
plot([0, 0], ylim, ':k', 'LineWidth', 2)
plot(-pi/2*[1, 1], ylim, ':r', 'LineWidth', 2)
xlabel('time')
ylabel('x(t)')
title('Stable, \tau < \pi/2', 'FontSize', 18)

ax2 = nexttile;
hold on 


for i = Ntau:-1:1
    plot(ax2,[-tau_high(i), 0], [1,1]*xh0, 'color', cl(Ntau+i, :), 'LineWidth', 2)
    plot(ax2,out_sim_high{i}.x, out_sim_high{i}.y, 'color', cl(Ntau+i, :), 'LineWidth', 2)
end

xlabel('time')
ylabel('x(t)')
title('Unstable, \tau > \pi/2', 'FontSize', 18)

linkaxes([ax1, ax2])
plot(ax1, [0, 0], ylim, ':k', 'LineWidth', 2)
plot(ax1, -pi/2*[1, 1], ylim, ':r', 'LineWidth', 2)

plot(ax2, [0, 0], ylim, ':k', 'LineWidth', 2)
plot(ax2, -pi/2*[1, 1], ylim, ':r', 'LineWidth', 2)
xlim([-pi, Tmax])
end