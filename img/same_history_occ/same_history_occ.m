% system
K0 = 2;
K1 = 3;
dydt = @(t, y, z)  -K0*y - K1*z;
tau = 1;
T = 5;

%histories are either +1 or -1
sol = cell(4,1);
history_f = {@(t) history_pure_p(t, tau),...
    @(t) history_pure_n(t, tau),...
    @(t) history_mix_p(t, tau),...
    @(t) history_mix_n(t, tau)};
% sol{1} = dde23(dydt, tau, @(t) history_pure_p(t,tau), [0, T]);
% sol{2} = dde23(dydt, tau, -1, [0, T]);

%% sample the trajectories 
N = 3000;
tsample = linspace(0, T, N);
thsample = linspace(-tau, 0, N/20);


xh = cell(4,1);
x_traj = cell(4, 1);
for i = 1:4    
    sol{i} = dde23(dydt, tau, @(t) history_f{i}(t), [0, T]);
    xh{i} = history_f{i}(thsample);
    x_traj{i} = deval(sol{i}, tsample);
end

%% perform plots

c = linspecer(2);

figure(100)
clf
ax1 = subplot(2,1,1);
hold on
for i = 1:2
    c_curr = c(i, :);
    plot(thsample, xh{i}, 'color', c_curr, 'linewidth', 2);
    plot(tsample, x_traj{i}, 'color', c_curr,'linewidth', 2);
end
ylabel('x(t)', 'fontsize', 12)
xlabel('t', 'fontsize', 12)
title('Pure History', 'fontsize', 16)

ax2 = subplot(2,1,2);
hold on
for i = 3:4
    c_curr = c(i-2, :);
    plot(thsample, xh{i}, 'color', c_curr, 'linewidth', 2);
    plot(tsample, x_traj{i}, 'color', c_curr,'linewidth', 2);
end

linkaxes([ax1, ax2])
ylabel('x(t)', 'fontsize', 12)
xlabel('t', 'fontsize', 12)
title('Mixed History', 'fontsize', 16)
%% history functions

function xh = history_pure_p(t, tau)
    xh = ones(size(t));
end

function xh = history_pure_n(t, tau)
    xh = -ones(size(t));
end

function xh = history_mix_p(t,tau)
    xh = ones(size(t));
    xh(t>=(-tau/2)) = -1;
end


function xh = history_mix_n(t,tau)
    xh = ones(size(t));
    xh(t<=(-tau/2)) = -1;
end
