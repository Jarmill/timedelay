%create three trajectories x1(t), x2(t), x3(t) of a 1d DDE system.
%form occupation measures by mixing trajectories, such as (t, x1(t),
%x2(t-tau)). This will obey the consistency constraints. Then see if the
%Liouville equation remains valid

%% generate weights for mixing measures


%% sample trajectories

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);

x0 = [2; 1; -1];

K0 = 2;
K1 = 3;
dydt = @(t, y, z)  -K0*y - K1*z;
tau = 1;
T = 1;
sol = cell(2, 1);
u = [-1; 0.5];
for i = 1:2
    sol{i} = ode23(@(t, y) dydt(t, y, u(i)), [0, T], x0(i));
end
sol_op = cell(2, 1);
for i = 1:2
    sol_op{i} = ode23(@(t, y) dydt(t, y, u(3-i)), [0, T], x0(i));
end


% N = 400;
% N = 2000;
N = 3000;
tsample = linspace(0, T, N);
% tn = linspace(-tau, 0, N/5);
% tz = linspace(0, T - tau, N*4/5);
% tp = linspace(T-tau, T, N/5);
x_traj  = zeros(2, N); %trajectories (x0(i), u(i))
x_op  = zeros(2, N);   %trajectories (x0(i), u(3-i))
% x_delay = zeros(3, N);

for i = 1:2
    x_traj(i, :) = deval(sol{i}, tsample);
    x_op(i, :) = deval(sol_op{i}, tsample);
end



%% evaluate monomial moments along mixed trajectories
d = 2;
dv = genPowGlopti(2,d);
% m_traj = zeros(3, 3, length(dv));
% for i = 1:3
%     for j = 1:3
%         m_traj(i, j, :) = monom_int(tsample, [x_traj(i, :); x_delay(j, :)], dv);
%     end
% end


% dv(2, :)
% 
% x0_monom = x0.^(dv(:, 2)')';
% xT_monom = T.^dv(:, 1).* (x_traj(:, end).^(dv(:, 2)')');
% 
% % mix_traj = m_traj.*mix;
% row_traj = squeeze(sum(mix_traj, 1))'./side_weights';
% col_traj = squeeze(sum(mix_traj, 2))'./side_weights';
% diag_traj = squeeze([m_traj(1, 1, :), m_traj(2, 2, :), m_traj(3, 3, :)])';

%% test liouville evaluation
dv_curr = [1; 1; 2];
i_curr = 50;

% dydt(tsample(i_curr), x_traj(i_curr), x_delay(i_curr))*

%% plot results


figure(1)
clf
subplot(2, 1, 1)
hold on
for i = 1:2
plot(tsample, x_traj(i, :))
end
title('Current time')
subplot(2, 1, 2)
hold on
for i = 1:2
plot(tsample, ones(N, 1)*u(i))
end
title('Input')

%% 3d plot of only sampled trajectories
figure(3)
clf
hold on
for i = 1:2
    plot3(tsample, x_traj(i, :), u(i)*ones(N, 1), 'LineWidth', 3)
end
view(3)
xlabel('t', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
zlabel('u(t)', 'interpreter', 'latex')
title('Input Embedding')


figure(4)
clf
hold on
for i = 1:2
    plot3(tsample, x_traj(i, :), u(i)*ones(N, 1), 'LineWidth', 3)
end

for i = 1:2
    for j = 1:(i-1)
        plot3(tsample, x_traj(i, :), u(j)*ones(N, 1), 'k',  'LineWidth', 1)
        plot3(tsample, x_traj(j, :), u(i)*ones(N, 1), 'k',  'LineWidth', 1)
    end
end
view(3)
xlabel('t', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
zlabel('u(t)', 'interpreter', 'latex')
title('Input Embedding (marginal confused)')

%% with opposite inputs
figure(5)
% figure(3)
clf
hold on
for i = 1:2
    plot3(tsample, x_traj(i, :), u(i)*ones(N, 1), 'LineWidth', 3)
end

for i = 1:2
    plot3(tsample, x_op(i, :), u(3-i)*ones(N, 1), 'LineWidth', 3)
end
view(3)
xlabel('t', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
zlabel('u(t)', 'interpreter', 'latex')
title('Input Embedding (with opposite input)')
