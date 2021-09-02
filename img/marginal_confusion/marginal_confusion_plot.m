%create three trajectories x1(t), x2(t), x3(t) of a 1d DDE system.
%form occupation measures by mixing trajectories, such as (t, x1(t),
%x2(t-tau)). This will obey the consistency constraints. Then see if the
%Liouville equation remains valid

%% generate weights for mixing measures

side_weights = [0.5; 0.3; 0.2];
diag_weights = [0.3; 0.15; 0.05];

A_row = kron(eye(3), [1,1,1]);
A_col = kron([1,1,1], eye(3));
A_diag = full(sparse([1,2,3], [1, 5, 9], [1,1,1]));

A = [A_row; A_col; A_diag];

% A = [kron(eye(3), [1,1,1]); ...
%     kron([1,1,1], eye(3)); ...
%     full(sparse([1,2,3], [1, 5, 9], [1,1,1]))];
b =  [side_weights; side_weights; diag_weights];

f = ones(9, 1);
mix = linprog(f, -eye(9), zeros(9, 1), A, b);
mix = reshape(mix, 3, 3);

%% sample trajectories

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);

x0 = [2; 1; -1];

K0 = 2;
K1 = 3;
dydt = @(t, y, z)  -K0*y - K1*z;
tau = 1;
T = 5;
sol = cell(3, 1);
for i = 1:3
    sol{i} = dde23(dydt, tau, @(t) x0(i), [0, T]);
end


% N = 400;
% N = 2000;
N = 3000;
% tsample = linspace(0, T, Nsample);
% tn = linspace(-tau, 0, N/5);
% tz = linspace(0, T - tau, N*4/5);
% tp = linspace(T-tau, T, N/5);
tsample = linspace(0, T, N);
tz = tsample(tsample <= T-tau);
tp = tsample(tsample >= T-tau);
x_traj  = zeros(3, N);
x_delay = zeros(3, N);

for i = 1:3
    x_traj(i, :) = deval(sol{i}, tsample);
    x_delay(i, (N/5)+1:end) = deval(sol{i}, tz);
    x_delay(i, 1:(N/5)) = x0(i);
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
for i = 1:3
plot(tsample, x_traj(i, :))
end
title('Current time')
subplot(2, 1, 2)
hold on
for i = 1:3
plot(tsample, x_delay(i, :))
end
title('Delayed time')

%% 3d plot of only sampled trajectories
figure(3)
clf
hold on
for i = 1:3
    plot3(tsample, x_traj(i, :), x_delay(i, :), 'LineWidth', 3)
end
view(3)
xlabel('t', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
zlabel('x(t-1)', 'interpreter', 'latex')
title('Delay Embedding')


figure(4)
clf
hold on
for i = 1:3
    plot3(tsample, x_traj(i, :), x_delay(i, :), 'LineWidth', 3)
end

for i = 1:3
    for j = 1:(i-1)
        plot3(tsample, x_traj(i, :), x_delay(j, :), 'k',  'LineWidth', 1)
        plot3(tsample, x_traj(j, :), x_delay(i, :), 'k',  'LineWidth', 1)
    end
end
view(3)
xlabel('t', 'interpreter', 'latex')
ylabel('x(t)', 'interpreter', 'latex')
zlabel('x(t-1)', 'interpreter', 'latex')
title('Delay Embedding')