SAMPLE = 0;
PLOT = 1;

rng(43, 'twister')

% T = 20;
T = 5;
% tau = 0.5;
tau = 0.75;

%initial set
C0 = [-0.75; 0];
R0 = 0.3;

% R0 = 1;
% C0 = [0; 0];


sample_x = @() R0*ball_sample(1, 2) + C0';


f_ode = @(t, x) [x(2)*t - 0.1*x(1) - x(1)*x(2);
    -x(1)*t - x(2) + x(1)^2];

f_dde = @(t, x, Z) [x(2)*t - 0.1*x(1) - Z(1)*Z(2);
    -x(1)*t - x(2) + x(1)*Z(1)];

%start with constant history
% f_history = @(t) x0;

% ode_options =   odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.01);
% [t_ode, x_ode] = ode45(f_ode, [0, T], x0, ode_options);

%history sampler
box_lim = 3;
Npts = 10;
box_event = @(t, x, Z) supp_event(box_lim, t, x, Z);
dde_options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], ...
                     'MaxStep', 0.1, 'Events', box_event);
% dde_options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], ...
%                      'MaxStep', 0.1);

%% Sample Trajectories
if SAMPLE
Nsample = 150;
out_dde = cell(Nsample, 1);
for i = 1:Nsample
    [f_history, t_pts, x_pts] = history_handle(sample_x, tau, Npts);
    dde_options.Jumps = t_pts;
    sol_dde = dde23(f_dde, tau, f_history, [0, T], dde_options);
    out_dde{i}.t= sol_dde.x;
    out_dde{i}.x= sol_dde.y';
    out_dde{i}.thist = t_pts;
    out_dde{i}.xhist = x_pts;
end
end
% [f_history, t_pts, x_pts] = history_handle(sample_x, tau, Npts);


%% Plot Trajectories

if PLOT
    
    if ~SAMPLE
        load('time_var_multi_sample.mat')
        Nsample= length(out_dde);
    end
        
figure(50)
clf
hold on
% plot(x_ode(:, 1), x_ode(:, 2))
i = 20;
for i = 1:Nsample
    plot(out_dde{i}.xhist(:, 1), out_dde{i}.xhist(:, 2), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(out_dde{i}.x(:, 1), out_dde{i}.x(:, 2), 'c');
end

theta = linspace(0, 2*pi, 100);
circ = R0 * [cos(theta); sin(theta)] + C0;

plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
% plot(t_dde, x_dde)
 
sol.obj_rec = 0.71826461805818;
plot(sol.obj_rec*[1,1], ylim, ':r', 'LineWidth', 2)

xlim([-1.1, 0.8])
% ylim([-1.25, 1.5])
pbaspect([diff(xlim), diff(ylim), 1])


xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title(['Order 5 bound: ', num2str(sol.obj_rec)], 'FontSize', 16)


%% 3d plot
figure(51)
clf
hold on
% plot(x_ode(:, 1), x_ode(:, 2))
for i = 1:Nsample
    plot3(out_dde{i}.thist, out_dde{i}.xhist(:, 1), out_dde{i}.xhist(:, 2),  'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot3(out_dde{i}.t, out_dde{i}.x(:, 1), out_dde{i}.x(:, 2), 'c');
end

bnd = [0.718264618058180];
ylim([-1.1, 0.8])
zlim([   -0.5000    1.0000])

xl = xlim;
zl = ylim;
patch(xl([1,1,2,2,1]), bnd*ones(1,5), zl([1,2,2,1,1]), 'r', 'EdgeColor', 'None', 'FaceAlpha', 0.5)



% theta = linspace(0, 2*pi, 100);
% circ = R0 * [cos(theta); sin(theta)] + C0;

% plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
% plot(t_dde, x_dde)
 
% sol.obj_rec = 1.162953787;
% plot(sol.obj_rec*[1,1], ylim, ':r', 'LineWidth', 2)

% ylim([-1.25, 1.5])
pbaspect([diff(xlim), diff(ylim), diff(zlim)])
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x_1$', 'interpreter', 'latex')
zlabel('$x_2$', 'interpreter', 'latex')
title(['Order 5 bound: ', num2str(sol.obj_rec)], 'FontSize', 16)


end

%% Function
function [hh, t_pts, x_pts] = history_handle(sampler, tau, Npts)

t0 = linspace(-tau, 0, Npts)';
t_pts = t0;
t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*4))*(2*rand(Npts-2,1)-1);
t_pts = sort(t_pts);
x_pts = zeros(Npts, 2);
for i = 1:Npts-1
    x_pts(i, :) = sampler();
end
x_pts(Npts, :) = x_pts(Npts-1, :);
hh = @(x) zoh(t_pts, x_pts, x);

end

function [event_eval, terminal, direction] = supp_event(boxlim , t, x, Z)
    %event function for @ode15 or other solver
    Npt = size(x, 2);
    event_eval = zeros(1, Npt);
    for i = 1:Npt
        xcurr = x(:, i);
%         tcurr = t(:, i);               

        event_eval(i) = all(abs(x) <= boxlim);
    end

    %stop integrating when the system falls outside support

    terminal = 1;
    direction = 0;                        
end