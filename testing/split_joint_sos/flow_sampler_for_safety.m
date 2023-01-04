SAMPLE = 1;
PLOT = 1;

rng(43, 'twister')

% T = 7;
% T = 5;
T = 8;
% tau = 0.5;
tau = 0.5;

C0 = [1.5; 0];
R0 = 0.4;

Cu = [-0.5; -1];
Ru = 0.5;

% R0 = 1;
% C0 = [0; 0];


sample_x = @() R0*ball_sample(1, 2) + C0';


f_ode = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

f_dde = @(t, x, Z) [x(2); -Z(1) + (1/3).*x(1).^3- x(2)];

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
Nsample = 100;
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
        load('flow_multi_sample.mat')
        Nsample= length(out_dde);
    end
        
figure(50)
clf
hold on
% plot(x_ode(:, 1), x_ode(:, 2))
% i = 20;
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

theta_u = linspace(3*pi/4, 7*pi/4, 100);
arc_u = Ru * [cos(theta_u); sin(theta_u)] + Cu;
patch(arc_u(1, :), arc_u(2, :), 'r', 'EdgeColor', 'none')


%time horizon T=5
% dist_rec = 0.304642524411440;

%time horizon T = 8 
% dist_rec = 0.157260043742384; %order 3.
dist_rec = 0.182053863415561; %order 4
x_dist = dist_contour(100, Ru, dist_rec);
x_dist_rot = [cos(pi/4), sin(pi/4); -sin(pi/4), cos(pi/4)] * x_dist + Cu;

plot(x_dist_rot(1, :), x_dist_rot(2, :), 'r', 'LineWidth', 3)

% xlim([-1.25, 2.5])
% ylim([-1.5, 1.1])

xlim([-1.5, 2.5])
ylim([-1.85, 1.1])
% plot(xlim, -sol.obj_rec*[1,1], ':r', 'LineWidth', 2)
pbaspect([diff(xlim), diff(ylim), 1])

% axis off

xlabel('x_1')
ylabel('x_2')
% title(['Order 4 bound: ', num2str(-sol.obj_rec)], 'FontSize', 16)

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

function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end