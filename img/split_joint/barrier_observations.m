%Sample trajectories from Section IV of https://web.mit.edu/~jadbabai/www/papers/CDCECC05_2434_MS.pdf

SAMPLE = 1;
PLOT = 1;

rng(43, 'twister')

% T = 20;
T = 5;
% tau = 0.5;
% tau = 0.75;
tau = 1;

C0 = [0; 1];
R0 = 0.1;

Ru = 0.1;
Cu = [0; 0.5];

% R0 = 1;
% C0 = [0; 0];


sample_x = @() R0*ball_sample(1, 2) + C0';

gamma = 0.9;
% f_ode = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];


f_dde = @(t, x, Z) [gamma*x(2) + (1-gamma)*Z(2); -gamma*x(1) - (1-gamma)*Z(1) - 1.5*x(2)];

%start with constant history

%history sampler
box_lim = 4;
Npts = 10;
box_event = @(t, x, Z) supp_event(box_lim, t, x, Z);
dde_options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], ...
                     'MaxStep', 0.1, 'Events', box_event);


%% Sample Trajectories
if SAMPLE
% Nsample = 10;
% Nsample = 200;
Nsample = 400;
out_dde = cell(Nsample, 1);
out_dde_pinhole = cell(Nsample, 1);
for i = 1:Nsample
    out_dde{i} = sample_hist_dde(@history_handle, f_dde, sample_x, Npts, T, tau, C0);
    out_dde_pinhole{i} = sample_hist_dde(@history_handle_pinhole, f_dde, sample_x, Npts, T, tau, C0);    
end
save('barrier_ex_multi.mat', 'T', 'tau', 'out_dde')

end
% [f_history, t_pts, x_pts] = history_handle(sample_x, tau, Npts);


%% Plot Trajectories

if PLOT
    
    if ~SAMPLE
        load('flow_multi_sample_choices.mat')
        Nsample= length(out_dde);
    end
        
figure(50)
clf

tiledlayout(1, 2)

%% jumping plot
nexttile;
hold on
% plot(x_ode(:, 1), x_ode(:, 2))
% i = 20;
% i = 300;
for i = 1:Nsample
    plot(out_dde{i}.xhist(:, 1), out_dde{i}.xhist(:, 2), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(out_dde{i}.x(:, 1), out_dde{i}.x(:, 2), 'c');
end

theta = linspace(0, 2*pi, 100);
circ = R0 * [cos(theta); sin(theta)] + C0;


circ_u = Ru * [cos(theta); sin(theta)] + Cu;

plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
patch(circ_u(1, :), circ_u(2, :), 'r', 'EdgeColor', 'none')
% plot(t_dde, x_dde)
 
sol.obj_rec = 1.162953787;

xlim([-0.5, 0.75])
ylim([-0.75, 1.25])
% plot(xlim, -sol.obj_rec*[1,1], ':r', 'LineWidth', 2)
pbaspect([diff(xlim), diff(ylim), 1])
axis off

%% pinhole plot
nexttile;
hold on
% plot(x_ode(:, 1), x_ode(:, 2))
% i = 20;
% i = 300;
for i = 1:Nsample
    plot(out_dde_pinhole{i}.xhist(:, 1), out_dde_pinhole{i}.xhist(:, 2), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(out_dde_pinhole{i}.x(:, 1), out_dde_pinhole{i}.x(:, 2), 'c');
end

theta = linspace(0, 2*pi, 100);
circ = R0 * [cos(theta); sin(theta)] + C0;


circ_u = Ru * [cos(theta); sin(theta)] + Cu;

plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
patch(circ_u(1, :), circ_u(2, :), 'r', 'EdgeColor', 'none')
% plot(t_dde, x_dde)

xlim([-0.5, 0.75])
ylim([-0.75, 1.25])
% plot(xlim, -sol.obj_rec*[1,1], ':r', 'LineWidth', 2)
pbaspect([diff(xlim), diff(ylim), 1])
axis off

% xlabel('$x_1$', 'interpreter', 'latex')
% ylabel('$x_2$', 'interpreter', 'latex')
% title(['Order 5 bound: ', num2str(-sol.obj_rec)], 'FontSize', 16)

end

%% History Functions
function [hh, t_pts, x_pts] = history_handle(sampler, tau, Npts, C0)
%arbitrarily jump within cylinder
t0 = linspace(-tau, 0, Npts)';
t_pts = t0;
t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*4))*(2*rand(Npts-2,1)-1);
t_pts = sort(t_pts);
x_pts = zeros(Npts, 2);
for i = 1:Npts-1
    x_pts(i, :) = sampler();
end
x_pts(Npts, :) = x_pts(Npts-1, :);
hh = @(t) zoh(t_pts, x_pts, t);

end

function [hh, t_pts, x_pts] = history_handle_pinhole(sampler, tau, Npts, C0)
%jump within cylinder and then go to the center

t0 = linspace(-tau, 0, Npts)';
t_pts = t0;
t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*4))*(2*rand(Npts-2,1)-1);
t_pts = sort(t_pts);
x_pts = zeros(Npts, 2);
for i = 1:Npts-1
    x_pts(i, :) = sampler();
end
x_pts(Npts, :) = C0';
hh = @(t) zoh(t_pts, x_pts, t);

end

function [hh, t_pts, x_pts] = history_handle_const(sampler, tau, Npts, C0)
%history is constant
t_pts = [-tau, 0]';

xchoice = sampler();
x_pts = [xchoice; xchoice];
hh = @(t) xchoice';

end

function [hh, t_pts, x_pts] = history_handle_constjump(sampler, tau, Npts, C0)
%history is set to C0 until a jump within a disk
t_pts = [-tau, 0]';

xchoice = sampler();
x_pts = [C0, xchoice']';
hh = @(t) xchoice;

end


%% Sampling Functions

function out_dde = sample_hist_dde(hist_handle, f_dde, sample_x, Npts, T, tau, C0)
    [f_history, t_pts, x_pts] = hist_handle(sample_x, tau, Npts, C0);
    dde_options.Jumps = t_pts;
    sol_dde = dde23(f_dde, tau, f_history, [0, T], dde_options);
    out_dde.t= sol_dde.x;
    out_dde.x= sol_dde.y';
    out_dde.thist = t_pts;
    out_dde.xhist = x_pts;
end

function [event_eval, terminal, direction] = supp_event(boxlim, t, x, Z)
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