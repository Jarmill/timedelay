SAMPLE = 1;
PLOT = 1;

rng(43, 'twister')

% T = 20;
T = 6;
% T = 5;
% tau = 0.5;
% tau = 0.75;
tau = 1;
tau2 = 2;

C0 = -1;
R0 = 0.2;

% R0 = 1;
% C0 = [0; 0];


sample_x = @() R0*ball_sample(1, 1) + C0';


f_ode = @(t, x) - 0.5*x(1);

f_dde = @(t, x, Z) 0.5*x(1) - Z(1);

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
    Nsample = 500;
    out_dde = cell(Nsample, 1);
    out_dde_const = cell(Nsample, 1);
    for i = 1:Nsample
        %jumping history
        [f_history, t_pts, x_pts] = history_handle(sample_x, tau2, Npts);
        dde_options.Jumps = t_pts;
        sol_dde = dde23(f_dde, [tau, tau2], f_history, [0, T], dde_options);
        out_dde{i}.t= sol_dde.x';
        out_dde{i}.x= sol_dde.y';
        out_dde{i}.thist = t_pts;
        out_dde{i}.xhist = x_pts;

        %constant history
        [f_history, t_pts, x_pts] = history_handle(sample_x, tau2, 1);        
        dde_options.Jumps = t_pts;
        sol_dde_const = dde23(f_dde, [tau, tau2], f_history, [0, T], dde_options);
        out_dde_const{i}.t= sol_dde_const.x';
        out_dde_const{i}.x= sol_dde_const.y';
        out_dde_const{i}.thist = t_pts;
        out_dde_const{i}.xhist = x_pts;
    end
end
% [f_history, t_pts, x_pts] = history_handle(sample_x, tau, Npts);


%% Plot Trajectories

if PLOT
    
    if ~SAMPLE
        load('univariate_multi_sample.mat')
        Nsample= length(out_dde);
    end
        
figure(50)
clf
tiledlayout(1, 2)
nexttile
hold on
% plot(x_ode(:, 1), x_ode(:, 2))

for i = 1:Nsample
    plot(out_dde{i}.thist(:, 1), out_dde{i}.xhist(:, 1), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(out_dde{i}.t(:, 1), out_dde{i}.x(:, 1), 'c');
end
 
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')
% title(['Order 5 bound: ', num2str(-sol.obj_rec)], 'FontSize', 16)
title('Jumping Histories', 'FontSize', 16)


nexttile
hold on
% plot(x_ode(:, 1), x_ode(:, 2))

for i = 1:Nsample
    plot(out_dde_const{i}.thist(:, 1), out_dde_const{i}.xhist(:, 1), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(out_dde_const{i}.t(:, 1), out_dde_const{i}.x(:, 1), 'c');
end
 
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')
title('Constant Histories', 'FontSize', 16)
% title(['Order 5 bound: ', num2str(-sol.obj_rec)], 'FontSize', 16)

end

%% Function
function [hh, t_pts, x_pts] = history_handle(sampler, tau, Npts)    

    if Npts > 1
        t0 = linspace(-tau, 0, Npts)';
        t_pts = t0;
        t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*4))*(2*rand(Npts-2,1)-1);
        t_pts = sort(t_pts);
        x_pts = zeros(Npts, 1);
        for i = 1:Npts-1
            x_pts(i, :) = sampler();
        end
        x_pts(Npts, :) = x_pts(Npts-1, :);
        hh = @(x) zoh(t_pts, x_pts, x);
    else
        xc = sampler();
        t_pts = [-tau; 0];
        x_pts = [xc; xc];
        hh = @(x) xc;
    end
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