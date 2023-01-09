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
out_dde_const = cell(Nsample, 1);
out_dde_constjump = cell(Nsample, 1);
out_dde_pinhole = cell(Nsample, 1);
for i = 1:Nsample
    out_dde{i} = sample_hist_dde(@history_handle, f_dde, sample_x, Npts, T, tau, C0);
    out_dde_const{i} = sample_hist_dde(@history_handle_const, f_dde, sample_x, Npts, T, tau, C0);
    out_dde_constjump{i} = sample_hist_dde(@history_handle_constjump, f_dde, sample_x, Npts, T, tau, C0);
    out_dde_pinhole{i} = sample_hist_dde(@history_handle_pinhole, f_dde, sample_x, Npts, T, tau, C0);    
end
save('time_var_multi_sample_choices.mat', 'T', 'tau', 'out_dde', 'out_dde_const', 'out_dde_constjump', 'out_dde_pinhole')

end
% [f_history, t_pts, x_pts] = history_handle(sample_x, tau, Npts);


%% Plot Trajectories

if PLOT
    
    if ~SAMPLE
        load('time_var_multi_sample_choices.mat')
        Nsample= length(out_dde);
    end
        
traj_cell = {out_dde, out_dde_pinhole, out_dde_const, out_dde_constjump};
traj_names = {'Free', 'Pinhole', 'Constant', 'Constant-Center Jump'};

figure(50)
clf
hold on
tccurr = out_dde_const;
% plot(x_ode(:, 1), x_ode(:, 2))
i = 20;
for i = 1:Nsample
    plot(tccurr{i}.xhist(:, 1), tccurr{i}.xhist(:, 2), 'color', 0.7*[1,1,1]);    
end

for i =1:Nsample
    plot(tccurr{i}.x(:, 1), tccurr{i}.x(:, 2), 'c');
end

theta = linspace(0, 2*pi, 100);
circ = R0 * [cos(theta); sin(theta)] + C0;

plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
% plot(t_dde, x_dde)
 
% sol.obj_rec = 0.71826461805818;
plot(sol.obj_rec*[1,1], ylim, ':r', 'LineWidth', 2)

xlim([-1.1, 0.8])
% ylim([-1.25, 1.5])
pbaspect([diff(xlim), diff(ylim), 1])


xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
% title(['Order 5 bound: ', num2str(sol.obj_rec)], 'FontSize', 16)


%% 3d plot comparison
figure(51)
tiledlayout(2, 2, 'TileSpacing','Compact','Padding','Compact')

tile_cell = cell(4, 1);
yl = [];
zl = [];
for k = 1:4
    tile_cell{k} = nexttile;
    hold on
    tc = traj_cell{k};
    for i = 1:Nsample
        plot3(tc{i}.thist_sample, tc{i}.xhist_sample(:, 1), tc{i}.xhist_sample(:, 2),  'color', 0.7*[1,1,1]);    
    end

    for i =1:Nsample
        plot3(tc{i}.t, tc{i}.x(:, 1), tc{i}.x(:, 2), 'c');
    end
    xlim([-tau, T])
    if k==1
        yl = ylim;
        zl = zlim;
    end
        ylim(yl);
        zlim(zl);
    
    view(-5, 77)
    xlabel('$t$', 'interpreter', 'latex')
    ylabel('$x_1$', 'interpreter', 'latex')
    zlabel('$x_2$', 'interpreter', 'latex')
    tile_cell{k}.YAxis.Visible='off';
    tile_cell{k}.ZAxis.Visible='off';
    title(traj_names{k});
%     title(['Order 5 bound: ', num2str(sol.obj_rec)], 'FontSize', 16)

end

%% show the plot

figure(52)
clf
% tiledlayout(2, 2, 'TileSpacing','Compact','Padding','Compact')

% tile_cell = cell(4, 1);
yl = [];
zl = [];
k=1;
%     tile_cell{k} = nexttile;
    hold on
    tc = traj_cell{k};
    for i = 1:Nsample
        plot3(tc{i}.thist_sample, tc{i}.xhist_sample(:, 1), tc{i}.xhist_sample(:, 2),  'color', 0.7*[1,1,1]);    
    end

    for i =1:Nsample
        plot3(tc{i}.t, tc{i}.x(:, 1), tc{i}.x(:, 2), 'c');
    end
    xlim([-tau, T])
    plot3(-tau*ones(length(circ), 1), circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
    plot3(0*ones(length(circ), 1), circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
    if k==1
        yl = ylim;
        zl = zlim;
    end
        ylim(yl);
        zlim(zl);
    xl =xlim;
    view(-5, 77)
    xlabel('$t$', 'interpreter', 'latex')
    ylabel('$x_1$', 'interpreter', 'latex')
    zlabel('$x_2$', 'interpreter', 'latex')
    tile_cell{k}.YAxis.Visible='off';
    tile_cell{k}.ZAxis.Visible='off';
    title(traj_names{k});

    bnd=0.71826;
    title(['Order 5 bound: ', num2str(bnd)], 'FontSize', 16)
    patch(xl([1,1,2,2,1]), bnd*ones(1,5), zl([1,2,2,1,1]), 'r', 'EdgeColor', 'None', 'FaceAlpha', 0.5)


% linkaxes([tile_cell{1}, tile_cell{2}, tile_cell{3}, tile_cell{4}], 'xyz')

% %% 3d plot
% figure(51)
% clf
% hold on
% % plot(x_ode(:, 1), x_ode(:, 2))
% for i = 1:Nsample
%     plot3(out_dde{i}.thist, out_dde{i}.xhist(:, 1), out_dde{i}.xhist(:, 2),  'color', 0.7*[1,1,1]);    
% end
% 
% for i =1:Nsample
%     plot3(out_dde{i}.t, out_dde{i}.x(:, 1), out_dde{i}.x(:, 2), 'c');
% end
% 
% bnd = [0.718264618058180];
% ylim([-1.1, 0.8])
% zlim([   -0.5000    1.0000])
% 
% xl = xlim;
% zl = ylim;
% patch(xl([1,1,2,2,1]), bnd*ones(1,5), zl([1,2,2,1,1]), 'r', 'EdgeColor', 'None', 'FaceAlpha', 0.5)
% 
% 
% pbaspect([diff(xlim), diff(ylim), diff(zlim)])
% xlabel('$t$', 'interpreter', 'latex')
% ylabel('$x_1$', 'interpreter', 'latex')
% zlabel('$x_2$', 'interpreter', 'latex')
% title(['Order 5 bound: ', num2str(sol.obj_rec)], 'FontSize', 16)


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
hh = @(t) (xchoice'*ones(1, length(t)))';

end

function [hh, t_pts, x_pts] = history_handle_constjump(sampler, tau, Npts, C0)
%history is set to C0 until a jump within a disk
t_pts = [-tau, 0]';

xchoice = sampler();
x_pts = [C0, xchoice']';
% hh = @(t) (C0*ones(1, length(t)))';
hh = @(t) (C0*(t<0))' + (xchoice'*(t==0))';

end


%% Sampling Functions

function [out_dde]= sample_hist_dde(hist_handle, f_dde, sample_x, Npts, T, tau, C0)
    [f_history, t_pts, x_pts] = hist_handle(sample_x, tau, Npts, C0);
    dde_options.Jumps = t_pts;
    sol_dde = dde23(f_dde, tau, f_history, [0, T], dde_options);
    out_dde.t= sol_dde.x;
    out_dde.x= sol_dde.y';
    out_dde.thist = t_pts;
    out_dde.xhist = x_pts;
    
    out_dde.thist_sample = linspace(-tau, 0, 500);
    out_dde.xhist_sample = f_history(out_dde.thist_sample);
%     out_dde.fhist = f_history;
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