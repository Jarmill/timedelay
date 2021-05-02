%show a set of trajectories with constant histories

rng(20, 'twister')

SAMPLE = 1;
PLOT = 1;

% K0 = 0;
% K1 = 1;
K0 = 2;
K1 = 2;
xh0 = 1;

hlim = [0.5, 1.5];
Nh = 1;

xh0_list = linspace(hlim(1), hlim(2), Nh);

tau = 1;
Tmax = 5;
% Npts = 3;

exp_delay = @(t,y,Z) -K0*y - K1*Z;
% exp_history = @(t) xh0;



sample_x = @() diff(hlim)*(2*rand()-1) + mean(hlim);

Npts_h = 5;
[exp_history, t_pts, x_pts]= history_handle(sample_x, tau, Npts_h);

t_h = linspace(-tau, 0, 100);
x_h = exp_history(t_h);
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'MaxStep', 0.1, 'Jumps', t_pts);
if SAMPLE

    out_sim= dde23(exp_delay, tau, exp_history, [0,Tmax], options);            
end

cl = linspecer(1);

%% Plot
figure(1)
clf
hold on 
plot(out_sim.x, out_sim.y, 'color', cl,'LineWidth', 2);

for i = 1:(Npts_h-1)
    plot(t_pts(i:(i+1)), x_pts(i)*[1,1], 'color', cl, 'LineWidth', 2)
end

text_vert = 1;
text_horz = tau/2 - 0.1;
FS_text = 14;
text(text_horz - tau-0.1, text_vert, '$x_h(t)$', 'interpreter', 'latex', 'FontSize', FS_text)
for i = 0:floor(Tmax/tau)        
    if i==0
        LW = 2;        
    else
        LW = 1;
    end
    plot(i*[1, 1], ylim, ':k', 'LineWidth', LW)
    text_cont = ['$C^', num2str(i), '$'];
    if i<Tmax
        text(text_horz + i*tau, text_vert, text_cont, 'interpreter', 'latex', 'FontSize', FS_text)
    end
end

xlabel('time')
ylabel('x(t)')
title('Increasing Continuity', 'FontSize', 14)


%% Function
function [hh, t_pts, x_pts] = history_handle(sampler, tau, Npts)

t0 = linspace(-tau, 0, Npts);
t_pts = t0;
t_pts(2:end-1) = t_pts(2:end-1) + (1/(Npts*2))*(2*rand(1, Npts-2)-1);

% t_pts = rand(Npts-2, 1)*(-tau);
% t_pts = [-tau; sort(t_pts); 0];
x_pts = zeros(Npts, 1);
for i = 1:Npts-1
    x_pts(i) = sampler();
end
x_pts(Npts) = x_pts(Npts-1);
hh = @(x) zoh(t_pts, x_pts, x);

end