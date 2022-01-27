% Peak estimation on trajectories of a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
% the objective is to maximize is p(x)
%
% Author: Jared Miller
%         Feb 2, 2021.

%% parameters
PLOT = 1;
PLOT_NONNEG = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
% tau = 0.4;  %delay x(t - tau)
% K0 = 1;     %gain in dynamics x(t)
% K1 = 4;     %gain in dynamics x(t-tau)
order = 4;
tau = 0.3;
K0 = 0.5;
K1 = 4;

p = @(x) x;     %objective in peak estimation


%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/T);
%ddesd(dynamics, delay, history, time range, options)
% sol = ddesd(@(t,y,z) -K*z, @(t,y) t-tau,@(t) xh0,[0,T], options);
% sol = dde23(@(t,y,z) -K*z, [tau],@(t) xh0,[0,T], options);
if tau == 0
    sol = ode45(@(t,y) -K0*y, [0, T], xh0, options);
else
    sol = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0,[0,T], options);
end

if PLOT
    figure(1)
    clf
    hold on
    plot([-tau sol.x], [-1 sol.y], 'DisplayName', 'trajectory')
    plot([-tau, T], [0, 0], ':k')
    xlim([-tau, T])
    hold off
    xlabel('t')
    ylabel('x(t)')
    title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(t-', num2str(tau),')$'], 'interpreter', 'latex', 'fontsize', 16)
end

%% Set up variables and measures
d = 2*order;


    

function m = history_mom(tau, x0, dv)
    %moments of (t: lebesgue of [-tau, 0]) times (x: delta at x0)
    alpha = dv(:, 1);
    beta  = dv(:, 2);
    
    tmom = -(-tau).^(alpha+1)./(alpha+1);
    xmom = x0.^beta;
    
    m = xmom .* tmom;
end

function em = monom_int(t, x, dv)
    %numerically integrate t^alpha x(t)^beta in span [0, T]
    %em: empirical moments
    %there is probably a more efficient way to implement this
    n_monom = size(dv, 1);
    em = zeros(n_monom, 1);
    for i = 1:n_monom        
        v_curr = (t.^dv(i, 1)) .* (x.^dv(i, 2));
%         em(i) = trapz(t, v_curr);
        em(i) = simps(t, v_curr);
    end
end