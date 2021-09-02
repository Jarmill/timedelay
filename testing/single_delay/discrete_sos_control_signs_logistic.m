% Optimal Control of a trajectory of a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
%
% Author: Jared Miller
%         June 9, 2021.

%This has a different sign convention and will hopefully work

%% parameters
SOLVE = 1;
PLOT = 1;
PLOT_NONNEG = 1;
T = 1;      %time horizon
xh0 = 0.5;   %constant history x(t) = xh0 for times [-tau, 0]

% T0 = 80;
T0 = 40;
% tau = 11/T0;
tau = 14/T0;
r = 0.15;
K = 1;
xstar = 0.75;
% order = 5;
order =4;
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/T);

sol_open = dde23(@(t,y,z) T0*r*y*(xstar-z/K), [tau],@(t) xh0,[0,T], options);

%% Set up variables and polynomials

if SOLVE
d = 2*order;

%variable declaration
t = sdpvar(1,1); %just to be safe
x = sdpvar(1,1);
x1 = sdpvar(1,1);
u  = sdpvar(1,1);

umax=1;
Xmax = 1.15;
X=struct('ineq', (Xmax-x)*x);
%Dynamics and costs
R = 0.1;
% R = 1;

J = 0.5*((x-1)^2 + R*u^2); 
JT = 0;

f0 = r*x*(xstar-x1/K);
f1 = 1;
% f1 = 0.5;
f = T0*(f0 + f1*u);

%auxiliary polynomials
[v, cv] = polynomial([t;x],d);
[phi0, c0] = polynomial([t;x],d);
[phi1, c1] = polynomial([t;x],d);


cons = [];
coeff = [cv; c0; c1];


%% get the objective
obj_v = replace(v, [t;x], [0; xh0]);

dv = monpowers(length([t;x]), d); 
yh = history_mom(tau/T, xh0, dv);


obj_phi = yh'*c1;

%% build the constraints

vT = replace(v, t, 1);
phi1joint = replace(phi1, x, x1); 

Lv = jacobian(v,t) + f*jacobian(v,x);

phi1_0 = replace(phi1, t, t+tau/T);



% terms
objective = -(obj_phi + obj_v);
%terminal
posterm = -vT+JT;
%joint occupation
posjoint = J-phi0 - phi1joint + Lv;
%times `zero' between 0 and T-tau 
pos0= (phi0+phi1_0);
%times `one' between T-tau and T
pos1= phi0;

%support sets
Xjoint = struct('ineq', [t*(1-t); (Xmax-[x;x1]).*[x; x1]; umax^2-u^2]);

X1 = struct('ineq', [(t-(T-tau)/T)*(1 -t); (Xmax-x)*x]);
X0 = struct('ineq', [t*((T-tau)/T -t); (Xmax-x)*x]);
%psatz expressions
[pterm, consterm, coeffterm] = constraint_psatz(posterm, X, [x], d);
[p1, cons1, coeff1] = constraint_psatz(pos1, X1, [t;x], d);
[p0, cons0, coeff0] = constraint_psatz(pos0, X0, [t;x], d);
[pjoint, consjoint, coeffjoint] = constraint_psatz(posjoint...
    , Xjoint, [t;x;x1;u], d);

nonneg=[posjoint; pos0; pos1];
cons = [cons; consjoint:'Joint'; cons0:'Time 0'; cons1:'Time 1'; consterm:'Terminal'];

coeff = [coeff; coeffjoint; coeff0; coeff1; coeffterm];

%% Solve problemv 

opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

% cons = nonneg;
% objective =0;

%I think this should be a minimizing objective, and am not sure why the
%+objective term seems to work (where the cost is negative). There is
%likely some interaction with yalmip's automatic dualization.
[sol, monom, Gram, residual] = solvesos(cons, objective, opts, coeff);
cost_val = value(objective);

%% recover the solved polynomials
[cv,mv] = coefficients(v, [t;x]);
v_eval = value(cv)'*mv;      
[c0,m0] = coefficients(phi0, [t;x]);
phi0_eval = value(c0)'*m0;      
[c1,m1] = coefficients(phi1, [t;x]);
phi1_eval = value(c1)'*m1;      
 
[cnn,mnn] = coefficients(nonneg,[t;x;x1;u]);
nn_eval = value(cnn)*mnn;

u_eval = jacobian(v_eval, x)*f1/(-R);

%% now create functions.
v_f = polyval_func(v_eval, [t;x]);
phi1_f = polyval_func(phi1_eval, [t;x]);
JT_f = polyval_func(JT, [x]);
J_f = polyval_func(J, [x; u]);
nonneg_f = polyval_func(nn_eval, [t;x;x1;u]);
uraw_f = polyval_func(u_eval, [t;x]);
u_f = @(te,xe) min(umax, max(-umax, uraw_f([te; xe])));

f_closed = @(t,y,z) T0*r*y*(xstar-z/K) + T0*f1*u_f(t,y);
sol_closed = dde23(f_closed, [tau],@(t) xh0,[0,T], options);

%interpolate the past trajectory
%only valid on constant history
ci = spline(sol_closed.x, sol_closed.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-tau ci.breaks];
ci.coefs = [zeros(1, size(ci.coefs, 2)-1) xh0; ci.coefs];
sol_closed.yd = ppval(ci, sol_closed.x - tau);

%analyze trajectory
sol_closed.u = zeros(size(sol_closed.x));
sol_closed.v = zeros(size(sol_closed.x));
sol_closed.nonneg = zeros(3,size(sol_closed.x,1));
sol_closed.dJ = zeros(size(sol_closed.x));
sol_closed.phi1 = zeros(size(sol_closed.x));
sol_closed.phi1_delay = zeros(size(sol_closed.x));
for i = 1:length(sol_closed.x)
    sol_closed.u(i) = u_f(sol_closed.x(i), sol_closed.y(i)); 
    sol_closed.v(i) = v_f([sol_closed.x(i); sol_closed.y(i)]);
    
    sol_closed.nonneg(:, i) = nonneg_f([sol_closed.x(i); sol_closed.y(i);sol_closed.yd(i); sol_closed.u(i)]);
    
    %component nonnegativity only valid in certain time ranges.
    if sol_closed.x(i) <= (T-tau)/T
        sol_closed.nonneg(3, i) = 0;
        sol_closed.phi1(i) = phi1_f([sol_closed.x(i)+tau/T; sol_closed.y(i)]);
    
    else
        sol_closed.nonneg(2, i) = 0;
        sol_closed.phi1(i) =0;
    
    end
    
    sol_closed.v(i) = v_f([sol_closed.x(i); sol_closed.y(i)]);
    
    sol_closed.dJ(i) = J_f([sol_closed.y(i); sol_closed.u(i)]);
    
    
    sol_closed.phi1_delay(i)= phi1_f([sol_closed.x(i); sol_closed.yd(i)]);
end
sol_closed.nonneg_term = JT_f(sol_closed.y(end)) - v_f([sol_closed.x(end); sol_closed.y(end)]);
sol_closed.J = simps(sol_closed.x, sol_closed.dJ);
sol_closed.JT = JT_f(sol_closed.y(end)); 
sol_closed.dphi1=sol_closed.phi1-sol_closed.phi1_delay;

sol_closed.phi1_accum = cumsimps(sol_closed.x, sol_closed.dphi1);


sol_open.dJ = 0.5*(sol_open.y-1).^2;
sol_open.J = simps(sol_open.x, sol_open.dJ);

sol_closed.value = sol_closed.v + sol_closed.phi1_accum + value(obj_phi);
end

%% create plots
if PLOT
    
    %%open vs closed loop trajectories
    figure(1)
    clf
    tiledlayout(2,1)
    nexttile;
    hold on
    plot([-tau sol_open.x], [xh0 sol_open.y], 'DisplayName', 'Open Loop')
    plot([-tau sol_closed.x], [xh0 sol_closed.y], 'DisplayName', 'Closed Loop')
    plot([-tau, T], [0, 0], ':k', 'HandleVisibility', 'off')
    xlim([-tau, T])
    hold off
    xlabel('t')
    ylabel('x(t)')
%     title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(t-', num2str(tau),') + u(t)$'], 'interpreter', 'latex', 'fontsize', 16)
    legend('location', 'northwest')

    nexttile;
    hold on
    plot(sol_closed.x, sol_closed.u)
    xlabel('t')
    ylabel('u(t)')
    plot([-tau, T], [0, 0], ':k')
    xlim([-tau, T])
   title('Control $u(t)$','interpreter', 'latex', 'fontsize', 16)
    
    
    figure(3)
    clf
    hold on
    plot(sol_closed.x, sol_closed.v)
    xlabel('t')
    ylabel('v(t)')
    plot([-tau, T], [0, 0], ':k')
     xlim([0, T])
    title('$v(t,x)$', 'interpreter', 'latex', 'fontsize', 14)
    
    
    figure(2)
    clf
    tiledlayout(3,1)
    ax1 = nexttile;
    plot(sol_closed.x, sol_closed.nonneg(1, :))
    title('$\mu: \quad J(x_0,u) + L_f v(t, x_0) - \phi_0(t, x_0) - \phi_1(t, x_1)$', 'interpreter', 'latex', 'fontsize', 14)
    
        ax2 = nexttile;
    plot(sol_closed.x, sol_closed.nonneg(2, :))
    title('$\nu_0: \quad \phi_0(t, x) + \phi_1(t + \tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    
        ax3 = nexttile;
      
    plot(sol_closed.x, sol_closed.nonneg(3, :))
     title('$\nu_1: \quad \phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    
    linkaxes([ax1,ax2,ax3],'x')
    
        figure(4)
    clf
    hold on
    plot(sol_closed.x, sol_closed.value)
    xlabel('t')
    title('$V(t, x(t)) = v(t,x(t)) + \int_{t}^{\min(t+\tau, T)} \phi_1(s, x(s-\tau))ds$', 'interpreter', 'latex', 'fontsize', 14)
    plot([0, T], [0, 0], ':k')
     xlim([0, T])
    ylabel('$V(t,x(t))$', 'interpreter', 'latex')
    

end
%% function definitions
function m = history_mom(tau, x0, dv)
    %moments of (t: lebesgue of [-tau, 0]) times (x: delta at x0)
    %Perform a shift in the time coordinate in order to accomplish this
    %simply, return to time in [0, tau]
    alpha = dv(:, 1);
    beta  = dv(:, 2);
    
    tmom = (tau).^(alpha+1)./(alpha+1);
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