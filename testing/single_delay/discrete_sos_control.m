% Optimal Control of a trajectory of a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
%
% Author: Jared Miller
%         May 28, 2021.

%% parameters
SOLVE = 1;
PLOT = 1;
PLOT_NONNEG = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]

tau = 0.25;
K0 = 3;
K1 = 5;
order =3;
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/T);

sol_open = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0,[0,T], options);

%% Set up variables and polynomials

if SOLVE
d = 2*order;

%variable declaration
t = sdpvar(1,1); %just to be safe
x = sdpvar(1,1);
x1 = sdpvar(1,1);
u  = sdpvar(1,1);

umax=1;
Xmax = 1;
X=struct('ineq', Xmax^2-x^2);
%Dynamics and costs
R = 0.01;
J = 0.5*(x^2 + R*u^2); 
JT = 0;

f0 = -K0*x-K1*x1 ;
f1 = 1;
f = T*(f0 + f1*u);

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

%maximizing objective
objective = (obj_phi + obj_v);

%% build the constraints
%terminal
vT = replace(v, t, 1);
posterm = vT+JT;
[pterm, consterm, coeffterm] = constraint_psatz(posterm, X, [x], d);

%joint occupation
Xjoint = struct('ineq', [t*(1-t); Xmax^2 - [x;x1].^2; umax^2-u^2]);
phi1joint = replace(phi1, x, x1); 

Lv = jacobian(v,t) + f*jacobian(v,x);
posjoint = J + phi0 + phi1joint - Lv;
[pjoint, consjoint, coeffjoint] = constraint_psatz(posjoint...
    , Xjoint, [t;x;x1;u], d);


%times `zero' between 0 and T-tau 
phi1_0 = replace(phi1, t, t+tau/T);
X0 = struct('ineq', [t*((T-tau)/T -t); Xmax^2-x^2]);
pos0=-(phi0+phi1_0);
[p0, cons0, coeff0] = constraint_psatz(pos0, X0, [t;x], d);

%times `one' between T-tau and T
X1 = struct('ineq', [(t-(T-tau)/T)*(1 -t); Xmax^2-x^2]);
pos1=-phi0;
[p1, cons1, coeff1] = constraint_psatz(pos1, X1, [t;x], d);


nonneg=[posjoint; pos0; pos1];
cons = [cons; consjoint:'Joint'; cons0:'Time 0'; cons1:'Time 1'; consterm:'Terminal'];

coeff = [coeff; coeffterm; coeffjoint; coeff0; coeff1];

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
[c1,m1] = coefficients(v, [t;x]);
phi1_eval = value(c1)'*m1;      
 
[cnn,mnn] = coefficients(nonneg,[t;x;x1;u]);
nn_eval = value(cnn)*mnn;

u_eval = jacobian(v_eval, x)*f1/(-R);


%% now create functions.
v_f = polyval_func(v_eval, [t;x]);
JT_f = polyval_func(JT, [x]);
J_f = polyval_func(J, [x; u]);
nonneg_f = polyval_func(nn_eval, [t;x;x1;u]);
uraw_f = polyval_func(u_eval, [t;x]);
u_f = @(te,xe) min(umax, max(-umax, uraw_f([te; xe])));

f_closed = @(t,y,z) -K0*y-K1*z + u_f(t,y);
sol_closed = dde23(f_closed, [tau],@(t) xh0,[0,T], options);
sol_closed.u = zeros(size(sol_closed.x));
sol_closed.v = zeros(size(sol_closed.x));
sol_closed.nonneg = zeros(3,size(sol_closed.x,1));
sol_closed.dJ = zeros(size(sol_closed.x));
for i = 1:length(sol_closed.x)
    sol_closed.u(i) = u_f(sol_closed.x(i), sol_closed.y(i)); 
    sol_closed.v(i) = v_f([sol_closed.x(i); sol_closed.y(i)]);
    
    sol_closed.nonneg(:, i) = nonneg_f([sol_closed.x(i); sol_closed.y(i);sol_closed.yp(i); sol_closed.u(i)]);
    if sol_closed.x(i) <= (T-tau)/T
        sol_closed.nonneg(3, i) = 0;
    else
        sol_closed.nonneg(2, i) = 0;
    end
    
    sol_closed.v(i) = v_f([sol_closed.x(i); sol_closed.y(i)]);
    
    sol_closed.dJ(i) = J_f([sol_closed.y(i); sol_closed.u(i)]);
    
    
end
sol_closed.nonneg_term = JT_f(sol_closed.y(end)) + v_f([sol_closed.x(end); sol_closed.y(end)]);
sol_closed.J = simps(sol_closed.x, sol_closed.dJ);
sol_closed.JT = JT_f(sol_closed.y(end)); 
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
    plot([-tau sol_closed.x], [xh0 sol_closed.y], 'DisplayName', 'Open Loop')
    plot([-tau, T], [0, 0], ':k')
    xlim([-tau, T])
    hold off
    xlabel('t')
    ylabel('x(t)')
    title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(t-', num2str(tau),') + u(t)$'], 'interpreter', 'latex', 'fontsize', 16)
    

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
    title('$\mu: \quad J(x_0,u)-L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1)$', 'interpreter', 'latex', 'fontsize', 14)
    
        ax2 = nexttile;
    plot(sol_closed.x, sol_closed.nonneg(2, :))
    title('$\nu_0: \quad \phi_0(t, x) + \phi_1(t + \tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    
        ax3 = nexttile;
      
    plot(sol_closed.x, sol_closed.nonneg(3, :))
     title('$\nu_1: \quad \phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    
    linkaxes([ax1,ax2,ax3],'x')

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