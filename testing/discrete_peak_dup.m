% Peak estimation on trajectories of a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
% the objective is to maximize is p(x)
% contains delayed states in the initial and peak measure
%
% Author: Jared Miller
%         Feb 2, 2021.

%% parameters
PLOT = 1;
PLOT_NONNEG = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
tau = 0.4;  %delay x(t - tau)
K0 = 1;     %gain in dynamics x(t)
K1 = 4;     %gain in dynamics x(t-tau)
order = 4;
% tau = 0.3;
% K0 = 0.5;
% K1 = 3;
% 
% T = 1.5;
% tau = 0.25;
% K0 = 3;
% K1 = 5;
% order = 8;      %relaxation order


% tau = 0.1;
% K0 = 1;
% K1 = 10;
% order = 8;      %relaxation order


% tau = 0;

p = @(x) x;     %objective in peak estimation


%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0);
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

mset clear
mpol('t','x0','x1');mu  = meas(t, x0, x1);  %joint occupation measure
mpol('tp', 'xp0', 'xp1');   mup = meas(tp, xp0, xp1);     %peak measure
mpol('tnn', 'xnn'); nun = meas(tnn, xnn);   %component -1
mpol('tnz', 'xnz'); nuz = meas(tnz, xnz);   %component  0
mpol('tnp', 'xnp'); nup = meas(tnp, xnp);   %component  1

%absolute continuity for peak estimation
mpol('tnpc', 'xn1c'); nunc = meas(tnpc, xn1c);   %component 1 [0, kappa T] complement

%support constraints
% tp == T
supp_con = [tp*(T-tp) >= 0; 
            tnz*((T - tau) - tnz) >= 0; 
            (tnp  - (T - tau)) * (T - tnp)  >= 0;
            tnpc * (T - tnpc) >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xp0^2 <= xh0^2;
            xp1^2 <= xh0^2;
            xnz^2 <= xh0^2;
            xnp^2 <= xh0^2;
            xn1c^2 <= xh0^2;];
%                     xnn == xh0; % redundant
%             (tnn + tau) * tnn >= 0; %redundant
        
%% monomials
yp = mmon([tp, xp0, xp1], d);
        
%% reference variables and measures
dv = genPowGlopti(3,d);
y0 = prod([0, xh0, xh0].^dv, 2); %initial measure
dv0 = genPowGlopti(2,d);
yh = history_mom(tau, xh0, dv0);


%% Affine relation
%Liouville Equation
v = mmon([t; x0; x1], d);
f = -K0*x0 -K1*x1; %dynamics
Ay = mom(diff(v, t) + diff(v, x0)*f); 
Liou = Ay + y0 - mom(yp);
% Liou = -(Ay + y0) + mom(yp);

%(t, x0) marginal
%test functions on components (t)^alpha x^beta
v0  = mmon([t; x0], d);

phi0 = mom(v0) - mom(mmon([tnz; xnz], d)) - mom(mmon([tnp; xnp], d));

%(t, x1) marginal (with shifting)
%test functions on components (t+tau)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, T-tau] component
mnz = mmon([tnz; xnz], d);
mnz_shift = subs(mnz, tnz, tnz + tau);  

%[-tau, 0] component]
mnn = mmon([tnn; xnn], d);
mnn_shift = subs(mnn, tnn, tnn + tau);

%absolute continuity
yc = mom(mmon([tnpc; xn1c], d));

phi1 = mom(v1) + yc - mom(mnz_shift) - mom(mnn_shift);



%% moment constraints
mom_con = [-Liou == 0;        %Liouville
           -phi0 == 0;        %x(t)
           -phi1 == 0;        %x(t - tau)
           mom(mnn) == yh];  %history given

%% Solve problem
objective = max(mom(p(xp0))); 
% objective = max(mass(nuz)); 
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
P = msdp(objective, ...
    mom_con, supp_con);

% solve LMIP moment problem
[status,obj_rec, m,dual_rec]= msol(P);
obj_rec


%peak extraction/estimation
Mp = double(mmat(mup));
Mp_1 = Mp(1:4, 1:4);
rp = rank(Mp_1, 1e-3);

%% Dual Functions

dual_rec_1 = dual_rec{1};
n_monom = length(Ay);
n_monom0 = length(phi0);
coeff_v    = dual_rec_1(1:n_monom);
coeff_phi0 = dual_rec_1(n_monom + (1:n_monom0));
coeff_phi1 = dual_rec_1(n_monom + n_monom0 + (1:n_monom0));

v_rec = coeff_v'*v;
Lv_rec = diff(v_rec, t) + f * diff(v_rec, x0);
v_f    = @(te, x0e, x1e) eval(coeff_v'*v, [t; x0; x1], [te; x0e; x1e]);
Lv_f    = @(te, x0e, x1e) eval(Lv_rec, [t; x0; x1], [te; x0e; x1e]);

phi0_f = @(te, xe) eval(coeff_phi0'*v1, [t; x1], [te; xe]);
phi1_f = @(te, xe) eval(coeff_phi1'*v1, [t; x1], [te; xe]);
%% Analyze trajectory

%trajectory interpolator
ci = spline(sol.x, sol.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-tau ci.breaks];
ci.coefs = [zeros(1, size(ci.coefs, 2)-1) xh0; ci.coefs];

Nt = 1000;
t_traj = linspace(0, T, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj - tau);
[peak_traj,  i_traj] = max(p(x0_traj));

% nonneg_T = v_f(T, x0_traj(end));

nonneg_cost = v_f(t_traj, x0_traj, x1_traj) - p(x0_traj);

nonneg_flow = -(Lv_f(t_traj, x0_traj, x1_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));

tnu0_traj = linspace(0, T-tau, Nt);
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj) + phi1_f(tnu0_traj + tau, xnu0_traj);

tnu1_traj = linspace(T-tau, T, floor(Nt/4));
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj);

nonneg_1c = -phi1_f(t_traj, x0_traj);

if PLOT
    hold on
    scatter(t_traj(i_traj), x0_traj(i_traj), 300, '*r', 'DisplayName', 'true peak')
    plot(xlim, obj_rec*[1, 1], '-.r', 'LineWidth', 3, 'DisplayName', 'recovered peak')
    hold off
    legend('location', 'southeast')
end
    
 if  PLOT_NONNEG    
    figure(2)
    clf
%     tiledlayout(3, 1)
    subplot(5,1,1)
    plot(t_traj, nonneg_flow)
    xlim([0, T])
    title('$\mu: \quad -(L_f v(t, x_0, x_1) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
%     nexttile
    subplot(5,1,2)
    plot(tnu0_traj, nonneg_0)
    xlim([0, T - tau])
    title('$\nu_0: \quad \phi_0(t, x) + \phi_1(t + \tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
        
    subplot(5,1,3)
    plot(t_traj, nonneg_1c)
    xlim([0, T])
    title('$\hat{\nu}_1: \quad -\phi_1(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    
    subplot(5,1,4)
    plot(tnu1_traj, nonneg_1)
    xlim([T-tau, T])
    title('$\nu_1: \quad \phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    subplot(5, 1, 5 )
    plot(t_traj, nonneg_cost)
    xlim([0, T])
    title('$\mu_p: \quad v(t, x0, x_1) - p(x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
 

end


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