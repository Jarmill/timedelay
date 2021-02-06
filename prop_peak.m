% peak estimation of a trajectory solving a delay differential equation
% x'(t) = -K0 x(t) - K1(kappa t) with some initial condition (no history)
% objective here is p(x) = x;
% 
% Author: Jared Miller
%         Feb 3, 2021.

%% parameters
PLOT = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
kappa = 0.5;  %delay x(t - tau)
K0 = 1;     %gain in dynamics x(t)
K1 = 3;     %gain in dynamics x(t-tau)

p = @(x) x;

order = 4;      %relaxation order

%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9);
%ddesd(dynamics, delay, history, time range, options)
% sol = ddesd(@(t,y,z) -K*z, @(t,y) t-tau,@(t) xh0,[0,T], options);
% sol = dde23(@(t,y,z) -K*z, [tau],@(t) xh0,[0,T], options);
if kappa == 0
    sol = ode45(@(t, y) -K0*y, [0, T], xh0, options);
else
    sol = ddesd(@(t, y, z) -K0*y - K1*z,@(t, y) kappa*t,@(t) xh0,[0,T], options);
end

if PLOT
    figure(1)
    clf
    hold on
    plot(sol.x, sol.y)
    plot([0 T], [0, 0], ':k')
    xlim([0 T])
    hold off
    xlabel('t')
    ylabel('x(t)')
    title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(', num2str(kappa),'t)$'], 'interpreter', 'latex', 'fontsize', 16)
end

%% Set up variables and measures
d = 2*order;

mset clear
mpol('t','x0','x1');mu  = meas(t, x0, x1);  %joint occupation measure
mpol('tp', 'xp');   mup = meas(tp, xp);     %peak measure
mpol('tn0', 'xn0'); nun = meas(tn0, xn0);   %component 0 [kappa T, T]
mpol('tn1', 'xn1'); nuz = meas(tn1, xn1);   %component 1 [0, kappa T]

%absolute continuity for peak estimation
mpol('tn1c', 'xn1c'); nunc = meas(tn1c, xn1c);   %component 1 [0, kappa T] complement
%support constraints
supp_con = [tp*(T - tp) >= 0; 
            tn1*(kappa*T - tn1) >= 0; 
            tn1c*(T - tn1c) >= 0; 
            (tn0 - kappa*T) * (T - tn0) >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xp^2 <= xh0^2;
            xn1^2 <= xh0^2;            
            xn1c^2 <= xh0^2;
            xn0^2 <= xh0^2];
 
%% monomials
yp = mmon([tp, xp], d);
        
%% reference variables and measures
dv = genPowGlopti(2,d);
y0 = prod([0, xh0].^dv, 2); %initial measure


%% Affine relation
%Liouville Equation
v0  = mmon([t; x0], d);
f = -K0*x0 -K1*x1; %dynamics
Ay = mom(diff(v0, t) + diff(v0, x0)*f); 
Liou = Ay + y0 - mom(yp);

%be very careful about signs
% Liou = -(Ay + y0) + mom(yp);

%(t, x0) marginal
%test functions on components (t)^alpha x^beta
phi0 = mom(v0) - mom(mmon([tn1; xn1], d)) - mom(mmon([tn0; xn0], d));

%(t, x1) marginal (with shifting)
%test functions on components (1/kappa) (t/kappa)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, kappa T] component scaled up
mn1 = mmon([tn1; xn1], d);
mn1_scale= (1/kappa)* subs(mn1, tn1, tn1/kappa);  

% phi1 = mom(v1) - mom(mn1_scale);

%includes absolute continuity
phi1 = mom(v1) + mom(mmon([tn1c; xn1c], d)) -  mom(mn1_scale);

%% moment constraints
mom_con = [-Liou == 0;        %Liouville
           -phi0 == 0;        %x(t)
           -phi1 == 0];       %x(kappa t)           

%% Solve problem
objective = max(mom(p(xp))); 

mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
P = msdp(objective, ...
    mom_con, supp_con);

% solve LMIP moment problem
[status,obj_rec, m,dual_rec]= msol(P);
obj_rec

%% Dual Functions

dual_rec_1 = dual_rec{1};
n_monom = length(Ay);
coeff_v    = dual_rec_1(1:n_monom);
coeff_phi0 = dual_rec_1(n_monom + (1:n_monom));
coeff_phi1 = dual_rec_1(2*n_monom + (1:n_monom));

v_f    = @(te, xe) eval(coeff_v'*v0, [t; x0], [te; xe]);
phi0_f = @(te, xe) eval(coeff_phi0'*v0, [t; x0], [te; xe]);
phi1_f = @(te, xe) eval(coeff_phi1'*v0, [t; x0], [te; xe]);

%% Analyze trajectory

%trajectory interpolator
ci = spline(sol.x, sol.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-1 ci.breaks];
ci.coefs = [zeros(1, size(ci.coefs, 2)-1) xh0; ci.coefs];


Nt = 400;
t_traj = linspace(0, T, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj*kappa);

%peak extraction/estimation
[peak_traj,  i_traj] = max(p(x0_traj));
Mp = double(mmat(mup));
Mp_1 = Mp(1:3, 1:3);
rp = rank(Mp_1);

%nonnegative functions (hopefully)
%mup
nonneg_cost = v_f(t_traj, x0_traj) - p(x0_traj);

%mu
nonneg_flow = (v_f(t_traj, x0_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));

% nonneg_flow = -(v_f(t_traj, x0_traj)) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj);
%nu0
tnu0_traj = linspace(kappa*T, T, floor(Nt/4));
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj);

%nu1
tnu1_traj = linspace(0, kappa*T, Nt);
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj) + phi1_f(tnu1_traj/kappa, xnu1_traj)/kappa;

%nu1c
nonneg_1c = -phi1_f(tnu1_traj/kappa, xnu1_traj)/kappa;

if PLOT
    figure(2)
    clf
%     tiledlayout(3, 1)
%     nexttile
    subplot(5,1,1)
    plot(t_traj, nonneg_flow)
    title('$\mu: \quad -(L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
%     nexttile
    subplot(5,1,2)
    plot(tnu1_traj, nonneg_1)
    title('$\nu_1: \quad \phi_0(t, x) + \kappa^{-1} \phi_1(\kappa^{-1} t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
        
    subplot(5,1,3)
    plot(tnu1_traj, nonneg_1c)
    title('$\hat{\nu}_1: \quad\phi_1(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    
    subplot(5,1,4)
    plot(tnu0_traj, nonneg_0)
    title('$\nu_0: \quad \phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    subplot(5, 1, 5 )
    plot(t_traj, nonneg_cost)
    title('$\mu_p: \quad v(t, x) - p(x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
%     figure(3)
%     semilogy((abs(m_traj - m_mom)), 'o')
%     title('Error in moment estimation')
%     xlabel('moment index')
%     ylabel('$\mid m_{\alpha \beta} - \hat{m}_{\alpha \beta} \mid$', 'interpreter', 'latex', 'fontsize', 14) 
%     grid on
end