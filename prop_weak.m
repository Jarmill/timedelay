% Estimate moments of a trajectory solving a delay differential equation
% x'(t) = -K0 x(t) - K1(kappa t) with some initial condition (no history)
%
% 
% Author: Jared Miller
%         Feb 3, 2021.

%% parameters
PLOT = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
kappa = 0.4;  %delay x(t - tau)
K0 = 1;     %gain in dynamics x(t)
K1 = 4;     %gain in dynamics x(t-tau)

% kappa = 0.2;  %delay x(t - tau)
% K0 = 1;     %gain in dynamics x(t)
% K1 = 2;     %gain in dynamics x(t-tau)
% 


% tau = 0;
order = 6;      %relaxation order

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
mpol('tT', 'xT');   muT = meas(tT, xT);     %final measure
mpol('tn0', 'xn0'); nun = meas(tn0, xn0);   %component 0 [kappa T, T]
mpol('tn1', 'xn1'); nuz = meas(tn1, xn1);   %component  1 [0, kappa T]

%support constraints
supp_con = [tT == T; 
            tn1*(kappa*T - tn1) >= 0; 
            (tn0 - kappa*T) * (T - tn0) >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xT^2 <= xh0^2;
            xn1^2 <= xh0^2;
            xn0^2 <= xh0^2];
 
%% monomials
yT = mmon([tT, xT], d);
        
%% reference variables and measures
dv = genPowGlopti(2,d);
y0 = prod([0, xh0].^dv, 2); %initial measure


%% Affine relation
%Liouville Equation
v0  = mmon([t; x0], d);
f = -K0*x0 -K1*x1; %dynamics
Ay = mom(diff(v0, t) + diff(v0, x0)*f); 
Liou = Ay + y0 - mom(yT);
% Liou = -(Ay + y0) + mom(yT);

%(t, x0) marginal
%test functions on components (t)^alpha x^beta
phi0 = mom(v0) - mom(mmon([tn1; xn1], d)) - mom(mmon([tn0; xn0], d));

%(t, x1) marginal (with shifting)
%test functions on components (1/kappa) (t/kappa)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, kappa T] component scaled up
mn1 = mmon([tn1; xn1], d);
mn1_scale= (1/kappa )* subs(mn1, tn1, tn1/kappa);  

phi1 = mom(v1) - mom(mn1_scale);

%% moment constraints
% mom_con = [Liou == 0;        %Liouville
%            phi0 == 0;        %x(t)
%            phi1 == 0];       %x(kappa t)       
mom_con = [Liou == 0;        %Liouville
           -phi0 == 0;        %x(t)
           -phi1 == 0];       %x(kappa t)       


%% Solve problem
objective = max(mass(muT)-1); 
% objective = max(mass(nuz)); 
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
P = msdp(objective, ...
    mom_con, supp_con);

% solve LMIP moment problem
[status,obj_rec, m,dual_rec]= msol(P);


%% Analyze moments

%estimation from trajectory
m_traj = monom_int(sol.x, sol.y, dv);
%(t, x0) marginal of joint occupation measure mu
m_mom  = double(mom(v0));
m_comp = [m_traj m_mom];
norm_diff = norm(m_traj-m_mom)

% disp(['Norm Gap between empirical and LMI moments = ', num2str(norm_diff)])

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

nonneg_flow = -(v_f(t_traj, x0_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));
tnu0_traj = linspace(kappa*T, T, floor(Nt/4));
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj);

tnu1_traj = linspace(0, kappa*T, Nt);
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj) + phi1_f(tnu1_traj/kappa, xnu1_traj)/kappa;


if PLOT
    figure(2)
    clf
%     tiledlayout(3, 1)
%     nexttile
    subplot(3,1,1)
    plot(t_traj, nonneg_flow)
    title('$-(L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    
%     nexttile
    subplot(3,1,2)
    plot(tnu1_traj, nonneg_1)
    title('$\phi_0(t, x) + \kappa^{-1} \phi_1(\kappa^{-1} t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    
    
%     nexttile
    subplot(3,1,3)
    plot(tnu0_traj, nonneg_0)
    title('$\phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    
    
    figure(3)
    semilogy((abs(m_traj - m_mom)), 'o')
    title('Error in moment estimation')
    xlabel('moment index')
    ylabel('$\mid m_{\alpha \beta} - \hat{m}_{\alpha \beta} \mid$', 'interpreter', 'latex', 'fontsize', 14) 
    grid on
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