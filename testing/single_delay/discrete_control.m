% Optimal Control of a trajectory of a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
%
% Author: Jared Miller
%         May 27, 2021.


%% parameters
PLOT = 0;
PLOT_NONNEG = 0;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]

tau = 0.25;
K0 = 3;
K1 = 5;
order = 3;
% order = 6;
% order = 5;
% order=7;


%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/T);
sol = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0,[0,T], options);

if PLOT
    figure(1)
    clf
    
    subplot(2,1,1)
    hold on
    plot([-tau sol.x], [xh0 sol.y], 'DisplayName', 'Open Loop')
%     plot([-tau, T], [0, 0], ':k')
    plot([-tau, T], [0, 0], ':k', 'HandleVisibility', 'off')
    xlim([-tau, T])
    hold off
    xlabel('t')
    ylabel('x(t)')
    title(['$\dot{x}(t) = -', num2str(K0), 'x(t) -', num2str(K1), 'x(t-', num2str(tau),') + u(t)$'], 'interpreter', 'latex', 'fontsize', 16)
end

%% Set up variables and measures
d = 2*order;

mset clear
mpol('t','x0','x1', 'u'); mu  = meas(t, x0, x1, u);  %joint occupation measure (controlled)
mpol('tp', 'xp');   mup = meas(tp, xp);     %peak measure
mpol('tnz', 'xnz'); nuz = meas(tnz, xnz);   %component  0
mpol('tnp', 'xnp'); nup = meas(tnp, xnp);   %component  1

%support constraints
% tp == T
% umax = 2;
umax = 1;
supp_con = [tp == T; 
            tnz*((T - tau) - tnz) >= 0; 
            (tnp  - (T - tau)) * (T - tnp)  >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xp^2 <= xh0^2;
            xnz^2 <= xh0^2;
            xnp^2 <= xh0^2;
            u^2 <= umax^2];       
%% monomials
yp = mmon([tp, xp], d);
        
%% reference variables and measures
dv = genPowGlopti(2,d);
y0 = prod([0, xh0].^dv, 2); %initial measure
yh = history_mom(tau, xh0, dv);
yh_shift = history_mom_shift(tau, xh0, dv);
%% Affine relation
%Liouville Equation
v0  = mmon([t; x0], d);
f0 = -K0*x0 -K1*x1 ;
f1 = 1;
f = f0 + f1*u; %dynamics
Ay = mom(diff(v0, t) + diff(v0, x0)*f); 
Liou = Ay + y0 - mom(yp);
% Liou = -(Ay + y0) + mom(yp);

%(t, x0) marginal
%test functions on components (t)^alpha x^beta
phi0 = mom(v0) - mom(mmon([tnz; xnz], d)) - mom(mmon([tnp; xnp], d));

%(t, x1) marginal (with shifting)
%test functions on components (t+tau)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, T-tau] component
mnz = mmon([tnz; xnz], d);
mnz_shift = subs(mnz, tnz, tnz + tau);  

%[-tau, 0] component
phi1 = mom(v1) - mom(mnz_shift) - yh_shift;

%costs in control
% R = 0.1;
R=0.01;
J =0.5*( x0^2 + R*u^2);
JT = 0;
%% moment constraints
mom_con = [Liou == 0;        %Liouville
           phi0 == 0;        %Consistency x(t)
           phi1 == 0];       %Consistency x(t - tau)
%% Solve problem
objective = min(mom(J)+mom(JT)); 
% objective = max(mass(nuz)); 
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
P = msdp(objective, ...
    mom_con, supp_con);

% solve LMI  moment problem
[status,obj_rec, m,dual_rec]= msol(P);
obj_rec;

%% Dual Functions

dual_rec_1 = dual_rec{1};
n_monom = length(Ay);
coeff_v    = dual_rec_1(1:n_monom);
coeff_phi0 = dual_rec_1(n_monom + (1:n_monom));
coeff_phi1 = dual_rec_1(2*n_monom + (1:n_monom));

phi1_int = coeff_phi1'*yh_shift;

v_rec = coeff_v'*v0; 
Lu_rec=diff(v_rec,x0)*f1;
Lv_rec = diff(v_rec, t) + f * diff(v_rec, x0);
v_f    = @(te, xe) eval(v_rec, [t; x0], [te; xe]);
Lv_f    = @(te, x0e, x1e, ue) eval(Lv_rec, [t; x0; x1; u], [te; x0e; x1e; ue]);
phi0_f = @(te, xe) eval(coeff_phi0'*v0, [t; x0], [te; xe]);
phi1_f = @(te, xe) eval(coeff_phi1'*v0, [t; x0], [te; xe]);
u_f =  @(te, x0e) eval(Lu_rec, [t; x0], [te; x0e])/(-R); 
u_clamp = @(te, x0e) min(umax, max(-umax, u_f(te, x0e))); 

J_f = @(x0e, ue) eval(J, [x0;u], [x0e;ue]);
JT_f = @(x0e) eval(JT, [x0], [x0e]);

obj_rec_dual = v_f(0, xh0) - phi1_int;

% fprintf('Primal: %0.5f, Dual: %0.5f\n', obj_rec, obj_rec_dual)
%% run closed loop trajectory

f_closed = @(t,y,z) -K0*y-K1*z + u_clamp(t,y);
sol_closed = dde23(f_closed, [tau],@(t) xh0,[0,T], options);

if PLOT
    hold on
    plot([-tau sol_closed.x], [xh0 sol_closed.y], 'DisplayName', 'Closed Loop')  
  legend('location', 'northwest')
end
     

%% Analyze trajectory

%trajectory interpolator
ci = spline(sol_closed.x, sol_closed.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-tau ci.breaks];
ci.coefs = [zeros(1, size(ci.coefs, 2)-1) xh0; ci.coefs];

Nt = 300;
t_traj = linspace(0, T, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj - tau);
u_traj = u_clamp(t_traj,x0_traj);

J_traj = J_f(x0_traj, u_traj);
cost_traj = simps(t_traj, J_traj) + JT_f(x0_traj(end));

fprintf('Primal:\t%0.5f\nDual:\t%0.5f\nTraj:\t%0.5f\n', obj_rec, obj_rec_dual, cost_traj)
% [peak_traj,  i_traj] = max(p(x0_traj));



nonneg_T = v_f(T, x0_traj(end)) + JT_f(x0_traj(end));
v_traj = v_f(t_traj, x0_traj);

phi1_traj = phi1_f(t_traj+tau, x0_traj);
phi1_delay_traj = phi1_f(t_traj, x1_traj);


phi1accum_traj = cumsimps(t_traj, phi1_traj-phi1_delay_traj) + phi1_int;

nonneg_flow = -Lv_f(t_traj, x0_traj, x1_traj, u_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj)...
    +J_f(t_traj, x0_traj);

tnu0_traj = linspace(0, T-tau, Nt);
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = -phi0_f(tnu0_traj, xnu0_traj) - phi1_f(tnu0_traj + tau, xnu0_traj);

tnu1_traj = linspace(T-tau, T, floor(Nt/4));
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = -phi0_f(tnu1_traj, xnu1_traj);
      

if PLOT
    subplot(2,1,2)
    hold on
    plot(t_traj, u_traj)
    xlabel('t')
    ylabel('u(t)')
    plot([-tau, T], [0, 0], ':k', 'HandleVisibility', 'off')
    xlim([-tau, T])
   title('Control $u(t)$','interpreter', 'latex', 'fontsize', 16)
    
end
%% Plot nonnegativity
 if  PLOT_NONNEG    
    figure(2)
    clf
%     tiledlayout(3, 1)
    subplot(3,1,1)
    plot(t_traj, nonneg_flow)
    xlim([0, T])
    title('$\mu: \quad J(x_0,u)-L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
%     nexttile
    subplot(3,1,2)
    plot(tnu0_traj, nonneg_0)
    xlim([0, T])
    title('$\nu_0: \quad -\phi_0(t, x) - \phi_1(t + \tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
       
    
    
    subplot(3,1,3)
    plot(tnu1_traj, nonneg_1)
    xlim([0, T])
    title('$\nu_1: \quad -\phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off


 figure(3) 
    clf
    plot(t_traj, v_traj)
    xlim([0, T])
    title('$v(t,x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
 figure(4) 
    clf
    plot(t_traj, v_traj - phi1accum_traj)
    xlim([0, T])
    title('$V(t, x(t)) = v(t,x(t)) - \int_{-\tau}^{0} \phi_1(t+s+\tau, x(t+s))ds$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
 end

%% function definitions
function m = history_mom(tau, x0, dv)
    %moments of (t: lebesgue of [-tau, 0]) times (x: delta at x0)
    alpha = dv(:, 1);
    beta  = dv(:, 2);
    
    tmom = -(-tau).^(alpha+1)./(alpha+1);
    xmom = x0.^beta;
    
    m = xmom .* tmom;
end
function m = history_mom_shift(tau, x0, dv)
    %moments of (t: lebesgue of [0, tau]) times (x: delta at x0)
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