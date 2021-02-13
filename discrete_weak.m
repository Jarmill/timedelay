% Estimate moments of a trajectory solving a delay differential equation
% x'(t) = -K0 x(t) -K1 x(t - tau) where the state history is constant
%
% Author: Jared Miller
%         Feb 2, 2021.

%% parameters
PLOT = 1;
PLOT_NONNEG = 1;
T = 1;      %time horizon
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
% x00 = xh0;  %discontinuity at time 0
x00 = 0;  %discontinuity at time 0

% tau = 0.4;  %delay x(t - tau)
% K0 = 1;     %gain in dynamics x(t)
% K1 = 4;     %gain in dynamics x(t-tau)

% tau = 0.3;
% K0 = 0.5;
% K1 = 3;

% T = 1.5;
tau = 0.25;
K0 = 3;
K1 = 5;

% tau = 0.1;
% K0 = 1;
% K1 = 10;
order = 6;      %relaxation order

%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0);
%ddesd(dynamics, delay, history, time range, options)
% sol = ddesd(@(t,y,z) -K*z, @(t,y) t-tau,@(t) xh0,[0,T], options);
% sol = dde23(@(t,y,z) -K*z, [tau],@(t) xh0,[0,T], options);
if tau == 0
    sol = ode45(@(t, y) -K0*y, [0, T], xh0, options);
else
    sol = dde23(@(t,y,z) -K0*y-K1*z, [tau],@(t) xh0*(t~=0) + x00*(t==0),[0,T], options);
end

if PLOT
    figure(1)
    clf
    hold on
    plot([-tau -1e-8 sol.x], [xh0 xh0 sol.y])
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
mpol('tT', 'xT');   muT = meas(tT, xT);     %final measure
mpol('tnn', 'xnn'); nun = meas(tnn, xnn);   %component -1
mpol('tnz', 'xnz'); nuz = meas(tnz, xnz);   %component  0
mpol('tnp', 'xnp'); nup = meas(tnp, xnp);   %component  1

%support constraints
supp_con = [tT == T; 
            tnz*((T - tau) - tnz) >= 0; 
            (tnp - (T - tau)) * (T - tnp) >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xT^2 <= xh0^2;
            xnz^2 <= xh0^2;
            xnp^2 <= xh0^2];
%                     xnn == xh0; % redundant
%             (tnn + tau) * tnn >= 0; %redundant
        
%% monomials
yT = mmon([tT, xT], d);
        
%% reference variables and measures
dv = genPowGlopti(2,d);
y0 = prod([0, x00].^dv, 2); %initial measure
yh = history_mom(tau, xh0, dv);


%% Affine relation
%Liouville Equation
v0  = mmon([t; x0], d);
f = -K0*x0 -K1*x1; %dynamics
Ay = mom(diff(v0, t) + diff(v0, x0)*f); 
Liou = Ay + y0 - mom(yT);
% Liou = -(Ay + y0) + mom(yT);

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

%[-tau, 0] component]
mnn = mmon([tnn; xnn], d);
mnn_shift = subs(mnn, tnn, tnn + tau);

phi1 = mom(v1) - mom(mnz_shift) - mom(mnn_shift);



%% moment constraints
mom_con = [-Liou == 0;        %Liouville
           -phi0 == 0;        %x(t)
           -phi1 == 0;        %x(t - tau)
           mom(mnn) == yh];  %history given

%% Solve problem
objective = max(mass(mnn)); 
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

v_rec = coeff_v'*v0;
Lv_rec = diff(v_rec, t) + f * diff(v_rec, x0);
v_f    = @(te, xe) eval(v_rec, [t; x0], [te; xe]);
Lv_f    = @(te, x0e, x1e) eval(Lv_rec, [t; x0; x1], [te; x0e; x1e]);

phi0_f = @(te, xe) eval(coeff_phi0'*v0, [t; x0], [te; xe]);
phi1_f = @(te, xe) eval(coeff_phi1'*v0, [t; x0], [te; xe]);

%% Analyze trajectory

%trajectory interpolator
ci = spline(sol.x, sol.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-tau ci.breaks];
ci.coefs = [zeros(1, size(ci.coefs, 2)-1) xh0; ci.coefs];


Nt = 400;
t_traj = linspace(0, T, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj - tau);

nonneg_T = v_f(T, x0_traj(end));

nonneg_flow = -(Lv_f(t_traj, x0_traj, x1_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));

tnu0_traj = linspace(0, T-tau, Nt);
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj) + phi1_f(tnu0_traj + tau, xnu0_traj);

tnu1_traj = linspace(T-tau, T, floor(Nt/4));
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj);


if PLOT
    
    figure(3)
    semilogy((abs(m_traj - m_mom)), 'o')
    title('Error in moment estimation')
    xlabel('moment index')
    ylabel('$\mid m_{\alpha \beta} - \hat{m}_{\alpha \beta} \mid$', 'interpreter', 'latex', 'fontsize', 14) 
    grid on
end
    
 if PLOT || PLOT_NONNEG    
    figure(2)
    clf
%     tiledlayout(3, 1)
    subplot(3,1,1)
%     nexttile
    plot(t_traj, nonneg_flow)
    title('$-(L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    xlim([0, T])
%     nexttile
    subplot(3,1,2)
    plot(tnu0_traj, nonneg_0)
    title('$\phi_0(t, x) + \phi_1(t+\tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlim([0, T-tau])
%     nexttile
    subplot(3,1,3)
    plot(tnu1_traj, nonneg_1)
    title('$\phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    xlim([T-tau, T])
    
    figure(4)
    plot(t_traj, v_f(t_traj, x0_traj))
    xlabel('time')
    title('Value function over time', 'fontsize', 16)
    ylabel('$v(t, x)$', 'interpreter', 'latex')

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