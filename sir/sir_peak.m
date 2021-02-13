%An SIR model of an epidemic for varying time delay
%peak espimation: bound maximum infection rate
%
%Only SI, R can be recovered later by integrating I
%Jared Miller, 2/12/2021

%% model parameters
beta = 0.4;
gamma = 0.1;

Tmax = 30;  %max time of simulation (days)
tau = 9;    %incubation period (days)

Trange = linspace(0, Tmax, 200);

I0 = 0.2;   %initial infection rate, conspant history
xh0 = [1-I0; I0];
x00 = xh0;

p = @(s, i) i; %objective to maximize

SOLVE       = 1;
PLOT_TRAJ   = 1;
PLOT_NONNEG = 1;

order = 3;

%% sample trajectory
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'Maxstep', 0.1);

sir_delay = @(t,y,Z) Tmax * [-beta*y(1)*y(2);
            beta*Z(1)*Z(2) - gamma*(y(2))];

sir_history = @(t) xh0;

sol = dde23(sir_delay, tau/Tmax, sir_history, Trange/Tmax, options);

%empirical moments


if PLOT_TRAJ
    figure(1)
    clf
    plot([-tau, Tmax * sol.x], [xh0(2), sol.y(2, :)], 'DisplayName', ['Delay=', num2str(tau)])
    
    title('Infection Rate of Epidemic', 'FontSize', 16)
    xlabel('time (days)')
    ylabel('infection rate')
end

%% Set up variables and measures
if SOLVE
    mset clear
    mpol('t','s0','s1', 'i0', 'i1');mu  = meas(t, s0, s1, i0, i1);  %joint occupation measure
    mpol('tp', 'sp', 'ip');   mup = meas(tp, sp, ip);     %final measure
    mpol('tnn', 'snn', 'inn'); nun = meas(tnn, snn, inn);   %component -1 [-tau, 0]
    mpol('tnz', 'snz', 'inz'); nuz = meas(tnz, snz, inz);   %component  0 [0, T-tau]
    mpol('tnp', 'snp', 'inp'); nup = meas(tnp, snp, inp);   %component  1 [T-tau, T]
    %absolute continuity for peak espimation
    mpol('tnpc', 'sn1c', 'in1c'); nunc = meas(tnpc, sn1c, in1c);   %component 1 complement [0, T]

    %support constraints
    taus = tau/Tmax;
    %dynamics scaled to t in [0, 1].
    supp_con = [t * (1-t) >= 0;
                tp * (1 - tp) >= 0;
                tnz * ( 1 - taus - tnz) >= 0;
                (tnp - (1-taus)) * (1 - tnp) >= 0;
                tnpc * (1-tnpc) >= 0;
                tri_con(s0, i0);
                tri_con(s1, i1);
                tri_con(snz, inz);
                tri_con(snp, inp);
                tri_con(sn1c, in1c);
    ];


    %% reference variables and measures
    d = 2*order;
    dv = genPowGlopti(3,d);


    y0 = prod([0; x00]'.^dv, 2); %initial measure
    yh = history_mom(taus, xh0, dv);

    %% Affine relations
    %Liouville Equation
    yp = mmon([tp, sp, ip], d);
    v0  = mmon([t; s0; i0], d);

    f = Tmax * [-beta*s0*i0; beta*s1*i1 - gamma*i0]; %dynamics
    Ay = mom(diff(v0, t) + diff(v0, [s0; i0])*f); 
    Liou = Ay + y0 - mom(yp);

    %(t, x0) marginal
    %test functions on components (t)^alpha x^beta
    phi0 = mom(v0) - mom(mmon([tnz; snz; inz], d)) - mom(mmon([tnp; snp; inp], d));

    %(t, x1) marginal (with shifting)
    %test functions on components (t+tau)^alpha x^beta
    % joint occupation measure
    v1  = mmon([t; s1; i1], d);

    %[0, T-tau] component
    mnz = mmon([tnz; snz; inz], d);
    mnz_shift = subs(mnz, tnz, tnz + taus);  

    %[-tau, 0] component
    mnn = mmon([tnn; snn; inn], d);
    mnn_shift = subs(mnn, tnn, tnn + taus);

    %absolute continuity
    yc = mom(mmon([tnpc; sn1c; in1c], d));
    
    phi1 = mom(v1) + yc - mom(mnz_shift) - mom(mnn_shift);

    %% moment constraints
    mom_con = [-Liou == 0;        %Liouville
               -phi0 == 0;        %x(t)
               -phi1 == 0;        %x(t - tau)
               mom(mnn) == yh];   %history given

    %% Solve problem
    objective = max(mom(p(sp, ip)));
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
end

%% Recover dual functions

dual_rec_1 = dual_rec{1};
n_monom = length(Ay);
coeff_v    = dual_rec_1(1:n_monom);
coeff_phi0 = dual_rec_1(n_monom + (1:n_monom));
coeff_phi1 = dual_rec_1(2*n_monom + (1:n_monom));

v_rec = coeff_v'*v0;
Lv_rec = diff(v_rec, t) + diff(v_rec, [s0; i0])*f;
v_f    = @(te, xe) eval(v_rec, [t; s0; i0], [te; xe]);
Lv_f    = @(te, x0e, x1e) eval(Lv_rec, [t; s0; i0; s1; i1], [te; x0e; x1e]);

phi0_f = @(te, xe) eval(coeff_phi0'*v0, [t; s0; i0], [te; xe]);
phi1_f = @(te, xe) eval(coeff_phi1'*v0, [t; s0; i0], [te; xe]);

%% Analyze trajectory
%trajectory interpolator
ci = spline(sol.x, sol.y);
ci.pieces = ci.pieces + 1;
ci.breaks = [-taus ci.breaks];
ci.coefs = [zeros(2, size(ci.coefs, 2)-1) xh0; ci.coefs];

%interpolated trajectories
Nt = 500;
t_traj = linspace(0, 1, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj - taus);
[peak_traj,  i_traj] = max(p(x0_traj(1, :), x0_traj(2, :)));

%nonnegative functions
nonneg_cost = v_f(t_traj, x0_traj) - p(x0_traj(1, :), x0_traj(2, :));

nonneg_flow = -(Lv_f(t_traj, x0_traj, x1_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));

tnu0_traj = linspace(0, 1-taus, Nt);
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj) + phi1_f(tnu0_traj + taus, xnu0_traj);

tnu1_traj = linspace(1-taus, 1, floor(Nt/4));
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj);

nonneg_1c = -phi1_f(t_traj, x0_traj);

if PLOT_TRAJ
    hold on
    scatter(t_traj(i_traj)*Tmax, x0_traj(2, i_traj), 300, '*r', 'DisplayName', 'true peak')
    plot(xlim, obj_rec*[1, 1], '-.r', 'LineWidth', 3, 'DisplayName', 'bound peak')
    if rp == 1
        scatter(Mp_1(1, 2)*Tmax, Mp_1(1, 4), 300, 'or', 'DisplayName', 'recovered peak')
    end
    hold off
    legend('location', 'northwest')
end

if PLOT_NONNEG
       figure(2)
 subplot(5,1,1)
    plot(Tmax*t_traj, nonneg_flow)
    xlim([0, Tmax])
    title('$\mu: \quad -(L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
%     nexttile
    subplot(5,1,2)
    plot(Tmax*tnu0_traj, nonneg_0)
    xlim([0, Tmax - tau])
    title('$\nu_0: \quad \phi_0(t, x) + \phi_1(t + \tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
        
    subplot(5,1,3)
    plot(Tmax*t_traj, nonneg_1c)
    xlim([0, Tmax])
    title('$\hat{\nu}_1: \quad -\phi_1(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    
    subplot(5,1,4)
    plot(Tmax*tnu1_traj, nonneg_1)
    xlim([Tmax-tau, Tmax])
    title('$\nu_1: \quad \phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    subplot(5, 1, 5 )
    plot(Tmax*t_traj, nonneg_cost)
    xlim([0, Tmax])
    title('$\mu_p: \quad v(t, x) - p(x)$', 'interpreter', 'latex', 'fontsize', 14)
    hold on
    plot(xlim, [0, 0], ':k')
    hold off
    
    figure(4)
    plot(Tmax * t_traj, v_f(t_traj, x0_traj))
    xlabel('time')
    title('Value function over time', 'fontsize', 16)
    ylabel('$v(t, x)$', 'interpreter', 'latex')
end

function m = history_mom(tau, x0, dv)
    %moments of (t: lebesgue of [-tau, 0]) times (x: delta at x0)
    alpha = dv(:, 1);
    beta  = dv(:, 2:end);
    
    tmom = -(-tau).^(alpha+1)./(alpha+1);
    xmom = prod(x0'.^beta, 2);
    
    m = xmom .* tmom;
end

function em = monom_int(t, x, dv)
    %numerically integrate t^alpha x(t)^beta in span [0, T]
    %em: empirical moments
    %there is probably a more efficient way to implement this
    n_monom = size(dv, 1);
    em = zeros(n_monom, 1);
    for i = 1:n_monom        
        v_curr = (t.^dv(i, 1)) .* prod(x.^(dv(i, 2:end)'), 1);
%         em(i) = trapz(t, v_curr);
        em(i) = simps(t, v_curr);
    end
end

function c = tri_con(s, i)  
    %support conspraint contained in a unip simplex
    c = [s >= 0; i>= 0; s+i <= 1];
end