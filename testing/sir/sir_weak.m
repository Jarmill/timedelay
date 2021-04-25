%An SIR model of an epidemic for varying time delay
%estimate moments of a single trajectory (weak solution)
%
%Only SI, R can be recovered later by integrating I
%Jared Miller, 2/12/2021

%% model parameters
beta = 0.4;
gamma = 0.1;

Tmax = 30;  %max time of simulation (days)
tau = 9;    %incubation period (days)

Trange = linspace(0, Tmax, 1000);


I0 = 0.2;   %initial infection rate, constant history

x00 = [1-I0; I0];
xh0 = x00;
% xh0 = [1; 0];


SOLVE       = 1;
PLOT_TRAJ   = 1;
PLOT_MOM    = 1;
PLOT_NONNEG = 1;

order = 3;

%% sample trajectory
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);

sir_delay = @(t,y,Z) Tmax * [-beta*y(1)*y(2);
            beta*Z(1)*Z(2) - gamma*(y(2))];

sir_history = @(t) xh0 * (t ~= 0) + x00 * (t == 0);

sol = dde23(sir_delay, tau/Tmax, sir_history, Trange/Tmax, options);

%empirical moments
d = 2*order;
dv = genPowGlopti(3,d);
m_traj = monom_int(sol.x, sol.y, dv);

if PLOT_TRAJ
    figure(1)
    clf
    plot([-tau, -1e-8, Tmax * sol.x], [xh0(2), xh0(2), sol.y(2, :)], 'DisplayName', ['Delay=', num2str(tau)])
    
    title('Infection Rate of Epidemic', 'FontSize', 16)
    xlabel('time (days)')
    ylabel('infection rate')
end

%% Set up variables and measures
if SOLVE
    mset clear
    mpol('t','s0','s1', 'i0', 'i1');mu  = meas(t, s0, s1, i0, i1);  %joint occupation measure
    mpol('tT', 'sT', 'iT');   muT = meas(tT, sT, iT);     %final measure
    mpol('tnn', 'snn', 'inn'); nun = meas(tnn, snn, inn);   %component -1 [-tau, 0]
    mpol('tnz', 'snz', 'inz'); nuz = meas(tnz, snz, inz);   %component  0 [0, T-tau]
    mpol('tnp', 'snp', 'inp'); nup = meas(tnp, snp, inp);   %component  1 [T-tau, T]

    %support constraints
    taus = tau/Tmax;
    %dynamics scaled to t in [0, 1].
    supp_con = [t * (1-t) >= 0;
                tT == 1;
                tnz * ( 1 - taus - tnz) >= 0;
                (tnp - (1-taus)) * (1 - tnp) >= 0;
                tri_con(s0, i0);
                tri_con(s1, i1);
                tri_con(snz, inz);
                tri_con(snp, inp);
    ];


    %% reference variables and measures


    y0 = prod([0; x00]'.^dv, 2); %initial measure
    yh = history_mom(taus, xh0, dv);

    %% Affine relations
    %Liouville Equation
    yT = mmon([tT, sT, iT], d);
    v0  = mmon([t; s0; i0], d);

    f = Tmax * [-beta*s0*i0; beta*s1*i1 - gamma*i0]; %dynamics
    Ay = mom(diff(v0, t) + diff(v0, [s0; i0])*f); 
    Liou = Ay + y0 - mom(yT);

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

    %[-tau, 0] component]
    mnn = mmon([tnn; snn; inn], d);
    mnn_shift = subs(mnn, tnn, tnn + taus);

    phi1 = mom(v1) - mom(mnz_shift) - mom(mnn_shift);

    %% moment constraints
    mom_con = [-Liou == 0;        %Liouville
               -phi0 == 0;        %x(t)
               -phi1 == 0;        %x(t - tau)
               mom(mnn) == yh];   %history given

    %% Solve problem
    objective = max(mass(mnn)); 
    % objective = max(mass(nuz)); 
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));
    P = msdp(objective, ...
        mom_con, supp_con);

    % solve LMIP moment problem
    [status,obj_rec, m,dual_rec]= msol(P);
end

%% Analyze moments

m_mom  = double(mom(v0));
m_comp = [m_traj m_mom];
norm_diff = norm(m_traj-m_mom)

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
Nt = 400;
t_traj = linspace(0, 1, Nt);
x0_traj = ppval(ci, t_traj);
x1_traj = ppval(ci, t_traj - taus);

%nonnegative functions
nonneg_T = v_f(1, x0_traj(:, end));
nonneg_flow = -(Lv_f(t_traj, x0_traj, x1_traj) + phi0_f(t_traj, x0_traj) + phi1_f(t_traj, x1_traj));

tnu0_traj = linspace(0, 1-taus, Nt);
xnu0_traj = ppval(ci, tnu0_traj);
nonneg_0    = phi0_f(tnu0_traj, xnu0_traj) + phi1_f(tnu0_traj + taus, xnu0_traj);

tnu1_traj = linspace(1-taus, 1, floor(Nt/4));
xnu1_traj = ppval(ci, tnu1_traj);
nonneg_1    = phi0_f(tnu1_traj, xnu1_traj);


if PLOT_MOM
    figure(3)
    semilogy((abs(m_traj - m_mom)), 'o')
    title('Error in moment estimation', 'FontSize', 16)
    xlabel('moment index')
    ylabel('$\mid m_{\alpha \beta} - \hat{m}_{\alpha \beta} \mid$', 'interpreter', 'latex', 'fontsize', 14) 
    grid on
end

if PLOT_NONNEG
       figure(2)
    clf
%     tiledlayout(3, 1)
    subplot(3,1,1)
%     nexttile
    plot(Tmax * t_traj, nonneg_flow)
    title('$-(L_f v(t, x_0) + \phi_0(t, x_0) + \phi_1(t, x_1))$', 'interpreter', 'latex', 'fontsize', 14)
    xlim([0, Tmax])
%     nexttile
    subplot(3,1,2)
    plot(Tmax * tnu0_traj, nonneg_0)
    title('$\phi_0(t, x) + \phi_1(t+\tau, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlim([0, Tmax-tau])
%     nexttile
    subplot(3,1,3)
    plot(Tmax * tnu1_traj, nonneg_1)
    title('$\phi_0(t, x)$', 'interpreter', 'latex', 'fontsize', 14)
    xlabel('time')
    xlim([Tmax-tau, Tmax])
    
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
    %support constraint contained in a unit simplex
    c = [s >= 0; i>= 0; s+i <= 1];
end