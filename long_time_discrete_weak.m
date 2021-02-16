% Discrete system with long discrete delays
% x[k+1] = -k0 x[k] - ktau x[k-tau]
% weak solution with time variable
% currently dual infeasible, there may be a bug?
%% parameters of system
% Tf = 120;
% tau = 30;

Tf = 30;
tau = 6;

k0 = 0.25;
ktau = 0.75;

taus = tau/Tf;

%% sample trajectory 
coeff = zeros(1, tau);
coeff(1) = k0;
coeff(end) = ktau;

A = [-coeff;
    eye(tau-1), zeros(tau-1, 1)];
C = [1 zeros(1, tau-1)];
sys = ss(A, [], C, [], 1);
xh = ones(1, tau);

[y, t, x] = initial(sys, xh, Tf-1);
eA = abs(eig(A));

%moment analysis
order = 6;
d = 2*order;
dv = genPowGlopti(2,d);
% m_traj = sum(y.^[0:d], 1)'; %moments
m_traj = monom_sum(t, y, dv, Tf);


%% Set up variables and measures 
mset clear
mpol('t', 'x0', 'x1');   mu = meas(t, x0, x1);
mpol('tT', 'xT');        muT = meas(tT, xT);
mpol('tnn', 'xnn');      nun = meas(tnn, xnn);
mpol('tnz', 'xnz');      nuz = meas(tnz, xnz);
mpol('tnp', 'xnp');      nup = meas(tnp, xnp);

%support constraints
supp_con = [tT == 1; 
            tnz*((1 - taus) - tnz) >= 0; 
            (tnp - (1 - taus)) * (1 - tnp) >= 0;
            t*(1 - t) >= 0;     
            x0^2 <= 1;
            x1^2 <= 1;
            xT^2 <= 1;
            xnz^2 <= 1;
            xnp^2 <= 1];
    
%% reference variables and measures 
y0 = prod([0, xh(1)].^dv, 2);
yh = monom_sum((-tau:-1), xh, dv, Tf);

%% Affine relations
%Liouville Equation
yT = mmon([tT, xT], d);
v0 = mmon([t, x0], d);
f  = -k0*x0 -ktau*x1;
push = subs(v0, [t, x0], [t + 1/Tf, f]);
Ay = mom(push - v0);
Liou = Ay + y0 - mom(yT);

%(t, x0) marginal
phi0 = mom(v0) - mom(mmon([tnz; xnz], d)) - mom(mmon([tnp; xnp], d));

%(t, x1) marginal (with shifting)
%test functions on components (t+tau)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, T-tau] component
mnz = mmon([tnz; xnz], d);
mnz_shift = subs(mnz, tnz, tnz + taus);  

%[-tau, 0] component]
mnn = mmon([tnn; xnn], d);
mnn_shift = subs(mnn, tnn, tnn + taus);

phi1 = mom(v1) - mom(mnz_shift) - mom(mnn_shift);

mom_con = [-Liou == 0;        %Liouville
           -phi0 == 0;        %x(t)
           -phi1 == 0;        %x(t - tau)
           mom(mnn) == yh];   %history given

objective = max(mass(mnn)); 
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
P = msdp(objective, ...
    mom_con, supp_con);

% solve LMIP moment problem
[status,obj_rec, m,dual_rec]= msol(P);

%% analyze moments
if status > 0
    m_mom  = double(mom(v0));
    m_comp = [m_traj m_mom];
    norm_diff = norm(m_traj-m_mom);
end


%% helper functions
function em = monom_sum(t, x, dv, T)
    %sum over t^alpha x(t)^beta in span [0, T]
    %em: empirical moments
    %there is probably a more efficient way to implement this
    if nargin < 4
        T = 1;
    end
    n_monom = size(dv, 1);
    em = zeros(n_monom, 1);
    for i = 1:n_monom        
        em(i) = sum(((t/T).^dv(i, 1)) .* (x.^dv(i, 2)))/T;
    end
end