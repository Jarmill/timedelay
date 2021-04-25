% Discrete system with long discrete delays
% weak solution without time variable (hopefully)
% does not work. will try with time variable
%% parameters of system
% Ts = 1;
% Tf = 120;
% z= tf([1 0], [1], Ts);
% 
% tau = 30;
% 
% k0 = 0.25;
% ktau = 0.75;

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

[y, t, x] = initial(sys, xh, Tf);
eA = abs(eig(A));

%moment analysis
order = 9;
d = 2*order;
m_traj = sum(y.^[0:d], 1)'; %moments
M_traj = hankel(m_traj(1:(order+1)), m_traj((order+1):end)); %moment matrix


%% Set up variables and measures 
mset clear
mpol('x0', 'x1'); mu = meas(x0, x1);
mpol('xT');       muT = meas(xT);
mpol('xnn');      nun = meas(xnn);
mpol('xnz');      nuz = meas(xnz);
mpol('xnp');      nup = meas(xnp);

%support constraints
supp_con = [x0^2 <= 1;
            x1^2 <= 1;
            xT^2 <= 1;
            xnz^2 <= 1;
            xnp^2 <= 1];
    
%% reference variables and measures 
y0 = (xh(1).^[0:d])';
yh = sum(xh(2:end)'.^[0:d], 1)';

%% Affine relations
%Liouville Equation
yT = mmon(xT, d);
v0 = mmon(x0, d);
f  = -k0*x0 -ktau*x1;
push = subs(v0, x0, -f);
Ay = mom(push - v0);
Liou = Ay + y0 - mom(yT);

%(t, x0) marginal
phi0 = mom(v0) - mom(mmon(xnz, d)) - mom(mmon(xnp, d));

%(t, x1) marginal
%no time, so no shifting needed
v1 = mmon(x1, d);
phi1 = mom(v1) - mom(mmon(xnz, d)) - mom(mmon(xnn, d));
mnn = mmon(xnn, d);

mom_con = [mass(mu) <= Tf;  %mass of joint occupation is total time
           -Liou == 0;        %Liouville
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
m_mom  = double(mom(v0));
m_comp = [m_traj m_mom]
norm_diff = norm(m_traj-m_mom)