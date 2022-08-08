
%dual formulation for peak estimation of dde

%% setup

%variables
t = sdpvar(1,1);
x0 = sdpvar(2,1);
x1 = sdpvar(2,1);

%parameters and dynamics 
Tmax = 5;
tau = 0.5;
ts = tau/Tmax;

f = [x0(2); -x1(1) + (1/3).*x0(1).^3- x0(2)];

%initial set
C0 = [1.5; 0];
R0 = 0.4;

%unsafe set
Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5; 
w_c = -[1; 1];
c1f = Ru^2 - (x0(1) - Cu(1)).^2 - (x0(2) - Cu(2)).^2;
% w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x0(1) - Cu(1)) + w_c(2) * (x0(2) - Cu(2)); 

cf = [c1f; c2f];

order = 3;
d=2*order;

%objective
p = x0(2);

%functions
[v, cv] = polynomial([t; x0], d);
[phi, cphi] = polynomial([t; x0], d);
[xi, cxi] = polynomial(t, d); 
gamma = sdpvar(1,1);


leb_mom = -(ts).^(1:(d+1)) ./ (1:(d+1));

Lv = jacobian(v, t) + Tmax * jacobian(v,x0)*f;

v0 = replace(v, t, 0);


%% support sets
% X0 = 
% Xu = 
box = [-1.25, 2.5; -1.25, 1.5];
X0 = [(x0 - box(:, 1)).*(box(:, 2) - x0)];
X1 = replace(X0, x0, x1);
tp = t*(1-t);
tn = (t + ts)*(-t);

t0 = t*((1-ts) - t);
t1 = (t - (1-ts))* (1-t);

H0 = R0^2 - sum((x0-C0).^2);
H = [H0; tn];

phi1 = replace(phi, x0, x1);
phishift = replace(phi, t, t+ts);


%% psatz constraints

[p0, cons0, coeff0] = constraint_psatz(v0 - gamma, H0, [x0], d);

[pc, consc, coeffc] = constraint_psatz(p - v, [tp; X0], [t; x0], d);

[ph, consh, coeffh] = constraint_psatz(phishift - xi, H0, [t; x0], d);

[plie0, conslie0, coefflie0] = constraint_psatz(Lv + phishift - phi1, [t0; X0; X1], [t; x0; x1], d+2);

[plie1, conslie1, coefflie1] = constraint_psatz(Lv + 0 - phi1, [t1; X0; X1], [t; x0; x1], d+2);

[pfree, consfree, coefffree] = constraint_psatz(-phi, [tp; X0], [t; x0], d);


%% accumulate constraints
objective = -(gamma + leb_mom*cxi);

coeff = [cv; cphi; cxi; gamma;
        coeff0; coeffc; coeffh; coefflie0;coefflie1; coefffree];
    
cons = [cons0; consc; consh; conslie0; conslie1; consfree];

%% solve
opts = sdpsettings('solver','mosek');

opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

